NUTS_method <- function(model,
                        sampling, warmup, delta = 0.8,
                        max_treedepth = 10, chain,
                        update_progress = NULL,
                        laplace, save_info = NULL) {

  iter <- sampling + warmup
  P_all <- length(model$env$last.par)

  nuts_core <- create_NUTS_core(model)
  fn_NUTS <- nuts_core$fn_NUTS
  gr_NUTS <- nuts_core$gr_NUTS
  calc_H <- nuts_core$calc_H
  safe_uturn <- nuts_core$safe_uturn
  BuildTree <- nuts_core$BuildTree
  FindReasonableEpsilon <- nuts_core$FindReasonableEpsilon

  q_fixed_init <- model$par
  P_fixed <- length(q_fixed_init)

  # --- CSV保存の準備 ---
  if (!is.null(save_info)) {
    backup_file <- file.path(save_info$dir, paste0(save_info$name, "-", chain, ".csv"))
    prog_file <- file.path(save_info$dir, paste0(save_info$name, "_progress_", chain, ".txt"))

    param_names <- names(q_fixed_init)
    if (is.null(param_names)) param_names <- paste0("V", 1:P_fixed)
    header <- c("iteration", "lp", "accept", "treedepth", "eps", param_names)
    write.table(t(header), file = backup_file, append = FALSE, sep = ",", col.names = FALSE, row.names = FALSE)
  }


  para_fixed <- array(NA, dim=c(iter, P_fixed))
  para_full  <- array(NA, dim=c(iter, P_all))

  para_fixed[1, ] <- q_fixed_init
  lp <- numeric(iter)
  lp[1] <- fn_NUTS(q_fixed_init)
  para_full[1, ] <- model$env$last.par

  treedepth_record <- numeric(iter)
  accept <- numeric(iter)
  accept[1] <- 1

  M_inv <- rep(1, P_fixed)
  eps <- FindReasonableEpsilon(q_fixed_init, M_inv)

  # Dual Averaging パラメータ
  mu_DA <- log(10 * eps)
  Hbar <- 0; gamma_DA <- 0.05; t0 <- 10; kappa <- 0.75
  eps_bar <- 1; log_eps_bar <- log(eps_bar)
  da_iter <- 1

  # --- Windowed Adaptation 設定 ---
  init_buffer <- 75
  term_buffer <- 50
  base_window <- 25

  if (warmup < 20) {
    init_buffer <- warmup
    term_buffer <- 0
    base_window <- 0
  } else if (warmup < init_buffer + term_buffer + base_window) {
    init_buffer <- floor(0.15 * warmup)
    term_buffer <- floor(0.10 * warmup)
    base_window <- warmup - init_buffer - term_buffer
  }

  next_window <- init_buffer + base_window
  window_size <- base_window

  # Welfordのアルゴリズム用状態変数
  w_mean <- rep(0, P_fixed)
  w_var  <- rep(0, P_fixed)
  w_n    <- 0

  for (i in 2:iter) {
    q_old <- para_fixed[i-1, ]
    p_old <- rnorm(P_fixed, mean = 0, sd = sqrt(1 / M_inv))
    gr_old <- gr_NUTS(q_old)
    H_old <- calc_H(q_old, p_old, M_inv)

    log_u <- log(runif(1)) + H_old
    q_minus <- q_old; q_plus <- q_old
    p_minus <- p_old; p_plus <- p_old
    gr_minus <- gr_old; gr_plus <- gr_old

    j <- 0; q_new <- q_old; n <- 1; s <- 1
    alpha_sum <- 0; n_alpha_sum <- 0

    while (s == 1 && j < max_treedepth) {
      v <- ifelse(runif(1) > 0.5, -1, 1)
      if (v == -1) {
        out <- BuildTree(p_minus, q_minus, gr_minus, log_u, v, j, eps, H_old, M_inv)
        p_minus <- out$p_minus; q_minus <- out$q_minus; gr_minus <- out$gr_minus
      } else {
        out <- BuildTree(p_plus, q_plus, gr_plus, log_u, v, j, eps, H_old, M_inv)
        p_plus <- out$p_plus; q_plus <- out$q_plus; gr_plus <- out$gr_plus
      }
      if (out$s == 1) {
        if (runif(1) < out$n / n) q_new <- out$q
      }
      n <- n + out$n
      s <- out$s * safe_uturn(q_plus, q_minus, p_minus, M_inv) * safe_uturn(q_plus, q_minus, p_plus, M_inv)
      alpha_sum <- alpha_sum + out$alpha; n_alpha_sum <- n_alpha_sum + out$n_alpha
      j <- j + 1
    }

    para_fixed[i, ] <- q_new
    lp[i] <- fn_NUTS(q_new)
    para_full[i, ] <- model$env$last.par

    treedepth_record[i] <- j
    accept_stat <- alpha_sum / n_alpha_sum
    if (is.na(accept_stat)) accept_stat <- 0
    accept[i] <- accept_stat

    if (i <= warmup) {
      Hbar <- (1 - 1/(da_iter+t0))*Hbar + (1/(da_iter+t0))*(delta - accept_stat)
      log_eps <- mu_DA - sqrt(da_iter)/gamma_DA * Hbar
      log_eps_bar <- da_iter^-kappa*log_eps + (1 - da_iter^-kappa)*log_eps_bar
      eps <- exp(log_eps)
      eps_bar <- exp(log_eps_bar)
      da_iter <- da_iter + 1

      if (base_window > 0 && i >= init_buffer && i <= warmup - term_buffer) {
        w_n <- w_n + 1
        delta_q <- q_new - w_mean
        w_mean <- w_mean + delta_q / w_n
        w_var <- w_var + delta_q * (q_new - w_mean)

        if (i == next_window || i == warmup - term_buffer) {
          new_M_inv <- w_var / (w_n - 1)
          new_M_inv[is.na(new_M_inv) | new_M_inv < 1e-8] <- 1e-8
          M_inv <- new_M_inv

          w_mean <- rep(0, P_fixed); w_var <- rep(0, P_fixed); w_n <- 0
          window_size <- window_size * 2
          next_window <- i + window_size
          if (next_window > warmup - term_buffer) next_window <- warmup - term_buffer

          eps <- FindReasonableEpsilon(q_new, M_inv)
          mu_DA <- log(10 * eps); Hbar <- 0; eps_bar <- 1; log_eps_bar <- log(eps_bar); da_iter <- 1
        }
      }
    } else {
      eps <- eps_bar
    }

    if (i %% 100 == 0) {
      msg <- paste0("chain ", chain, ": iter ", i, ifelse(i <= warmup, " warmup", " sampling"))
      if (!is.null(update_progress)) update_progress(msg)
      else cat(msg, "\n")

      # --- 100回ごとのファイル書き出し処理 ---
      if (!is.null(save_info)) {
        writeLines(as.character(i), con = prog_file) # 進捗を上書き

        start_idx <- i - 99
        if (start_idx == 1) start_idx <- 2
        len_idx <- length(start_idx:i)

        backup_data <- matrix(NA, nrow = len_idx, ncol = 5 + P_fixed)
        backup_data[, 1] <- start_idx:i
        backup_data[, 2] <- lp[start_idx:i]
        backup_data[, 3] <- accept[start_idx:i]
        backup_data[, 4] <- treedepth_record[start_idx:i]
        backup_data[, 5] <- eps
        backup_data[, 6:(5 + P_fixed)] <- para_fixed[start_idx:i, ]

        write.table(backup_data, file = backup_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
      }
    }
  }

  # --- ループ終了後の端数の保存処理 ---
  if (!is.null(save_info)) {
    remainder <- iter %% 100
    if (remainder > 0) {
      start_idx <- iter - remainder + 1
      len_idx <- length(start_idx:iter)
      backup_data <- matrix(NA, nrow = len_idx, ncol = 5 + P_fixed)
      backup_data[, 1] <- start_idx:iter
      backup_data[, 2] <- lp[start_idx:iter]
      backup_data[, 3] <- accept[start_idx:iter]
      backup_data[, 4] <- treedepth_record[start_idx:iter]
      backup_data[, 5] <- eps
      backup_data[, 6:(5 + P_fixed)] <- para_fixed[start_idx:iter, , drop = FALSE]

      write.table(backup_data, file = backup_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
      writeLines(as.character(iter), con = prog_file)
    }
  }

  return(list(
    para_fixed = para_fixed, para_full = para_full,
    lp = lp, accept = accept, treedepth = treedepth_record, eps = eps
  ))
}

create_NUTS_core <- function(ad_obj) {
  fn_NUTS <- function(q) -ad_obj$fn(q)
  gr_NUTS <- function(q) -ad_obj$gr(q)

  calc_H <- function(q, p, M_inv) {
    fn_NUTS(q) - sum((p^2) * M_inv) / 2
  }

  safe_uturn <- function(q_plus, q_minus, p, M_inv) {
    dq <- q_plus - q_minus
    dot_prod <- sum(dq * (p * M_inv))
    if (!is.finite(dot_prod)) return(0L)
    if (dot_prod >= 0) 1L else 0L
  }

  BuildTree <- function(p, q, gr, log_u, v, j, eps, H0, M_inv) {
    if (j == 0L) {
      eps_v <- v * eps
      p_prime <- p + (eps_v / 2) * gr
      q_prime <- q + eps_v * M_inv * p_prime
      gr_prime <- gr_NUTS(q_prime)
      p_prime <- p_prime + (eps_v / 2) * gr_prime

      H_prime <- calc_H(q_prime, p_prime, M_inv)

      if (!is.finite(H_prime)) {
        n_prime <- 0L; s_prime <- 0L; alpha_prime <- 0
      } else {
        n_prime <- if (isTRUE(log_u <= H_prime)) 1L else 0L
        s_prime <- if (isTRUE(log_u - 1000 < H_prime)) 1L else 0L
        alpha_raw <- exp(H_prime - H0)
        alpha_prime <- if (is.na(alpha_raw)) 0 else min(1, alpha_raw)
      }

      return(list(
        p_minus = p_prime, q_minus = q_prime, gr_minus = gr_prime,
        p_plus  = p_prime, q_plus  = q_prime, gr_plus  = gr_prime,
        q       = q_prime,
        n       = n_prime,
        s       = s_prime,
        alpha   = alpha_prime,
        n_alpha = 1L
      ))
    }

    out1 <- BuildTree(p, q, gr, log_u, v, j - 1L, eps, H0, M_inv)
    # isTRUEで確実に論理値にする
    if (isTRUE(out1$s == 0L)) return(out1)

    if (v == -1L) {
      out2 <- BuildTree(out1$p_minus, out1$q_minus, out1$gr_minus, log_u, v, j - 1L, eps, H0, M_inv)
      out1$p_minus <- out2$p_minus; out1$q_minus <- out2$q_minus; out1$gr_minus <- out2$gr_minus
    } else {
      out2 <- BuildTree(out1$p_plus, out1$q_plus, out1$gr_plus, log_u, v, j - 1L, eps, H0, M_inv)
      out1$p_plus <- out2$p_plus; out1$q_plus <- out2$q_plus; out1$gr_plus <- out2$gr_plus
    }

    total_n <- out1$n + out2$n
    if (total_n > 0L) {
      prob <- out2$n / total_n
      if (isTRUE(runif(1) < prob)) out1$q <- out2$q
    }

    dq <- out1$q_plus - out1$q_minus
    sum1 <- sum(dq * (out1$p_minus * M_inv))
    sum2 <- sum(dq * (out1$p_plus * M_inv))

    # 判定結果にNAが混入しないようisTRUEを使用
    s_uturn <- isTRUE(sum1 >= 0) && isTRUE(sum2 >= 0)

    out1$n <- total_n

    # 掛け算を使わず、論理積で確実に 0L か 1L にする
    out1$s <- if (isTRUE(out1$s == 1L) && isTRUE(out2$s == 1L) && s_uturn) 1L else 0L

    out1$alpha <- out1$alpha + out2$alpha
    out1$n_alpha <- out1$n_alpha + out2$n_alpha

    return(out1)
  }

  FindReasonableEpsilon <- function(q, M_inv) {
    eps_try <- 1
    p <- rnorm(length(q), mean = 0, sd = sqrt(1 / M_inv))
    gr <- gr_NUTS(q)
    H_old <- calc_H(q, p, M_inv)

    p_new <- p + (eps_try / 2) * gr
    q_new <- q + eps_try * M_inv * p_new
    p_new <- p_new + (eps_try / 2) * gr_NUTS(q_new)
    H_new <- calc_H(q_new, p_new, M_inv)

    ratio <- if (!is.finite(H_new) || !is.finite(H_old)) 0 else exp(H_new - H_old)
    a <- if (isTRUE(ratio > 0.5)) 1 else -1

    while (isTRUE(ratio^a > 2^(-a))) {
      eps_try <- (2^a) * eps_try
      p_new <- p + (eps_try / 2) * gr
      q_new <- q + eps_try * M_inv * p_new
      p_new <- p_new + (eps_try / 2) * gr_NUTS(q_new)
      H_new <- calc_H(q_new, p_new, M_inv)
      ratio <- if (!is.finite(H_new) || !is.finite(H_old)) 0 else exp(H_new - H_old)
      if (eps_try > 1e7 || eps_try < 1e-7) break
    }
    eps_try
  }

  return(list(
    fn_NUTS = fn_NUTS, gr_NUTS = gr_NUTS, calc_H = calc_H,
    safe_uturn = safe_uturn, BuildTree = BuildTree,
    FindReasonableEpsilon = FindReasonableEpsilon
  ))
}

