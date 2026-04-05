NUTS_method <- function(model,
                        sampling, warmup, delta = 0.8,
                        max_treedepth = 10, chain,
                        update_progress = NULL,
                        laplace) {

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
  da_iter <- 1 # DA専用のカウンタ

  # --- Windowed Adaptation 設定 ---
  init_buffer <- 75
  term_buffer <- 50
  base_window <- 25

  # ウォームアップが短すぎる場合の調整
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
      s <- out$s *
        safe_uturn(q_plus, q_minus, p_minus, M_inv) *
        safe_uturn(q_plus, q_minus, p_plus, M_inv)
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

    # --- ウォームアップ時の適応 ---
    if (i <= warmup) {
      # 1. Dual Averaging によるステップサイズ (eps) の更新
      Hbar <- (1 - 1/(da_iter+t0))*Hbar + (1/(da_iter+t0))*(delta - accept_stat)
      log_eps <- mu_DA - sqrt(da_iter)/gamma_DA * Hbar
      log_eps_bar <- da_iter^-kappa*log_eps + (1 - da_iter^-kappa)*log_eps_bar
      eps <- exp(log_eps)
      eps_bar <- exp(log_eps_bar)
      da_iter <- da_iter + 1

      # 2. Windowed Adaptation による質量行列 (M_inv) の更新
      if (base_window > 0 && i >= init_buffer && i <= warmup - term_buffer) {
        # Welfordのオンライン分散更新
        w_n <- w_n + 1
        delta_q <- q_new - w_mean
        w_mean <- w_mean + delta_q / w_n
        w_var <- w_var + delta_q * (q_new - w_mean)

        # ウィンドウの終端に達した場合の処理
        if (i == next_window || i == warmup - term_buffer) {
          # 分散の計算と極端な値のクリッピング
          new_M_inv <- w_var / (w_n - 1)
          new_M_inv[is.na(new_M_inv) | new_M_inv < 1e-8] <- 1e-8
          M_inv <- new_M_inv

          # Welfordの状態をリセット
          w_mean <- rep(0, P_fixed)
          w_var  <- rep(0, P_fixed)
          w_n    <- 0

          # 次のウィンドウのサイズを倍増させて設定
          window_size <- window_size * 2
          next_window <- i + window_size
          if (next_window > warmup - term_buffer) {
            next_window <- warmup - term_buffer
          }

          # 空間のスケールが変わったため、epsを再探索し Dual Averaging をリセット
          eps <- FindReasonableEpsilon(q_new, M_inv)
          mu_DA <- log(10 * eps)
          Hbar <- 0
          eps_bar <- 1
          log_eps_bar <- log(eps_bar)
          da_iter <- 1
        }
      }
    } else {
      # サンプリング期間は適応済みの eps_bar を使用
      eps <- eps_bar
    }

    if (i %% 100 == 0) {
      if (!is.null(update_progress)) update_progress()
      else cat(paste0("chain ", chain, ": iter ", i,
                      ifelse(i <= warmup, " warmup", " sampling"), "\n"))
    }
  }

  return(list(
    para_fixed = para_fixed,
    para_full  = para_full,
    lp         = lp,
    accept     = accept,
    treedepth  = treedepth_record,
    eps        = eps
  ))
}
