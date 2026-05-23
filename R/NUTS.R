NUTS_method <- function(model,
                        sampling, warmup, delta = 0.8,
                        max_treedepth = 10, chain,
                        update_progress = NULL,
                        laplace, save_info = NULL,
                        nuts_variant = c("multinomial", "slice"),
                        metric = c("diag", "dense"),
                        metric_init = c("identity", "hessian"),
                        metric_adaptation = c("stan_window", "cumulative", "window"),
                        metric_regularization = TRUE,
                        metric_shrinkage = 5,
                        metric_min = 1e-6,
                        metric_max = 1e6) {

  nuts_variant <- match.arg(nuts_variant)
  metric <- match.arg(metric)
  metric_init <- match.arg(metric_init)
  metric_adaptation <- match.arg(metric_adaptation)
  metric_regularization <- isTRUE(metric_regularization)
  if (!is.numeric(metric_shrinkage) || length(metric_shrinkage) != 1L || metric_shrinkage < 0) {
    stop("'metric_shrinkage' must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(metric_min) || length(metric_min) != 1L || metric_min <= 0) {
    stop("'metric_min' must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(metric_max) || length(metric_max) != 1L || metric_max <= metric_min) {
    stop("'metric_max' must be larger than 'metric_min'.", call. = FALSE)
  }

  iter <- sampling + warmup
  P_all <- length(model$env$last.par)

  nuts_core <- create_NUTS_core(model)
  fn_NUTS <- nuts_core$fn_NUTS
  gr_NUTS <- nuts_core$gr_NUTS
  calc_H <- nuts_core$calc_H
  safe_uturn <- nuts_core$safe_uturn
  BuildTreeSlice <- nuts_core$BuildTreeSlice
  RunMultinomialTransition <- nuts_core$RunMultinomialTransition
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
  n_leapfrog_record <- numeric(iter)
  divergent_record <- logical(iter)
  energy_record <- numeric(iter)
  eps_record <- numeric(iter)
  metric_condition_record <- numeric(iter)
  warmup_phase_record <- rep(NA_character_, iter)
  warmup_window_record <- rep(NA_integer_, iter)
  metric_updated_record <- rep(FALSE, iter)
  accept <- numeric(iter)
  accept[1] <- 1

  M_inv <- if (identical(metric_init, "hessian")) {
    nuts_core$hessian_metric(q_fixed_init, metric, metric_min, metric_max)
  } else {
    NULL
  }
  if (is.null(M_inv)) M_inv <- if (identical(metric, "dense")) diag(P_fixed) else rep(1, P_fixed)
  eps <- FindReasonableEpsilon(q_fixed_init, M_inv)

  # Dual Averaging parameters
  mu_DA <- log(10 * eps)
  Hbar <- 0; gamma_DA <- 0.05; t0 <- 10; kappa <- 0.75
  eps_bar <- 1; log_eps_bar <- log(eps_bar)
  da_iter <- 1

  # --- Windowed Adaptation setting ---
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
  window_id <- 0L

  w_mean <- rep(0, P_fixed)
  w_var  <- if (identical(metric, "dense")) matrix(0, P_fixed, P_fixed) else rep(0, P_fixed)
  w_n    <- 0
  c_mean <- rep(0, P_fixed)
  c_var  <- if (identical(metric, "dense")) matrix(0, P_fixed, P_fixed) else rep(0, P_fixed)
  c_n    <- 0

  msg_start <- paste0("chain ", chain, " started...")
  if (!is.null(update_progress)) {
    update_progress(msg = msg_start, amt = 1)
  } else {
    cat(msg_start, "\n")
  }

  for (i in 2:iter) {
    if (!is.null(model$rtmb_pd_error_state)) {
      model$rtmb_pd_error_state$post_warmup <- i > warmup
    }
    q_old <- para_fixed[i-1, ]
    p_old <- nuts_core$rmomentum(M_inv)
    gr_old <- gr_NUTS(q_old)
    H_old <- calc_H(q_old, p_old, M_inv)
    energy_record[i] <- -H_old
    eps_record[i] <- eps
    if (i <= warmup) {
      warmup_phase_record[i] <- if (base_window > 0 && i >= init_buffer && i <= warmup - term_buffer) {
        "slow"
      } else if (i < init_buffer) {
        "initial"
      } else {
        "terminal"
      }
      warmup_window_record[i] <- if (identical(warmup_phase_record[i], "slow")) window_id + 1L else NA_integer_
    }

    if (identical(nuts_variant, "multinomial")) {
      transition <- RunMultinomialTransition(q_old, p_old, gr_old, H_old, eps, max_treedepth, M_inv)
      q_new <- transition$q
      j <- transition$treedepth
      n_leapfrog <- transition$n_leapfrog
      divergent <- transition$divergent
      accept_stat <- transition$accept_stat
    } else {
      log_u <- log(runif(1)) + H_old
      q_minus <- q_old; q_plus <- q_old
      p_minus <- p_old; p_plus <- p_old
      gr_minus <- gr_old; gr_plus <- gr_old

      j <- 0; q_new <- q_old; n <- 1; s <- 1
      n_leapfrog <- 0L
      divergent <- FALSE
      alpha_sum <- 0; n_alpha_sum <- 0

      while (s == 1 && j < max_treedepth) {
        v <- ifelse(runif(1) > 0.5, -1, 1)
        if (v == -1) {
          out <- BuildTreeSlice(p_minus, q_minus, gr_minus, log_u, v, j, eps, H_old, M_inv)
          p_minus <- out$p_minus; q_minus <- out$q_minus; gr_minus <- out$gr_minus
        } else {
          out <- BuildTreeSlice(p_plus, q_plus, gr_plus, log_u, v, j, eps, H_old, M_inv)
          p_plus <- out$p_plus; q_plus <- out$q_plus; gr_plus <- out$gr_plus
        }
        if (out$s == 1) {
          if (runif(1) < out$n / n) q_new <- out$q
        }
        n <- n + out$n
        s <- out$s * safe_uturn(q_plus, q_minus, p_minus, M_inv) * safe_uturn(q_plus, q_minus, p_plus, M_inv)
        n_leapfrog <- n_leapfrog + out$n_leapfrog
        divergent <- divergent || isTRUE(out$divergent)
        alpha_sum <- alpha_sum + out$alpha; n_alpha_sum <- n_alpha_sum + out$n_alpha
        j <- j + 1
      }

      accept_stat <- alpha_sum / n_alpha_sum
      if (is.na(accept_stat)) accept_stat <- 0
    }

    para_fixed[i, ] <- q_new
    lp[i] <- fn_NUTS(q_new)
    para_full[i, ] <- model$env$last.par

    treedepth_record[i] <- j
    n_leapfrog_record[i] <- n_leapfrog
    divergent_record[i] <- divergent
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
        if (identical(metric, "dense")) {
          w_var <- w_var + tcrossprod(delta_q, q_new - w_mean)
        } else {
          w_var <- w_var + delta_q * (q_new - w_mean)
        }

        c_n <- c_n + 1
        delta_q_c <- q_new - c_mean
        c_mean <- c_mean + delta_q_c / c_n
        if (identical(metric, "dense")) {
          c_var <- c_var + tcrossprod(delta_q_c, q_new - c_mean)
        } else {
          c_var <- c_var + delta_q_c * (q_new - c_mean)
        }

        if (i == next_window || i == warmup - term_buffer) {
          metric_stats <- if (identical(metric_adaptation, "cumulative")) {
            list(var = c_var, n = c_n)
          } else {
            list(var = w_var, n = w_n)
          }
          min_metric_samples <- if (identical(metric, "dense")) max(5L, P_fixed + 1L) else 3L
          if (metric_stats$n < min_metric_samples) {
            metric_stats <- NULL
          }
          if (!is.null(metric_stats)) {
          M_inv <- nuts_core$regularize_metric(
            metric_stats$var, metric_stats$n,
            metric = metric,
            metric_regularization = metric_regularization,
            metric_shrinkage = metric_shrinkage,
            metric_min = metric_min,
            metric_max = metric_max
          )
          metric_updated_record[i] <- TRUE
          }

          w_mean <- rep(0, P_fixed)
          w_var <- if (identical(metric, "dense")) matrix(0, P_fixed, P_fixed) else rep(0, P_fixed)
          w_n <- 0
          window_id <- window_id + 1L
          window_size <- window_size * 2
          next_window <- i + window_size
          if (next_window > warmup - term_buffer) next_window <- warmup - term_buffer

          if (!is.null(metric_stats)) {
            eps <- FindReasonableEpsilon(q_new, M_inv)
            mu_DA <- log(10 * eps); Hbar <- 0; eps_bar <- 1; log_eps_bar <- log(eps_bar); da_iter <- 1
          }
        }
      }
      metric_condition_record[i] <- nuts_core$metric_condition(M_inv)
    } else {
      eps <- eps_bar
      metric_condition_record[i] <- nuts_core$metric_condition(M_inv)
    }
    if (i %% 200 == 0) {
      msg <- paste0("chain ", chain, ": iter ", i, ifelse(i <= warmup, " warmup", " sampling"))
      if (!is.null(update_progress)) {
        update_progress(amt = 1, msg = msg)
      } else {
        cat(msg, "\n")
      }
    }
  }

  if (exists("para_fixed")) para_fixed <- Re(para_fixed)
  if (exists("para_full") && !is.null(para_full)) para_full <- Re(para_full)
  warmup_diagnostics <- if (warmup >= 2L) {
    idx <- 2:warmup
    data.frame(
      iteration = idx,
      phase = warmup_phase_record[idx],
      window = warmup_window_record[idx],
      accept = accept[idx],
      treedepth = treedepth_record[idx],
      n_leapfrog = n_leapfrog_record[idx],
      divergent = divergent_record[idx],
      energy = Re(energy_record[idx]),
      eps = eps_record[idx],
      metric_condition = metric_condition_record[idx],
      metric_updated = metric_updated_record[idx]
    )
  } else {
    data.frame()
  }

  return(list(
    para_fixed = para_fixed, para_full = para_full,
    lp = Re(lp), accept = accept, treedepth = treedepth_record, eps = eps,
    n_leapfrog = n_leapfrog_record, divergent = divergent_record,
    energy = Re(energy_record), metric = M_inv, metric_type = metric,
    metric_init = metric_init, metric_adaptation = metric_adaptation,
    nuts_variant = nuts_variant,
    warmup_diagnostics = warmup_diagnostics,
    pd_error_count = if (!is.null(model$rtmb_pd_error_state)) model$rtmb_pd_error_state$count else 0L
  ))
}

create_NUTS_core <- function(ad_obj) {
  fn_NUTS <- function(q) -ad_obj$fn(q)
  gr_NUTS <- function(q) as.numeric(-ad_obj$gr(q))

  metric_mul <- function(M_inv, p) {
    p <- as.numeric(p)
    if (is.matrix(M_inv)) drop(M_inv %*% p) else M_inv * p
  }

  metric_quad <- function(M_inv, p) {
    p <- as.numeric(p)
    sum(p * metric_mul(M_inv, p))
  }

  metric_condition <- function(M_inv) {
    if (!is.matrix(M_inv)) {
      vals <- M_inv[is.finite(M_inv) & M_inv > 0]
      if (length(vals) == 0L) return(NA_real_)
      return(max(vals) / min(vals))
    }
    vals <- tryCatch(eigen(M_inv, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
    vals <- vals[is.finite(vals) & vals > 0]
    if (length(vals) == 0L) return(NA_real_)
    max(vals) / min(vals)
  }

  rmomentum <- function(M_inv) {
    z <- rnorm(if (is.matrix(M_inv)) nrow(M_inv) else length(M_inv))
    if (!is.matrix(M_inv)) return(z / sqrt(M_inv))
    R <- tryCatch(chol(M_inv), error = function(e) NULL)
    if (is.null(R)) return(z)
    as.numeric(backsolve(R, z))
  }

  regularize_metric <- function(w_var, w_n, metric,
                                metric_regularization,
                                metric_shrinkage,
                                metric_min,
                                metric_max) {
    if (!identical(metric, "dense")) {
      new_var <- w_var / max(1, w_n - 1)
      new_var[!is.finite(new_var)] <- 1
      new_var <- pmax(new_var, 1e-8)
      if (metric_regularization && metric_shrinkage > 0) {
        sample_weight <- w_n / (w_n + metric_shrinkage)
        jitter_weight <- metric_shrinkage / (w_n + metric_shrinkage)
        new_var <- sample_weight * new_var + 1e-3 * jitter_weight
      }
      return(pmin(pmax(new_var, metric_min), metric_max))
    }

    P <- nrow(w_var)
    cov_mat <- w_var / max(1, w_n - 1)
    cov_mat[!is.finite(cov_mat)] <- 0
    cov_mat <- (cov_mat + t(cov_mat)) / 2
    diag(cov_mat) <- pmax(diag(cov_mat), 1e-8)
    if (metric_regularization && metric_shrinkage > 0) {
      sample_weight <- w_n / (w_n + metric_shrinkage)
      jitter_weight <- metric_shrinkage / (w_n + metric_shrinkage)
      cov_mat <- sample_weight * cov_mat
      diag(cov_mat) <- diag(cov_mat) + 1e-3 * jitter_weight
    }

    eig <- tryCatch(eigen(cov_mat, symmetric = TRUE), error = function(e) NULL)
    if (is.null(eig) || any(!is.finite(eig$values))) {
      diag_var <- pmin(pmax(diag(cov_mat), metric_min), metric_max)
      return(diag(diag_var, P))
    }
    vals <- pmin(pmax(eig$values, metric_min), metric_max)
    cov_mat <- eig$vectors %*% (vals * t(eig$vectors))
    cov_mat <- (cov_mat + t(cov_mat)) / 2
    chol_ok <- !inherits(try(chol(cov_mat), silent = TRUE), "try-error")
    if (!chol_ok) {
      diag_var <- pmin(pmax(diag(cov_mat), metric_min), metric_max)
      cov_mat <- diag(diag_var, P)
    }
    cov_mat
  }

  hessian_metric <- function(q, metric, metric_min, metric_max) {
    H <- tryCatch({
      if (!is.null(ad_obj$he)) ad_obj$he(q) else stats::optimHess(q, ad_obj$fn, ad_obj$gr)
    }, error = function(e) NULL)
    if (is.null(H)) return(NULL)
    H <- as.matrix(H)
    if (nrow(H) != length(q) || ncol(H) != length(q) || any(!is.finite(H))) return(NULL)
    H <- (H + t(H)) / 2

    if (!identical(metric, "dense")) {
      d <- diag(H)
      if (any(!is.finite(d) | d <= 0)) return(NULL)
      return(pmin(pmax(1 / d, metric_min), metric_max))
    }

    eig <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
    if (is.null(eig) || any(!is.finite(eig$values)) || any(eig$values <= 0)) return(NULL)
    vals <- pmin(pmax(1 / eig$values, metric_min), metric_max)
    cov_mat <- eig$vectors %*% (vals * t(eig$vectors))
    cov_mat <- (cov_mat + t(cov_mat)) / 2
    if (inherits(try(chol(cov_mat), silent = TRUE), "try-error")) return(NULL)
    cov_mat
  }

  calc_H <- function(q, p, M_inv) {
    fn_NUTS(q) - metric_quad(M_inv, p) / 2
  }

  safe_uturn <- function(q_plus, q_minus, p, M_inv) {
    dq <- q_plus - q_minus
    dot_prod <- sum(dq * metric_mul(M_inv, p))
    if (!is.finite(dot_prod)) return(0L)
    if (dot_prod >= 0) 1L else 0L
  }

  sharp_momentum <- function(p, M_inv) {
    metric_mul(M_inv, p)
  }

  safe_multinomial_uturn <- function(rho, p_sharp) {
    dot_prod <- sum(as.numeric(rho) * as.numeric(p_sharp))
    if (!is.finite(dot_prod)) return(0L)
    if (dot_prod > 0) 1L else 0L
  }

  multinomial_criterion <- function(p_sharp_left, p_sharp_right, rho) {
    isTRUE(safe_multinomial_uturn(rho, p_sharp_left) == 1L) &&
      isTRUE(safe_multinomial_uturn(rho, p_sharp_right) == 1L)
  }

  log_sum_exp2 <- function(a, b) {
    if (!is.finite(a)) return(b)
    if (!is.finite(b)) return(a)
    m <- max(a, b)
    m + log(exp(a - m) + exp(b - m))
  }

  BuildTreeSlice <- function(p, q, gr, log_u, v, j, eps, H0, M_inv) {
    if (j == 0L) {
      eps_v <- v * eps
      p_prime <- p + (eps_v / 2) * gr
      q_prime <- q + eps_v * metric_mul(M_inv, p_prime)
      gr_prime <- gr_NUTS(q_prime)
      p_prime <- p_prime + (eps_v / 2) * gr_prime

      H_prime <- calc_H(q_prime, p_prime, M_inv)
      divergent_prime <- !is.finite(H_prime) || isTRUE(H0 - H_prime > 1000)

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
        n_alpha = 1L,
        n_leapfrog = 1L,
        divergent = divergent_prime
      ))
    }

    out1 <- BuildTreeSlice(p, q, gr, log_u, v, j - 1L, eps, H0, M_inv)
    if (isTRUE(out1$s == 0L)) return(out1)

    if (v == -1L) {
      out2 <- BuildTreeSlice(out1$p_minus, out1$q_minus, out1$gr_minus, log_u, v, j - 1L, eps, H0, M_inv)
      out1$p_minus <- out2$p_minus; out1$q_minus <- out2$q_minus; out1$gr_minus <- out2$gr_minus
    } else {
      out2 <- BuildTreeSlice(out1$p_plus, out1$q_plus, out1$gr_plus, log_u, v, j - 1L, eps, H0, M_inv)
      out1$p_plus <- out2$p_plus; out1$q_plus <- out2$q_plus; out1$gr_plus <- out2$gr_plus
    }

    total_n <- out1$n + out2$n
    if (isTRUE(out2$s == 1L) && total_n > 0L) {
      prob <- out2$n / total_n
      if (isTRUE(runif(1) < prob)) out1$q <- out2$q
    }

    s_uturn <- isTRUE(safe_uturn(out1$q_plus, out1$q_minus, out1$p_minus, M_inv) == 1L) &&
      isTRUE(safe_uturn(out1$q_plus, out1$q_minus, out1$p_plus, M_inv) == 1L)

    out1$n <- total_n

    out1$s <- if (isTRUE(out1$s == 1L) && isTRUE(out2$s == 1L) && s_uturn) 1L else 0L

    out1$alpha <- out1$alpha + out2$alpha
    out1$n_alpha <- out1$n_alpha + out2$n_alpha
    out1$n_leapfrog <- out1$n_leapfrog + out2$n_leapfrog
    out1$divergent <- out1$divergent || out2$divergent

    return(out1)
  }

  BuildTreeMultinomial <- function(p, q, gr, log_u, v, j, eps, H0, M_inv) {
    if (j == 0L) {
      eps_v <- v * eps
      p_prime <- p + (eps_v / 2) * gr
      q_prime <- q + eps_v * metric_mul(M_inv, p_prime)
      gr_prime <- gr_NUTS(q_prime)
      p_prime <- p_prime + (eps_v / 2) * gr_prime

      H_prime <- calc_H(q_prime, p_prime, M_inv)
      divergent_prime <- !is.finite(H_prime) || isTRUE(H0 - H_prime > 1000)

      log_weight_prime <- if (is.finite(H_prime)) H_prime - H0 else -Inf
      s_prime <- if (isTRUE(!divergent_prime)) 1L else 0L
      alpha_prime <- if (!is.finite(log_weight_prime)) {
        0
      } else if (log_weight_prime > 0) {
        1
      } else {
        exp(log_weight_prime)
      }
      p_sharp_prime <- sharp_momentum(p_prime, M_inv)

      return(list(
        p_beg = p_prime, q_beg = q_prime, gr_beg = gr_prime,
        p_end = p_prime, q_end = q_prime, gr_end = gr_prime,
        p_sharp_beg = p_sharp_prime,
        p_sharp_end = p_sharp_prime,
        q = q_prime,
        log_weight = log_weight_prime,
        rho = p_prime,
        valid = isTRUE(s_prime == 1L),
        s = s_prime,
        alpha = alpha_prime,
        n_alpha = 1L,
        n_leapfrog = 1L,
        divergent = divergent_prime
      ))
    }

    out1 <- BuildTreeMultinomial(p, q, gr, log_u, v, j - 1L, eps, H0, M_inv)
    if (!isTRUE(out1$valid)) return(out1)

    out2 <- BuildTreeMultinomial(out1$p_end, out1$q_end, out1$gr_end, log_u, v, j - 1L, eps, H0, M_inv)

    alpha_total <- out1$alpha + out2$alpha
    n_alpha_total <- out1$n_alpha + out2$n_alpha
    n_leapfrog_total <- out1$n_leapfrog + out2$n_leapfrog
    divergent_total <- out1$divergent || out2$divergent

    if (!isTRUE(out2$valid)) {
      out1$alpha <- alpha_total
      out1$n_alpha <- n_alpha_total
      out1$n_leapfrog <- n_leapfrog_total
      out1$divergent <- divergent_total
      out1$valid <- FALSE
      out1$s <- 0L
      return(out1)
    }

    log_weight_total <- log_sum_exp2(out1$log_weight, out2$log_weight)
    q_propose <- out1$q
    if (out2$log_weight > log_weight_total) {
      q_propose <- out2$q
    } else if (is.finite(log_weight_total) &&
               isTRUE(runif(1) < exp(out2$log_weight - log_weight_total))) {
      q_propose <- out2$q
    }

    rho_total <- out1$rho + out2$rho
    s_uturn <- multinomial_criterion(out1$p_sharp_beg, out2$p_sharp_end, rho_total) &&
      multinomial_criterion(out1$p_sharp_beg, out2$p_sharp_beg, out1$rho + out2$p_beg) &&
      multinomial_criterion(out1$p_sharp_end, out2$p_sharp_end, out2$rho + out1$p_end)

    list(
      p_beg = out1$p_beg, q_beg = out1$q_beg, gr_beg = out1$gr_beg,
      p_end = out2$p_end, q_end = out2$q_end, gr_end = out2$gr_end,
      p_sharp_beg = out1$p_sharp_beg,
      p_sharp_end = out2$p_sharp_end,
      q = q_propose,
      log_weight = log_weight_total,
      rho = rho_total,
      valid = isTRUE(s_uturn),
      s = if (isTRUE(s_uturn)) 1L else 0L,
      alpha = alpha_total,
      n_alpha = n_alpha_total,
      n_leapfrog = n_leapfrog_total,
      divergent = divergent_total
    )
  }

  RunMultinomialTransition <- function(q_old, p_old, gr_old, H0, eps, max_treedepth, M_inv) {
    q_bck <- q_old; q_fwd <- q_old
    p_bck <- p_old; p_fwd <- p_old
    gr_bck <- gr_old; gr_fwd <- gr_old

    q_new <- q_old
    log_weight <- 0
    rho <- p_old
    p_sharp_old <- sharp_momentum(p_old, M_inv)
    p_fwd_fwd <- p_old
    p_fwd_bck <- p_old
    p_bck_fwd <- p_old
    p_bck_bck <- p_old
    p_sharp_fwd_fwd <- p_sharp_old
    p_sharp_fwd_bck <- p_sharp_old
    p_sharp_bck_fwd <- p_sharp_old
    p_sharp_bck_bck <- p_sharp_old

    j <- 0L
    n_leapfrog <- 0L
    divergent <- FALSE
    alpha_sum <- 0

    while (j < max_treedepth) {
      v <- ifelse(runif(1) > 0.5, -1L, 1L)
      if (v == -1L) {
        rho_fwd <- rho
        p_fwd_bck <- p_bck_bck
        p_sharp_fwd_bck <- p_sharp_bck_bck
        out <- BuildTreeMultinomial(p_bck, q_bck, gr_bck, NULL, v, j, eps, H0, M_inv)
        q_bck <- out$q_end; p_bck <- out$p_end; gr_bck <- out$gr_end
        p_bck_fwd <- out$p_beg
        p_bck_bck <- out$p_end
        p_sharp_bck_fwd <- out$p_sharp_beg
        p_sharp_bck_bck <- out$p_sharp_end
        rho_bck <- out$rho
      } else {
        rho_bck <- rho
        p_bck_fwd <- p_fwd_fwd
        p_sharp_bck_fwd <- p_sharp_fwd_fwd
        out <- BuildTreeMultinomial(p_fwd, q_fwd, gr_fwd, NULL, v, j, eps, H0, M_inv)
        q_fwd <- out$q_end; p_fwd <- out$p_end; gr_fwd <- out$gr_end
        p_fwd_bck <- out$p_beg
        p_fwd_fwd <- out$p_end
        p_sharp_fwd_bck <- out$p_sharp_beg
        p_sharp_fwd_fwd <- out$p_sharp_end
        rho_fwd <- out$rho
      }

      n_leapfrog <- n_leapfrog + out$n_leapfrog
      divergent <- divergent || isTRUE(out$divergent)
      alpha_sum <- alpha_sum + out$alpha

      if (!isTRUE(out$valid)) break

      j <- j + 1L

      if (out$log_weight > log_weight) {
        q_new <- out$q
      } else if (is.finite(out$log_weight) &&
                 isTRUE(runif(1) < exp(out$log_weight - log_weight))) {
        q_new <- out$q
      }
      log_weight <- log_sum_exp2(log_weight, out$log_weight)
      rho <- rho_bck + rho_fwd

      persist <- multinomial_criterion(p_sharp_bck_bck, p_sharp_fwd_fwd, rho)
      persist <- persist &&
        multinomial_criterion(p_sharp_bck_bck, p_sharp_fwd_bck, rho_bck + p_fwd_bck)
      persist <- persist &&
        multinomial_criterion(p_sharp_bck_fwd, p_sharp_fwd_fwd, rho_fwd + p_bck_fwd)
      if (!isTRUE(persist)) break
    }

    accept_stat <- alpha_sum / n_leapfrog
    if (is.na(accept_stat)) accept_stat <- 0

    list(
      q = q_new,
      treedepth = j,
      n_leapfrog = n_leapfrog,
      divergent = divergent,
      accept_stat = accept_stat
    )
  }

  FindReasonableEpsilon <- function(q, M_inv) {
    eps_try <- 1
    p <- rmomentum(M_inv)
    gr <- gr_NUTS(q)
    H_old <- calc_H(q, p, M_inv)

    p_new <- p + (eps_try / 2) * gr
    q_new <- q + eps_try * metric_mul(M_inv, p_new)
    p_new <- p_new + (eps_try / 2) * gr_NUTS(q_new)
    H_new <- calc_H(q_new, p_new, M_inv)

    ratio <- if (!is.finite(H_new) || !is.finite(H_old)) 0 else exp(H_new - H_old)
    a <- if (isTRUE(ratio > 0.5)) 1 else -1

    while (isTRUE(ratio^a > 2^(-a))) {
      eps_try <- (2^a) * eps_try
      p_new <- p + (eps_try / 2) * gr
      q_new <- q + eps_try * metric_mul(M_inv, p_new)
      p_new <- p_new + (eps_try / 2) * gr_NUTS(q_new)
      H_new <- calc_H(q_new, p_new, M_inv)
      ratio <- if (!is.finite(H_new) || !is.finite(H_old)) 0 else exp(H_new - H_old)
      if (eps_try > 1e7 || eps_try < 1e-7) break
    }
    eps_try
  }

  return(list(
    fn_NUTS = fn_NUTS, gr_NUTS = gr_NUTS, calc_H = calc_H,
    safe_uturn = safe_uturn,
    safe_multinomial_uturn = safe_multinomial_uturn,
    metric_condition = metric_condition,
    rmomentum = rmomentum,
    regularize_metric = regularize_metric,
    hessian_metric = hessian_metric,
    BuildTreeSlice = BuildTreeSlice,
    BuildTreeMultinomial = BuildTreeMultinomial,
    RunMultinomialTransition = RunMultinomialTransition,
    FindReasonableEpsilon = FindReasonableEpsilon,
    log_sum_exp2 = log_sum_exp2
  ))
}
