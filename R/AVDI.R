#' Automatic Differentiation Variational Inference (ADVI)
#'
#' @param model RTMB ad_obj
#' @param par_list Parameter definitions
#' @param pl_full Full parameter definitions
#' @param iter Maximum number of iterations
#' @param tol_rel_obj Relative tolerance for convergence
#' @param window_size Window size for median smoothing
#' @param num_samples Number of posterior draws
#' @param chains Number of chains to format the output
#' @param alpha Learning rate for Adam
#' @param laplace Logical; whether Laplace approximation is used
ADVI_method <- function(model, par_list, pl_full,
                        iter = 10000, min_iter = 1000, tol_rel_obj = 0.001,
                        window_size = 100, num_samples = 1000,
                        chains = 1, alpha = 0.01, laplace = FALSE,
                        print_freq = 500) {

  P_fixed <- length(model$par)

  # 変分パラメータの初期化 (平均場近似)
  # mu は初期値、omega (log_sd) は小さめの値(-2)からスタート
  mu <- model$par
  omega <- rep(-2, P_fixed)

  # Adamオプティマイザの初期化
  m_mu <- rep(0, P_fixed); v_mu <- rep(0, P_fixed)
  m_omega <- rep(0, P_fixed); v_omega <- rep(0, P_fixed)
  beta1 <- 0.9; beta2 <- 0.999; eps_adam <- 1e-8

  elbo_history <- numeric(iter)
  converged <- FALSE

  cat("Starting ADVI optimization with Adam...\n")

  for (t in 1:iter) {
    # 1. サンプリング (Reparameterization Trick)
    eta <- rnorm(P_fixed)
    zeta <- mu + exp(omega) * eta

    # 2. 勾配と対数尤度の評価
    # model$fn: -log p(zeta, Y), model$gr: -nabla log p(zeta, Y)
    fn_val <- tryCatch(model$fn(zeta), error = function(e) NA)
    gr_val <- tryCatch(model$gr(zeta), error = function(e) NA)

    if (any(is.na(fn_val)) || any(is.na(gr_val))) {
      warning(sprintf("Iteration %d: NaN in objective or gradient. Skipping update.", t))
      if (t > 1) elbo_history[t] <- elbo_history[t-1]
      next
    }

    # 3. 目的関数 (負のELBO) の勾配計算
    grad_mu <- gr_val
    grad_omega <- gr_val * exp(omega) * eta - 1

    # 4. Adamによるパラメータ更新 (mu)
    m_mu <- beta1 * m_mu + (1 - beta1) * grad_mu
    v_mu <- beta2 * v_mu + (1 - beta2) * (grad_mu^2)
    m_mu_hat <- m_mu / (1 - beta1^t)
    v_mu_hat <- v_mu / (1 - beta2^t)
    mu <- mu - alpha * m_mu_hat / (sqrt(v_mu_hat) + eps_adam)

    # 5. Adamによるパラメータ更新 (omega)
    m_omega <- beta1 * m_omega + (1 - beta1) * grad_omega
    v_omega <- beta2 * v_omega + (1 - beta2) * (grad_omega^2)
    m_omega_hat <- m_omega / (1 - beta1^t)
    v_omega_hat <- v_omega / (1 - beta2^t)
    omega <- omega - alpha * m_omega_hat / (sqrt(v_omega_hat) + eps_adam)

    # ELBOの計算と記録 (エントロピーの定数項は省略)
    elbo_history[t] <- -fn_val + sum(omega)

    if (print_freq > 0 && t %% print_freq == 0) {
      cat(sprintf("Iter %d: Approx ELBO = %.2f\n", t, elbo_history[t]))
    }

    check_start <- max(min_iter, 2 * window_size)

    # 6. 収束判定 (移動ウィンドウを用いたメディアンの変化量)
    if (t > check_start && t %% 10 == 0) {
      med_prev <- median(elbo_history[(t - 2 * window_size + 1):(t - window_size)])
      med_curr <- median(elbo_history[(t - window_size + 1):t])

      rel_change <- abs(med_curr - med_prev) / (abs(med_curr) + 1e-8)

      if (rel_change < tol_rel_obj) {
        cat(sprintf("Converged at iteration %d (Relative change: %.5f < %.5f)\n",
                    t, rel_change, tol_rel_obj))
        converged <- TRUE
        break
      }
    }
  }

  if (!converged) {
    warning("ADVI did not converge within the maximum number of iterations.")
  }

  # 実行したイテレーションまでの履歴に切り詰める（以降の 0 を削除）
  elbo_history <- elbo_history[1:t]

  # 最後のステップの平均値（直近 window_size 分）を計算
  calc_window <- min(t, window_size)
  elbo_final <- median(elbo_history[(t - calc_window + 1):t])
  # --- ここまで ---

  # --- 近似事後分布からのサンプリング ---
  cat("Generating posterior samples from variational distribution...\n")
  samples_per_chain <- ceiling(num_samples / chains)
  actual_num_samples <- samples_per_chain * chains

  random_flags <- sapply(par_list, function(x) isTRUE(x$random))
  if (laplace && any(random_flags)) {
    pl_fixed  <- parse_parameters(par_list[!random_flags])
    pl_random <- parse_parameters(par_list[random_flags])
  } else {
    pl_fixed  <- pl_full
    pl_random <- NULL
  }

  fixed_idx  <- which(pl_full$names %in% pl_fixed$names)
  random_idx <- if (!is.null(pl_random)) which(pl_full$names %in% pl_random$names) else integer(0)

  P_all_true <- length(pl_full$names)
  para_final <- array(NA, dim = c(actual_num_samples, P_all_true))
  lp_final <- numeric(actual_num_samples)

  sd_vec <- exp(omega)
  for (i in 1:actual_num_samples) {
    # 独立正規分布からのサンプル
    z <- rnorm(P_fixed)
    zeta_sample <- mu + sd_vec * z
    lp_final[i] <- -model$fn(zeta_sample)

    # 空間の復元と周辺化された変量効果の取得
    if (laplace && length(model$env$random) > 0) {
      para_list_res <- model$env$parList()
    } else {
      para_list_res <- model$env$parList(x = zeta_sample)
    }

    con_list <- to_constrained(para_list_res, par_list)
    para_final[i, ] <- unlist(con_list, use.names = FALSE)
  }

  # MCMC_Fit互換の配列に整形
  fit <- array(NA, dim = c(samples_per_chain, chains, length(pl_fixed$names) + 1))
  dimnames(fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = c("lp", pl_fixed$names))

  random_fit <- NULL
  if (!is.null(pl_random) && length(pl_random$names) > 0) {
    random_fit <- array(NA, dim = c(samples_per_chain, chains, length(pl_random$names)))
    dimnames(random_fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = pl_random$names)
  }

  idx <- 1
  for (c in 1:chains) {
    for (i in 1:samples_per_chain) {
      fit[i, c, 1] <- lp_final[idx]
      fit[i, c, -1] <- para_final[idx, fixed_idx]
      if (!is.null(random_fit)) {
        random_fit[i, c, ] <- para_final[idx, random_idx]
      }
      idx <- idx + 1
    }
  }

  return(list(
    fit          = fit,
    random_fit   = random_fit,
    elbo_history = elbo_history,
    elbo_final   = elbo_final
  ))
}
