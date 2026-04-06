#' Automatic Differentiation Variational Inference (ADVI)
#'
#' @param model An RTMB objective function object (`ad_obj`).
#' @param par_list A list defining the structure of parameters to be estimated.
#' @param pl_full A list defining the full structure of parameters including random effects.
#' @param min_iter Integer; minimum number of iterations to run before checking for convergence. Default is 1000.
#' @param iter Integer; maximum number of iterations for the optimization. Default is 10000.
#' @param tol_rel_obj Numeric; relative tolerance for the ELBO to determine convergence. Default is 0.01.
#' @param window_size Integer; size of the moving window to calculate the median ELBO for the convergence check. Default is 100.
#' @param num_samples Integer; number of posterior draws to generate after optimization. Default is 1000.
#' @param chains Integer; number of chains to structure the output array (for compatibility). Default is 1.
#' @param alpha Numeric; learning rate (step size) for the Adam optimizer. Default is 0.01.
#' @param laplace Logical; whether Laplace approximation is used. Default is FALSE.
#' @param print_freq Integer; frequency of printing progress to the console. Set to 0 to disable. Default is 500.
#' @param fullrank Logical; whether to use a full-rank multivariate normal distribution for the approximation. Default is FALSE.
#' @return A list containing `fit`, `random_fit`, `elbo_history`, and `elbo_final`.
ADVI_method <- function(model, par_list, pl_full,
                        iter = 10000, min_iter = 1000, tol_rel_obj = 0.01,
                        window_size = 100, num_samples = 1000,
                        chains = 1, alpha = 0.01, laplace = FALSE,
                        print_freq = 500,
                        fullrank = FALSE) {

  # --- 初期化 ---
  pl_fixed <- par_list
  pl_random <- pl_full[!(names(pl_full) %in% names(par_list))]
  P <- length(model$par)
  mu <- model$par

  # Adamのハイパーパラメータ
  beta1 <- 0.9
  beta2 <- 0.999
  epsilon <- 1e-8

  if (fullrank) {
    # Full-rank用のパラメータ: コレスキー分解の下三角行列 L
    L_diag <- rep(0, P) # 対角成分 (正の値を保証するため対数スケール)
    L_off <- rep(0, P * (P - 1) / 2) # 非対角成分

    m_mu <- rep(0, P); v_mu <- rep(0, P)
    m_diag <- rep(0, P); v_diag <- rep(0, P)
    m_off <- rep(0, length(L_off)); v_off <- rep(0, length(L_off))
  } else {
    # Mean-field用のパラメータ
    omega <- rep(0, P)
    m_mu <- rep(0, P); v_mu <- rep(0, P)
    m_omega <- rep(0, P); v_omega <- rep(0, P)
  }

  elbo_history <- numeric(iter)
  converged <- FALSE

  # --- 最適化ループ ---
  for (t in 1:iter) {
    eps <- rnorm(P)

    # 1. パラメータのサンプリング
    if (fullrank) {
      L <- matrix(0, P, P)
      L[lower.tri(L)] <- L_off
      diag(L) <- exp(L_diag)
      theta <- mu + as.vector(L %*% eps)
    } else {
      sigma <- exp(omega)
      theta <- mu + sigma * eps
    }

    # 2. 尤度と勾配の評価 (RTMBは負の対数尤度を返す前提)
    fn_val <- -model$fn(theta)
    gr_val <- -model$gr(theta) # 勾配 (d log_p / d theta)

    # 3. ELBOと勾配の計算
    if (fullrank) {
      elbo_history[t] <- fn_val + sum(L_diag)
      grad_mu <- gr_val
      grad_L_mat <- tcrossprod(gr_val, eps)
      grad_diag <- diag(grad_L_mat) * exp(L_diag) + 1
      grad_off <- grad_L_mat[lower.tri(grad_L_mat)]

    } else {
      elbo_history[t] <- fn_val + sum(omega)
      grad_mu <- gr_val
      grad_omega <- gr_val * (eps * sigma) + 1
    }

    # 4. Adamによるパラメータ更新
    if (fullrank) {
      # muの更新
      m_mu <- beta1 * m_mu + (1 - beta1) * grad_mu
      v_mu <- beta2 * v_mu + (1 - beta2) * (grad_mu^2)
      m_hat_mu <- m_mu / (1 - beta1^t)
      v_hat_mu <- v_mu / (1 - beta2^t)
      mu <- mu + alpha * m_hat_mu / (sqrt(v_hat_mu) + epsilon)

      # L_diagの更新
      m_diag <- beta1 * m_diag + (1 - beta1) * grad_diag
      v_diag <- beta2 * v_diag + (1 - beta2) * (grad_diag^2)
      m_hat_diag <- m_diag / (1 - beta1^t)
      v_hat_diag <- v_diag / (1 - beta2^t)
      L_diag <- L_diag + alpha * m_hat_diag / (sqrt(v_hat_diag) + epsilon)

      # L_offの更新
      if (length(L_off) > 0) {
        m_off <- beta1 * m_off + (1 - beta1) * grad_off
        v_off <- beta2 * v_off + (1 - beta2) * (grad_off^2)
        m_hat_off <- m_off / (1 - beta1^t)
        v_hat_off <- v_off / (1 - beta2^t)
        L_off <- L_off + alpha * m_hat_off / (sqrt(v_hat_off) + epsilon)
      }
    } else {
      # 既存のMean-fieldの更新処理
      m_mu <- beta1 * m_mu + (1 - beta1) * grad_mu
      v_mu <- beta2 * v_mu + (1 - beta2) * (grad_mu^2)
      m_hat_mu <- m_mu / (1 - beta1^t)
      v_hat_mu <- v_mu / (1 - beta2^t)
      mu <- mu + alpha * m_hat_mu / (sqrt(v_hat_mu) + epsilon)

      m_omega <- beta1 * m_omega + (1 - beta1) * grad_omega
      v_omega <- beta2 * v_omega + (1 - beta2) * (grad_omega^2)
      m_hat_omega <- m_omega / (1 - beta1^t)
      v_hat_omega <- v_omega / (1 - beta2^t)
      omega <- omega + alpha * m_hat_omega / (sqrt(v_hat_omega) + epsilon)
    }

    # --- 収束判定とprint_freqの処理 (既存のコードそのまま) ---
    if (print_freq > 0 && t %% print_freq == 0) {
      cat(sprintf("Iter %d: Approx ELBO = %.2f\n", t, elbo_history[t]))
    }

    check_start <- max(min_iter, 2 * window_size)
    if (t > check_start && t %% 10 == 0) {
      med_prev <- median(elbo_history[(t - 2 * window_size + 1):(t - window_size)])
      med_curr <- median(elbo_history[(t - window_size + 1):t])
      if (abs(med_curr - med_prev) / abs(med_prev) < tol_rel_obj) {
        cat(sprintf("Converged at iteration %d\n", t))
        converged <- TRUE
        break
      }
    }
  }

  elbo_history <- elbo_history[1:t]
  calc_window <- min(t, window_size)
  elbo_final <- median(elbo_history[(t - calc_window + 1):t])

  # --- 事後サンプルの生成 ---
  cat("Generating posterior samples from variational distribution...\n")
  fit_matrix <- matrix(NA, nrow = num_samples, ncol = P)

  if (fullrank) {
    L <- matrix(0, P, P)
    L[lower.tri(L)] <- L_off
    diag(L) <- exp(L_diag)
    for (i in 1:num_samples) {
      fit_matrix[i, ] <- mu + as.vector(L %*% rnorm(P))
    }
  } else {
    sigma <- exp(omega)
    for (i in 1:num_samples) {
      fit_matrix[i, ] <- mu + sigma * rnorm(P)
    }
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
