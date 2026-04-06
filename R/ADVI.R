#' Automatic Differentiation Variational Inference (ADVI)
#'
#' @param model An RTMB objective function object (`ad_obj`).
#' @param par_list A list defining the structure of parameters to be estimated.
#' @param pl_full A list defining the full structure of parameters including random effects.
#' @param iter Integer; maximum number of iterations for the optimization. Default is 10000.
#' @param tol_rel_obj Numeric; relative tolerance for the ELBO to determine convergence. Default is 0.01.
#' @param window_size Integer; size of the moving window to calculate the median ELBO for the convergence check. Default is 100.
#' @param num_samples Integer; number of posterior draws to generate after optimization. Default is 4000.
#' @param chains Integer; number of chains to structure the output array (for compatibility). Default is 1.
#' @param alpha Numeric; learning rate (step size) for the Adam optimizer. Default is 0.01.
#' @param laplace Logical; whether Laplace approximation is used. Default is FALSE.
#' @param print_freq Integer; frequency of printing progress to the console. Set to 0 to disable. Default is 500.
#' @param min_iter Integer; minimum number of iterations to run before checking for convergence. Default is 1000.
#' @param fullrank Logical; whether to use a full-rank multivariate normal distribution for the approximation. Default is FALSE.
#' @return A list containing `fit`, `random_fit`, `elbo_history`, and `elbo_final`.
ADVI_method <- function(model, par_list, pl_full,
                        iter = 10000, min_iter = 1000, tol_rel_obj = 0.001,
                        window_size = 100, num_samples = 1000,
                        chains = 1, alpha = 0.01, laplace = FALSE,
                        print_freq = 500,
                        fullrank = FALSE) {

  # --- 1. 変数の初期化 ---
  P <- length(model$par)
  mu <- model$par
  par_names <- names(model$par)

  # サンプル数の計算
  samples_per_chain <- ceiling(num_samples / chains)
  actual_num_samples <- samples_per_chain * chains

  # --- 2. Adamと近似分布のハイパーパラメータ初期化 ---
  beta1 <- 0.9
  beta2 <- 0.999
  epsilon <- 1e-8
  entropy_const <- (P / 2) * log(2 * pi * exp(1)) # エントロピーの定数項を追加

  cat("Starting ADVI optimization with Adam...\n")

  if (fullrank) {
    # 修正: 初期分散を小さくする
    L_diag <- rep(-2, P)
    L_off <- rep(0, P * (P - 1) / 2)

    m_mu <- rep(0, P); v_mu <- rep(0, P)
    m_diag <- rep(0, P); v_diag <- rep(0, P)
    m_off <- rep(0, length(L_off)); v_off <- rep(0, length(L_off))
  } else {
    # 修正: 初期分散を小さくする
    omega <- rep(-2, P)
    m_mu <- rep(0, P); v_mu <- rep(0, P)
    m_omega <- rep(0, P); v_omega <- rep(0, P)
  }

  elbo_history <- numeric(iter)
  converged <- FALSE

  # --- 3. 最適化ループ ---
  for (t in 1:iter) {
    eps <- rnorm(P)

    # 3-1. パラメータのサンプリング
    if (fullrank) {
      L <- matrix(0, P, P)
      L[lower.tri(L)] <- L_off
      diag(L) <- exp(L_diag)
      theta <- mu + as.vector(L %*% eps)
    } else {
      sigma <- exp(omega)
      theta <- mu + sigma * eps
    }

    # 3-2. 尤度と勾配の評価
    fn_val <- -model$fn(theta)
    gr_val <- as.vector(-model$gr(theta)) # 修正: 確実にベクトル化する

    # NA/NaNハンドリングの追加
    if (!is.finite(fn_val) || any(!is.finite(gr_val))) {
      # 崩壊を防ぐため、このイテレーションの更新をスキップ
      elbo_history[t] <- if(t > 1) elbo_history[t-1] else 0
      next
    }

    # 3-3. ELBOと勾配の計算
    if (fullrank) {
      elbo_history[t] <- fn_val + sum(L_diag) + entropy_const
      grad_mu <- gr_val

      grad_L_mat <- outer(gr_val, eps)
      grad_diag <- diag(grad_L_mat) * exp(L_diag) + 1
      grad_off <- grad_L_mat[lower.tri(grad_L_mat)]
    } else {
      elbo_history[t] <- fn_val + sum(omega) + entropy_const
      grad_mu <- gr_val
      grad_omega <- gr_val * (eps * sigma) + 1
    }

    # 3-4. Adamによるパラメータ更新
    if (fullrank) {
      m_mu <- beta1 * m_mu + (1 - beta1) * grad_mu
      v_mu <- beta2 * v_mu + (1 - beta2) * (grad_mu^2)
      mu <- mu + alpha * (m_mu / (1 - beta1^t)) / (sqrt(v_mu / (1 - beta2^t)) + epsilon)

      m_diag <- beta1 * m_diag + (1 - beta1) * grad_diag
      v_diag <- beta2 * v_diag + (1 - beta2) * (grad_diag^2)
      L_diag <- L_diag + alpha * (m_diag / (1 - beta1^t)) / (sqrt(v_diag / (1 - beta2^t)) + epsilon)

      if (length(L_off) > 0) {
        m_off <- beta1 * m_off + (1 - beta1) * grad_off
        v_off <- beta2 * v_off + (1 - beta2) * (grad_off^2)
        L_off <- L_off + alpha * (m_off / (1 - beta1^t)) / (sqrt(v_off / (1 - beta2^t)) + epsilon)
      }
    } else {
      m_mu <- beta1 * m_mu + (1 - beta1) * grad_mu
      v_mu <- beta2 * v_mu + (1 - beta2) * (grad_mu^2)
      mu <- mu + alpha * (m_mu / (1 - beta1^t)) / (sqrt(v_mu / (1 - beta2^t)) + epsilon)

      m_omega <- beta1 * m_omega + (1 - beta1) * grad_omega
      v_omega <- beta2 * v_omega + (1 - beta2) * (grad_omega^2)
      omega <- omega + alpha * (m_omega / (1 - beta1^t)) / (sqrt(v_omega / (1 - beta2^t)) + epsilon)
    }

    # 3-5. 進行状況の出力
    if (print_freq > 0 && t %% print_freq == 0) {
      cat(sprintf("Iter %d: Approx ELBO = %.2f\n", t, elbo_history[t]))
    }

    # 3-6. 収束判定
    check_start <- max(min_iter, 2 * window_size)
    if (t > check_start && t %% 10 == 0) {
      med_prev <- median(elbo_history[(t - 2 * window_size + 1):(t - window_size)])
      med_curr <- median(elbo_history[(t - window_size + 1):t])
      # 修正: 分母に1e-8を足してゼロ除算を防ぐ
      if (abs(med_curr - med_prev) / (abs(med_prev) + 1e-8) < tol_rel_obj) {
        cat(sprintf("Converged at iteration %d\n", t))
        converged <- TRUE
        break
      }
    }
  }

  if (!converged) {
    warning("ADVI did not converge within the maximum number of iterations.")
  }

  elbo_history <- elbo_history[1:t]
  calc_window <- min(t, window_size)
  elbo_final <- median(elbo_history[(t - calc_window + 1):t])

  # --- 4. 事後サンプルの生成と整形 ---
  cat("Generating posterior samples from variational distribution...\n")

  # 近似分布からのサンプリング
  fit_matrix <- matrix(NA, nrow = actual_num_samples, ncol = P)

  if (fullrank) {
    L <- matrix(0, P, P)
    L[lower.tri(L)] <- L_off
    diag(L) <- exp(L_diag)
    for (i in 1:actual_num_samples) {
      fit_matrix[i, ] <- mu + as.vector(L %*% rnorm(P))
    }
  } else {
    sigma <- exp(omega)
    for (i in 1:actual_num_samples) {
      fit_matrix[i, ] <- mu + sigma * rnorm(P)
    }
  }

  # パラメータ名から固定効果と変量効果の「フラットな名前」を確実に取り出す処理
  base_names_full <- sub("\\[.*\\]", "", pl_full$names)
  base_names_fixed <- names(par_list)

  fixed_names <- pl_full$names[base_names_full %in% base_names_fixed]
  random_names <- pl_full$names[!(base_names_full %in% base_names_fixed)]

  fixed_idx  <- which(pl_full$names %in% fixed_names)
  random_idx <- if (length(random_names) > 0) which(pl_full$names %in% random_names) else integer(0)

  P_all_true <- length(pl_full$names)
  para_final <- array(NA, dim = c(actual_num_samples, P_all_true))
  lp_final <- numeric(actual_num_samples)

  # fit_matrixのサンプルを使って対数尤度計算と変量効果の周辺化を行う
  for (i in 1:actual_num_samples) {
    zeta_sample <- fit_matrix[i, ]
    lp_final[i] <- -model$fn(zeta_sample)

    if (laplace && length(model$env$random) > 0) {
      para_list_res <- model$env$parList()
    } else {
      para_list_res <- model$env$parList(x = zeta_sample)
    }

    con_list <- to_constrained(para_list_res, par_list)
    para_final[i, ] <- unlist(con_list, use.names = FALSE)
  }

  # MCMC_Fit互換の配列に整形
  fit <- array(NA, dim = c(samples_per_chain, chains, length(fixed_names) + 1))
  dimnames(fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = c("lp", fixed_names))

  random_fit <- NULL
  if (length(random_names) > 0) {
    random_fit <- array(NA, dim = c(samples_per_chain, chains, length(random_names)))
    dimnames(random_fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = random_names)
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
