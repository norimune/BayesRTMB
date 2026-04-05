#' Preconditioned MALA (Metropolis-Adjusted Langevin Algorithm)
#'
#' @param ad_obj AD object from RTMB
#' @param mu MAP estimates (mean vector)
#' @param Sigma Covariance matrix from Laplace approximation
#' @param sampling Number of sampling iterations
#' @param warmup Number of warmup iterations
#' @param delta Target acceptance rate (Default: 0.574 for MALA)
#' @param chain Chain ID (for progress reporting)
#' @return A list containing samples and diagnostic statistics
MALA_method <- function(ad_obj, mu, Sigma, sampling, warmup, delta = 0.574, chain = 1) {
  iter <- sampling + warmup
  P <- length(mu)

  # 共分散行列のコレスキー分解 (下三角行列 L)
  L <- tryCatch(t(chol(Sigma)), error = function(e) NULL)
  use_chol <- !is.null(L)
  if (!use_chol) {
    warning("Cholesky factorization failed. Falling back to eigenvalue decomposition.")
    eig <- eigen(Sigma, symmetric = TRUE)
    val <- pmax(eig$values, 1e-10) # ゼロ割を防ぐ
    L <- eig$vectors %*% diag(sqrt(val))
    L_inv <- diag(1 / sqrt(val)) %*% t(eig$vectors)
  }

  # Dual Averaging のハイパーパラメータ
  epsilon <- 0.1 / (P^(1/6)) # 初期ステップサイズを次元数に応じて小さく開始
  mu_da <- log(10 * epsilon)
  gamma_da <- 0.05
  t0_da <- 10
  kappa_da <- 0.75
  epsilon_bar <- 1.0
  H_bar <- 0

  # 保存用配列
  draws <- matrix(NA, nrow = iter, ncol = P)
  accept_stat <- numeric(iter)
  lp <- numeric(iter)

  # 初期状態の設定
  theta <- mu
  U <- tryCatch(ad_obj$fn(theta), error = function(e) Inf)

  # 境界値等でMAPの目的関数が計算できない場合の退避処理
  if (!is.finite(U)) {
    cat(sprintf("Chain %d: Initial objective is non-finite. Perturbing from MAP...\n", chain))
    for (try_idx in 1:50) {
      theta <- mu + rnorm(P, mean = 0, sd = 0.1)
      U <- tryCatch(ad_obj$fn(theta), error = function(e) Inf)
      if (is.finite(U)) break
    }
  }

  grad_U <- tryCatch(ad_obj$gr(theta), error = function(e) rep(0, P))

  draws[1, ] <- theta
  lp[1] <- if (is.finite(U)) -U else NA
  accept_stat[1] <- 1

  cat(sprintf("Chain %d: Warmup & Sampling (Preconditioned MALA)...\n", chain))
  pb <- txtProgressBar(min = 0, max = iter, style = 3)

  for (i in 2:iter) {
    # --- 1. 提案候補の生成 ---
    mu_prop <- theta - 0.5 * (epsilon^2) * as.vector(Sigma %*% grad_U)

    z <- rnorm(P)
    d_fwd <- as.vector(L %*% z)
    theta_star <- mu_prop + epsilon * d_fwd

    # --- 2. 提案位置での評価 ---
    U_star <- tryCatch(ad_obj$fn(theta_star), error = function(e) Inf)

    if (is.infinite(U_star) || is.na(U_star)) {
      log_alpha <- -Inf
    } else {
      grad_U_star <- tryCatch(ad_obj$gr(theta_star), error = function(e) rep(0, P))
      mu_prop_star <- theta_star - 0.5 * (epsilon^2) * as.vector(Sigma %*% grad_U_star)

      # --- 3. 提案確率の計算 ---
      log_q_fwd <- -0.5 * sum(z^2)

      d_rev <- theta - mu_prop_star

      if (use_chol) {
        z_rev <- forwardsolve(L, d_rev) / epsilon
      } else {
        z_rev <- as.vector(L_inv %*% d_rev) / epsilon
      }

      log_q_rev <- -0.5 * sum(z_rev^2)

      # 採択確率 (対数スケール)
      log_alpha <- -U_star + U + log_q_rev - log_q_fwd
    }

    alpha <- min(1, exp(log_alpha))
    if (is.na(alpha)) alpha <- 0

    # --- 4. 採択判定 ---
    if (runif(1) < alpha) {
      theta <- theta_star
      U <- U_star
      grad_U <- grad_U_star
    }

    # 結果の記録
    draws[i, ] <- theta
    lp[i] <- if (is.finite(U)) -U else NA
    accept_stat[i] <- alpha

    # --- 5. Dual Averaging によるステップサイズ (epsilon) の調整 ---
    if (i <= warmup) {
      eta <- 1 / (i + t0_da)
      H_bar <- (1 - eta) * H_bar + eta * (delta - alpha)
      log_epsilon <- mu_da - sqrt(i) / gamma_da * H_bar
      epsilon <- exp(log_epsilon)

      eta_bar <- i^(-kappa_da)
      epsilon_bar <- exp((1 - eta_bar) * log(epsilon_bar) + eta_bar * log_epsilon)
    } else if (i == warmup + 1) {
      epsilon <- epsilon_bar
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("\n")

  return(list(
    fit = draws,
    lp = lp,
    accept = accept_stat,
    epsilon = epsilon
  ))
}
