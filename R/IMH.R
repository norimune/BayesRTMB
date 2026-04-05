#' Independent Metropolis-Hastings Algorithm
#'
#' @param ad_obj AD object from RTMB
#' @param mu MAP estimates (mean vector)
#' @param Sigma Covariance matrix from Laplace approximation
#' @param sampling Number of sampling iterations
#' @param warmup Number of warmup iterations
#' @param chain Chain ID (for progress reporting)
#' @param update_progress Callback function for progress reporting
#' @return A list containing samples and diagnostic statistics
IMH_method <- function(ad_obj, mu, Sigma, sampling, warmup, chain = 1, update_progress = NULL) {
  iter <- sampling + warmup
  P <- length(mu)

  # 提案分布のための事前計算
  L <- tryCatch(t(chol(Sigma)), error = function(e) NULL)
  if (is.null(L)) {
    warning("Cholesky factorization failed. Falling back to eigenvalue decomposition.")
    eig <- eigen(Sigma, symmetric = TRUE)
    val <- pmax(eig$values, 1e-10) # ゼロ割を防ぐ
    L <- eig$vectors %*% diag(sqrt(val))
    Sigma_inv <- eig$vectors %*% diag(1 / val) %*% t(eig$vectors)
  } else {
    Sigma_inv <- chol2inv(t(L))
  }

  # 提案分布の対数密度計算関数 (定数項は相殺されるため二次形式のみ)
  calc_log_q <- function(x) {
    diff <- x - mu
    -0.5 * as.numeric(t(diff) %*% Sigma_inv %*% diff)
  }

  # 保存用配列
  draws <- matrix(NA, nrow = iter, ncol = P)
  accept_stat <- numeric(iter)
  lp <- numeric(iter)

  # 初期状態の設定 (MAP推定値からスタート)
  theta <- mu
  U <- tryCatch(ad_obj$fn(theta), error = function(e) Inf)

  # MAP位置で目的関数が計算できない場合の退避処理
  if (!is.finite(U)) {
    for (try_idx in 1:50) {
      theta <- mu + as.vector(L %*% rnorm(P)) * 0.05
      U <- tryCatch(ad_obj$fn(theta), error = function(e) Inf)
      if (is.finite(U)) break
    }
  }

  log_q_theta <- calc_log_q(theta)

  draws[1, ] <- theta
  lp[1] <- if (is.finite(U)) -U else NA
  accept_stat[1] <- 1

  for (i in 2:iter) {
    # --- 1. 提案候補の生成 ---
    theta_star <- mu + as.vector(L %*% rnorm(P))

    # --- 2. 提案位置での評価 ---
    U_star <- tryCatch(ad_obj$fn(theta_star), error = function(e) Inf)

    if (is.infinite(U_star) || is.na(U_star)) {
      log_alpha <- -Inf
    } else {
      # --- 3. 提案確率の計算 ---
      log_q_star <- calc_log_q(theta_star)
      log_alpha <- -U_star + U + log_q_theta - log_q_star
    }

    alpha <- min(1, exp(log_alpha))
    if (is.na(alpha)) alpha <- 0

    # --- 4. 採択判定 ---
    if (runif(1) < alpha) {
      theta <- theta_star
      U <- U_star
      log_q_theta <- log_q_star
    }

    # 結果の記録
    draws[i, ] <- theta
    lp[i] <- if (is.finite(U)) -U else NA
    accept_stat[i] <- alpha

    # NUTSと同様のプログレスレポート
    if (i %% 100 == 0) {
      if (!is.null(update_progress)) {
        update_progress()
      } else {
        cat(paste0("chain ", chain, ": iter ", i,
                   ifelse(i <= warmup, " warmup", " sampling"), "\n"))
      }
    }
  }

  return(list(
    para_fixed = draws,
    lp = lp,
    accept = accept_stat,
    treedepth = rep(NA, iter), # NUTS互換用ダミー
    eps = NA # NUTS互換用ダミー
  ))
}
