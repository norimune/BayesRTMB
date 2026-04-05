#' Hybrid IMH and Preconditioned MALA sampling using Laplace approximation
#'
#' This function performs a mixture of Independent Metropolis-Hastings (IMH)
#' and Preconditioned Metropolis-Adjusted Langevin Algorithm (pMALA).
#' The IMH step uses a multivariate t-distribution centered at the MAP estimate.
#' The pMALA step uses gradient information preconditioned by the inverse Hessian.
#'
#' @param model An RTMB automatic differentiation objective object (`ad_obj`).
#' @param sampling Integer. Number of sampling iterations.
#' @param warmup Integer. Number of warmup iterations.
#' @param chain Integer. The ID of the current MCMC chain.
#' @param update_progress An optional callback function used to update the progress bar. Default is NULL.
#' @param df Numeric. Degrees of freedom for the multivariate t-distribution proposal in IMH. Default is 4.
#' @param mix_ratio Numeric. The probability of choosing the IMH step. Default is 0.8 (80% IMH, 20% MALA).
#' @param eps_mala Numeric. The step size parameter for the pMALA step. Default is 0.5.
#'
#' @return A list containing MCMC results.
#' @export
IMH_method <- function(model, sampling, warmup, chain, update_progress = NULL, df = 4, mix_ratio = 0.8, eps_mala = 0.1) {
  iter <- sampling + warmup
  P <- length(model$par)

  # ターゲット関数（負の対数事後確率）
  fn <- function(x) {
    val <- tryCatch(model$fn(x), error = function(e) Inf)
    if (is.na(val) || is.nan(val)) return(Inf)
    return(val)
  }

  # 勾配関数（負の対数事後確率の勾配を反転させ、対数事後確率の勾配にする）
  gr_fn <- function(x) {
    g <- tryCatch(model$gr(x), error = function(e) rep(NaN, P))
    if (any(is.nan(g))) return(rep(NaN, P))
    return(-g)
  }

  # MAP推定とヘッシアンの計算
  opt <- nlminb(start = model$par, objective = fn, gradient = model$gr)
  sdr <- tryCatch(RTMB::sdreport(model), error = function(e) NULL)

  Map <- opt$par
  if (!is.null(sdr) && !is.null(sdr$cov.fixed)) {
    Cov <- sdr$cov.fixed
  } else {
    Hess <- optimHess(Map, fn, model$gr)
    Cov <- solve(Hess)
  }

  # 正定値性の保証
  Cov <- Cov + diag(1e-8, P)
  L <- t(chol(Cov))
  invL <- solve(L) # MALAの後退確率計算用

  para_fixed <- array(NA, dim = c(iter, P))
  para_fixed[1, ] <- Map

  lp <- numeric(iter)
  lp[1] <- -fn(Map)

  # 勾配の初期状態
  gr_curr <- gr_fn(Map)

  accept <- numeric(iter)
  accept[1] <- 1

  # t分布の対数密度計算関数
  calc_log_q_t <- function(d2, P, df) {
    -0.5 * (df + P) * log(1 + d2 / df)
  }

  z_curr_raw <- solve(L, para_fixed[1, ] - Map)
  d2_curr_imh <- sum(z_curr_raw^2)
  lq_curr_imh <- calc_log_q_t(d2_curr_imh, P, df)

  for (i in 2:iter) {
    current <- para_fixed[i-1, ]
    lp_curr <- lp[i-1]

    # 手法のランダム選択
    is_imh <- runif(1) < mix_ratio

    if (is_imh) {
      # ==========================================
      # 1. IMH Step
      # ==========================================
      z <- rnorm(P)
      u <- rchisq(1, df = df)
      scale_factor <- sqrt(df / u)

      z_scaled <- z * scale_factor*1.5
      propose <- Map + as.vector(L %*% z_scaled)

      d2_prop <- sum(z_scaled^2)
      lq_prop_imh <- calc_log_q_t(d2_prop, P, df)

      lp_prop <- -fn(propose)

      if (is.infinite(lp_prop) || lp_prop < -1e20) {
        accept[i] <- 0
      } else {
        log_ratio <- (lp_prop - lq_prop_imh) - (lp_curr - lq_curr_imh)

        if (log(runif(1)) < log_ratio) {
          current <- propose
          lp_curr <- lp_prop
          lq_curr_imh <- lq_prop_imh
          gr_curr <- gr_fn(propose) # MALA用に勾配を更新
          accept[i] <- 1
        } else {
          accept[i] <- 0
        }
      }

    } else {
      # ==========================================
      # 2. pMALA Step
      # ==========================================
      # 勾配が計算できない場合は棄却
      if (any(is.nan(gr_curr))) {
        accept[i] <- 0
        para_fixed[i, ] <- current
        lp[i] <- lp_curr
        next
      }

      # 前進提案 (Forward proposal)
      z <- rnorm(P)
      mu_fwd <- current + 0.5 * (eps_mala^2) * as.vector(Cov %*% gr_curr)
      propose <- mu_fwd + eps_mala * as.vector(L %*% z)

      lp_prop <- -fn(propose)

      if (is.infinite(lp_prop) || lp_prop < -1e20) {
        accept[i] <- 0
      } else {
        gr_prop <- gr_fn(propose)

        if (any(is.nan(gr_prop))) {
          accept[i] <- 0
        } else {
          # 前進遷移確率: q(x_prop | x_curr)
          # マハラノビス距離の平方部分はそのまま sum(z^2)
          log_q_fwd <- -0.5 * sum(z^2)

          # 後退遷移確率: q(x_curr | x_prop)
          mu_bwd <- propose + 0.5 * (eps_mala^2) * as.vector(Cov %*% gr_prop)
          z_bwd <- as.vector(invL %*% (current - mu_bwd)) / eps_mala
          log_q_bwd <- -0.5 * sum(z_bwd^2)

          # MH比の計算
          log_ratio <- (lp_prop - lp_curr) + (log_q_bwd - log_q_fwd)

          if (log(runif(1)) < log_ratio) {
            current <- propose
            lp_curr <- lp_prop
            gr_curr <- gr_prop

            # IMH用に現在位置の提案密度を更新
            z_curr_raw <- as.vector(invL %*% (current - Map))
            lq_curr_imh <- calc_log_q_t(sum(z_curr_raw^2), P, df)

            accept[i] <- 1
          } else {
            accept[i] <- 0
          }
        }
      }
    }

    para_fixed[i, ] <- current
    lp[i] <- lp_curr

    if (i %% 100 == 0) {
      if (!is.null(update_progress)) update_progress()
      else cat(paste0("chain ", chain, ": iter ", i, ifelse(i <= warmup, " warmup", " sampling"), "\n"))
    }
  }

  return(list(
    para_fixed = para_fixed,
    lp         = lp,
    accept     = accept,
    treedepth  = rep(0, iter),
    eps        = 1
  ))
}
