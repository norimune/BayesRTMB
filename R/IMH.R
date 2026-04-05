#' Independent Metropolis-Hastings (IMH) sampling using Laplace approximation
#'
#' This function performs IMH sampling using a multivariate t-distribution
#' as the proposal distribution. The proposal distribution is centered at the
#' MAP estimate, and its scale matrix is derived from the inverse Hessian.
#' This is particularly effective when the random effects are marginalized out
#' via the Laplace approximation in RTMB.
#'
#' @param model An RTMB automatic differentiation objective object (`ad_obj`).
#' @param sampling Integer. Number of sampling iterations.
#' @param warmup Integer. Number of warmup iterations.
#' @param chain Integer. The ID of the current MCMC chain.
#' @param update_progress An optional callback function used to update the progress bar. Default is NULL.
#' @param df Numeric. Degrees of freedom for the multivariate t-distribution proposal. Lower values result in heavier tails. Default is 4.
#'
#' @return A list containing:
#'   \item{para_fixed}{A matrix of posterior samples for the fixed parameters.}
#'   \item{lp}{A numeric vector of log-posterior values at each iteration.}
#'   \item{accept}{A numeric vector indicating whether the proposal was accepted (1) or rejected (0).}
#'   \item{treedepth}{A numeric vector of zeros, included for compatibility with `NUTS_method` outputs.}
#'   \item{eps}{A scalar zero, included for compatibility with `NUTS_method` outputs.}
#' @export
IMH_method <- function(model, sampling, warmup, chain, update_progress = NULL, df = 4) {
  iter <- sampling + warmup
  P <- length(model$par)

  # ターゲット関数（負の対数事後確率）
  fn <- function(x) {
    val <- tryCatch(model$fn(x), error = function(e) Inf)
    if (is.na(val) || is.nan(val)) return(Inf)
    return(val)
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

  para_fixed <- array(NA, dim = c(iter, P))
  para_fixed[1, ] <- Map # 初期値はMAP推定値

  lp <- numeric(iter)
  lp[1] <- -fn(Map)

  accept <- numeric(iter)
  accept[1] <- 1

  calc_log_q <- function(d2, P, df) {
    -0.5 * (df + P) * log(1 + d2 / df)
  }

  z_curr <- solve(L, para_fixed[1, ] - Map)
  d2_curr <- sum(z_curr^2)
  lq_curr <- calc_log_q(d2_curr, P, df)

  for (i in 2:iter) {
    current <- para_fixed[i-1, ]
    lp_curr <- lp[i-1]

    # t分布からのサンプリング
    z <- rnorm(P)
    u <- rchisq(1, df = df)
    scale_factor <- sqrt(df / u)

    z_scaled <- z * scale_factor
    propose <- Map + as.vector(L %*% z_scaled)

    d2_prop <- sum(z_scaled^2)
    lq_prop <- calc_log_q(d2_prop, P, df)

    lp_prop <- -fn(propose)

    if (is.infinite(lp_prop) || lp_prop < -1e20) {
      accept[i] <- 0
    } else {
      log_ratio <- (lp_prop - lq_prop) - (lp_curr - lq_curr)

      if (log(runif(1)) < log_ratio) {
        current <- propose
        lp_curr <- lp_prop
        lq_curr <- lq_prop
        accept[i] <- 1
      } else {
        accept[i] <- 0
      }
    }

    para_fixed[i, ] <- current
    lp[i]     <- lp_curr

    if (i %% 100 == 0) {
      if (!is.null(update_progress)) update_progress()
      else cat(paste0("chain ", chain, ": iter ", i, ifelse(i <= warmup, " warmup", " sampling"), "\n"))
    }
  }

  # NUTS_methodと戻り値の構造を揃える
  return(list(
    para_fixed = para_fixed,
    lp         = lp,
    accept     = accept,
    treedepth  = rep(0, iter),
    eps        = 0
  ))
}
