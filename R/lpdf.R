#' Log-density and log-mass functions
#'
#' `lpdf` is an environment containing log-density / log-mass functions
#' used to build `log_prob()` functions.
#'
#' @export
lpdf <- new.env(parent = baseenv())

lpdf$normal <- function(x, mean, sd) {
  sum(dnorm(x, mean = mean, sd = sd, log = TRUE))
}

lpdf$lognormal <- function(x, meanlog, sdlog) {
  sum(suppressWarnings(dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE)))
}

lpdf$exponential <- function(x, rate) {
  sum(dexp(x, rate,log=TRUE))
}

lpdf$beta <- function(x,a,b){
  sum(dbeta(x,a,b,log=TRUE))
}

lpdf$gamma <- function(x, shape, rate) {
  sum(suppressWarnings(dgamma(x, shape = shape, rate = rate, log = TRUE)))
}

lpdf$inverse_gamma <- function(x, shape, scale) {
  sum(shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x)
}

lpdf$cauchy <- function(x, location, scale) {
  sum(suppressWarnings(dcauchy(x, location = location, scale = scale, log = TRUE)))
}

lpdf$student_t <- function(x, df, mu = 0, sigma = 1) {
  sum(suppressWarnings(dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma)))
}

lpdf$laplace <- function(x, location = 0, scale = 1) {
  sum(-log(2 * scale) - abs(x - location) / scale)
}

lpdf$weibull <- function(x, shape, scale) {
  sum(suppressWarnings(dweibull(x, shape = shape, scale = scale, log = TRUE)))
}

lpdf$uniform <- function(x,a,b){
  sum(dunif(x,a,b,log=TRUE))
}

lpdf$lkj_corr <- function(Omega, eta = 1) {
  if (eta == 1) return(0)
  safe_Omega <- Omega + diag(1e-8, nrow(Omega))
  U <- chol(safe_Omega)
  return((eta - 1) * 2 * sum(log(diag(U))))
}

lpdf$lkj_CF_corr <- function(CF_Omega, eta = 1) {
  if (eta == 1) return(0)
  return((eta - 1) * 2 * sum(log(diag(CF_Omega))))
}


lpdf$bernoulli <- function(x,prob){
  sum(dbinom(x,1,prob,log=TRUE))
}

lpdf$binomial <- function(x,size,prob) {
  sum(dbinom(x,size,prob,log=TRUE))
}

lpdf$bernoulli_logit <- function(x, mu) {
  sum(x * (-math$log1p_exp(-mu)) + (1 - x) * (-math$log1p_exp(mu)))
}

lpdf$binomial_logit <- function(x, size, mu) {
  sum(lchoose(size, x) + x * (-math$log1p_exp(-mu)) + (size - x) * (-math$log1p_exp(mu)))
}

lpdf$poisson <- function(x,mean) {
  sum(dpois(x,mean,log=TRUE))
}

lpdf$neg_binomial <- function(x, size, prob) {
  sum(suppressWarnings(dnbinom(x, size = size, prob = prob, log = TRUE)))
}

lpdf$neg_binomial_2 <- function(x, mu, size) {
  sum(suppressWarnings(dnbinom(x, size = size, mu = mu, log = TRUE)))
}

lpdf$categorical <- function(x, prob) {
  if (is.matrix(prob)) {
    sum(log(prob[cbind(seq_along(x), x)]))
  } else {
    sum(log(prob[x]))
  }
}

lpdf$multinomial <- function(x, size, prob) {
  sum(suppressWarnings(dmultinom(x, size = size, prob = prob, log = TRUE)))
}

lpdf$ordered_logistic <- function(x, eta, cutpoints) {
  N <- length(x)
  K <- length(cutpoints) + 1

  if (length(eta) == 1 && N > 1) {
    eta_vec <- eta[1] * rep(1, N)
  } else {
    eta_vec <- eta
  }

  lp <- 0
  for (i in 1:N) {
    y_i <- x[i]
    eta_i <- eta_vec[i]
    if (y_i == 1) {
      lp <- lp - math$log1p_exp(-(cutpoints[1] - eta_i))
    } else if (y_i == K) {
      lp <- lp - math$log1p_exp(cutpoints[K - 1] - eta_i)
    } else {
      A <- -math$log1p_exp(-(cutpoints[y_i] - eta_i))
      B <- -math$log1p_exp(-(cutpoints[y_i-1] - eta_i))
      lp <- lp + A + math$log1m_exp(B - A)
    }
  }
  return(lp)
}

lpdf$dirichlet <- function(x, alpha) {
  sum((alpha - 1) * log(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha))
}

lpdf$multi_normal_CF <- function(x, mean, sd, CF_Omega) {
  L_Sigma <- diag(sd) %*% CF_Omega

  N <- if (is.matrix(x)) nrow(x) else 1
  K <- ncol(L_Sigma)

  # math環境の関数を使用 (math$log_det_chol 等が登録済みである前提)
  log_det <- 2 * sum(log(diag(L_Sigma)))

  if (is.matrix(x)) {
    resid_t <- t(x) - mean
    z <- solve(L_Sigma, resid_t)
    quad_form <- sum(z^2)
    return(-0.5 * (N * K * log(2 * pi) + N * log_det + quad_form))
  } else {
    z <- solve(L_Sigma, x - mean)
    quad_form <- sum(z^2)
    return(-0.5 * (K * log(2 * pi) + log_det + quad_form))
  }
}

lpdf$multi_normal <- function(x, mean, Sigma) {
  K <- nrow(Sigma)
  eps_mat <- diag(diag(Sigma) * 1e-6 + 1e-8, K)
  safe_Sigma <- Sigma + eps_mat

  U <- chol(safe_Sigma)
  log_det <- math$log_det_chol(U)
  L <- t(U)

  if (is.matrix(x)) {
    n <- nrow(x)
    k <- ncol(x)
    resid_t <- t(x) - mean
    z <- solve(L, resid_t)
    quad_form <- sum(z^2)
    return(-0.5 * (n * k * log(2 * pi) + n * log_det + quad_form))
  } else {
    k <- length(x)
    # ベクトルに対する計算
    quad_form <- math$quad_form_chol(x - mean, L)
    return(-0.5 * (k * log(2 * pi) + log_det + quad_form))
  }
}

lpdf$sum_to_zero_mvnorm <- function(x, sigma = 1, K = NULL) {
  if(is.null(K)) K <- length(x)
  lp <- -0.5*(K - 1)*log(2 * pi) -(K - 1)*log(sigma) - 0.5*sum(x^2)/sigma^2
  return(lp)
}

lpdf$centered_tri_mvnormal <- function(x, sigma = 1) {
  R <- nrow(x)
  C <- ncol(x)
  max_d <- min(C, R - 1)

  lp <- 0
  if (length(sigma) == 1) {
    df <- sum(R - (1:max_d))
    lp <- -0.5 * df * log(2 * pi) - df * log(sigma) - 0.5 * sum(x^2) / sigma^2
  } else {
    # 次元ごとに異なる sigma を指定する場合
    for (d in 1:max_d) {
      x_d <- x[d:R, d]
      K <- length(x_d)
      sd_d <- sigma[d]
      lp <- lp - 0.5*(K-1)*log(2 * pi) - (K-1)*log(sd_d) - 0.5*sum(x_d^2)/sd_d^2
    }
  }
  return(lp)
}
