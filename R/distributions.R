#' @description Normal log-probability density function.
#' @param x Vector of quantiles.
#' @param mean Vector of means.
#' @param sd Vector of standard deviations.
#' @return The sum of the log-density.
#' @export
normal_lpdf <- function(x, mean, sd) {
  sum(dnorm(x, mean = mean, sd = sd, log = TRUE))
}

#' @description Lognormal log-probability density function.
#' @param x Vector of quantiles.
#' @param meanlog Mean of the distribution on the log scale.
#' @param sdlog Standard deviation of the distribution on the log scale.
#' @return The sum of the log-density.
#' @export
lognormal_lpdf <- function(x, meanlog, sdlog) {
  sum(suppressWarnings(dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE)))
}

#' @description Exponential log-probability density function.
#' @param x Vector of quantiles.
#' @param rate Vector of rates.
#' @return The sum of the log-density.
#' @export
exponential_lpdf <- function(x, rate) {
  sum(dexp(x, rate, log = TRUE))
}

#' @description Beta log-probability density function.
#' @param x Vector of quantiles.
#' @param a Shape parameter alpha.
#' @param b Shape parameter beta.
#' @return The sum of the log-density.
#' @export
beta_lpdf <- function(x, a, b) {
  sum(dbeta(x, a, b, log = TRUE))
}

#' @description Gamma log-probability density function.
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param rate Rate parameter.
#' @return The sum of the log-density.
#' @export
gamma_lpdf <- function(x, shape, rate) {
  sum(suppressWarnings(dgamma(x, shape = shape, rate = rate, log = TRUE)))
}

#' @description Inverse-gamma log-probability density function.
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @export
inverse_gamma_lpdf <- function(x, shape, scale) {
  sum(shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x)
}

#' @description Cauchy log-probability density function.
#' @param x Vector of quantiles.
#' @param location Location parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @export
cauchy_lpdf <- function(x, location, scale) {
  sum(suppressWarnings(dcauchy(x, location = location, scale = scale, log = TRUE)))
}

#' @description Student-t log-probability density function.
#' @param x Vector of quantiles.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#' @return The sum of the log-density.
#' @export
student_t_lpdf <- function(x, df, mu = 0, sigma = 1) {
  sum(suppressWarnings(dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma)))
}

#' @description Laplace log-probability density function.
#' @param x Vector of quantiles.
#' @param location Location parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @export
laplace_lpdf <- function(x, location = 0, scale = 1) {
  sum(-log(2 * scale) - abs(x - location) / scale)
}

#' @description Weibull log-probability density function.
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @export
weibull_lpdf <- function(x, shape, scale) {
  sum(suppressWarnings(dweibull(x, shape = shape, scale = scale, log = TRUE)))
}

#' @description Uniform log-probability density function.
#' @param x Vector of quantiles.
#' @param a Lower limit of the distribution.
#' @param b Upper limit of the distribution.
#' @return The sum of the log-density.
#' @export
uniform_lpdf <- function(x, a, b) {
  sum(dunif(x, a, b, log = TRUE))
}

#' @description LKJ correlation log-probability density function.
#' @param Omega Correlation matrix.
#' @param eta Shape parameter.
#' @return The log-density.
#' @export
lkj_corr_lpdf <- function(Omega, eta = 1) {
  if (eta == 1) return(0)
  safe_Omega <- Omega + diag(1e-8, nrow(Omega))
  U <- chol(safe_Omega)
  return((eta - 1) * 2 * sum(log(diag(U))))
}

#' @description LKJ correlation log-probability density function for Cholesky factors.
#' @param CF_Omega Cholesky factor of a correlation matrix.
#' @param eta Shape parameter.
#' @return The log-density.
#' @export
lkj_CF_corr_lpdf <- function(CF_Omega, eta = 1) {
  if (eta == 1) return(0)
  return((eta - 1) * 2 * sum(log(diag(CF_Omega))))
}

#' @description Bernoulli log-probability mass function.
#' @param x Vector of binary outcomes (0 or 1).
#' @param prob Probability of success.
#' @return The sum of the log-probability.
#' @export
bernoulli_lpdf <- function(x, prob) {
  sum(dbinom(x, 1, prob, log = TRUE))
}

#' @description Binomial log-probability mass function.
#' @param x Vector of quantiles (number of successes).
#' @param size Number of trials.
#' @param prob Probability of success.
#' @return The sum of the log-probability.
#' @export
binomial_lpdf <- function(x, size, prob) {
  sum(dbinom(x, size, prob, log = TRUE))
}

#' @description Bernoulli log-probability mass function with logit parameterization.
#' @param x Vector of binary outcomes (0 or 1).
#' @param mu Linear predictor (logit probability).
#' @return The sum of the log-probability.
#' @export
bernoulli_logit_lpdf <- function(x, mu) {
  log1p_exp <- function(x) {
    max_val <- (x + sqrt(x^2 + 1e-10)) / 2
    return(max_val + log(exp(x - max_val) + exp(-max_val)))
  }
  sum(x * (-log1p_exp(-mu)) + (1 - x) * (-log1p_exp(mu)))
}

#' @description Binomial log-probability mass function with logit parameterization.
#' @param x Vector of quantiles (number of successes).
#' @param size Number of trials.
#' @param mu Linear predictor (logit probability).
#' @return The sum of the log-probability.
#' @export
binomial_logit_lpdf <- function(x, size, mu) {
  log1p_exp <- function(x) {
    max_val <- (x + sqrt(x^2 + 1e-10)) / 2
    return(max_val + log(exp(x - max_val) + exp(-max_val)))
  }
  sum(lchoose(size, x) + x * (-log1p_exp(-mu)) + (size - x) * (-log1p_exp(mu)))
}

#' @description Poisson log-probability mass function.
#' @param x Vector of quantiles.
#' @param mean Expected value.
#' @return The sum of the log-probability.
#' @export
poisson_lpdf <- function(x, mean) {
  sum(dpois(x, mean, log = TRUE))
}

#' @description Negative binomial log-probability mass function.
#' @param x Vector of quantiles.
#' @param size Target for number of successful trials.
#' @param prob Probability of success in each trial.
#' @return The sum of the log-probability.
#' @export
neg_binomial_lpdf <- function(x, size, prob) {
  sum(suppressWarnings(dnbinom(x, size = size, prob = prob, log = TRUE)))
}

#' @description Negative binomial log-probability mass function (alternative parameterization).
#' @param x Vector of quantiles.
#' @param mu Mean parameter.
#' @param size Dispersion parameter.
#' @return The sum of the log-probability.
#' @export
neg_binomial_2_lpdf <- function(x, mu, size) {
  sum(suppressWarnings(dnbinom(x, size = size, mu = mu, log = TRUE)))
}

#' @description Categorical log-probability mass function.
#' @param x Vector of categorical outcomes.
#' @param prob Vector or matrix of probabilities.
#' @return The sum of the log-probability.
#' @export
categorical_lpdf <- function(x, prob) {
  if (is.matrix(prob)) {
    sum(log(prob[cbind(seq_along(x), x)]))
  } else {
    sum(log(prob[x]))
  }
}

#' @description Multinomial log-probability mass function.
#' @param x Vector of counts.
#' @param size Total number of trials.
#' @param prob Vector of probabilities for each category.
#' @return The sum of the log-probability.
#' @export
multinomial_lpdf <- function(x, size, prob) {
  sum(suppressWarnings(dmultinom(x, size = size, prob = prob, log = TRUE)))
}

#' @description Ordered logistic log-probability mass function.
#' @param x Vector of ordered categorical outcomes.
#' @param eta Linear predictor.
#' @param cutpoints Vector of cutpoints.
#' @return The sum of the log-probability.
#' @export
ordered_logistic_lpdf <- function(x, eta, cutpoints) {
  N <- length(x)
  K <- length(cutpoints) + 1

  log1p_exp <- function(x) {
    max_val <- (x + sqrt(x^2 + 1e-10)) / 2
    return(max_val + log(exp(x - max_val) + exp(-max_val)))
  }
  log1m_exp <- function(x) {
    return(log(1 - exp(x) + 1e-10))
  }

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
      lp <- lp - log1p_exp(-(cutpoints[1] - eta_i))
    } else if (y_i == K) {
      lp <- lp - log1p_exp(cutpoints[K - 1] - eta_i)
    } else {
      A <- -log1p_exp(-(cutpoints[y_i] - eta_i))
      B <- -log1p_exp(-(cutpoints[y_i-1] - eta_i))
      lp <- lp + A + log1m_exp(B - A)
    }
  }
  return(lp)
}

#' @description Dirichlet log-probability density function.
#' @param x Vector or matrix of simplexes.
#' @param alpha Vector of concentration parameters.
#' @return The sum of the log-density.
#' @export
dirichlet_lpdf <- function(x, alpha) {
  sum((alpha - 1) * log(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha))
}

#' @description Multivariate normal log-probability density function parameterized by Cholesky factor of correlation matrix.
#' @param x Vector or matrix of quantiles.
#' @param mean Vector or matrix of means.
#' @param sd Vector of standard deviations.
#' @param CF_Omega Cholesky factor of the correlation matrix.
#' @return The sum of the log-density.
#' @export
multi_normal_CF_lpdf <- function(x, mean, sd, CF_Omega) {
  L_Sigma <- diag(sd) %*% CF_Omega

  N <- if (is.matrix(x)) nrow(x) else 1
  K <- ncol(L_Sigma)

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

#' @description Multivariate normal log-probability density function.
#' @param x Vector or matrix of quantiles.
#' @param mean Vector or matrix of means.
#' @param Sigma Covariance matrix.
#' @return The sum of the log-density.
#' @export
multi_normal_lpdf <- function(x, mean, Sigma) {

  log_det_chol <- function(L) {
    return(2 * sum(log(diag(L))))
  }
  quad_form_chol <- function(x, L) {
    z <- solve(L, x)
    return(sum(z^2))
  }

  K <- nrow(Sigma)
  eps_mat <- diag(diag(Sigma) * 1e-6 + 1e-8, K)
  safe_Sigma <- Sigma + eps_mat

  U <- chol(safe_Sigma)
  log_det <- log_det_chol(U)
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
    quad_form <- quad_form_chol(x - mean, L)
    return(-0.5 * (k * log(2 * pi) + log_det + quad_form))
  }
}

#' @description Sum-to-zero multivariate normal log-probability density function.
#' @param x Vector of quantiles.
#' @param sigma Standard deviation parameter.
#' @param K Dimension. If NULL, inferred from the length of x.
#' @return The log-density.
#' @export
sum_to_zero_mvnorm_lpdf <- function(x, sigma = 1, K = NULL) {
  if (is.null(K)) K <- length(x)
  lp <- -0.5 * (K - 1) * log(2 * pi) - (K - 1) * log(sigma) - 0.5 * sum(x^2) / sigma^2
  return(lp)
}

#' @description Centered triangular multivariate normal log-probability density function.
#' @param x Matrix of quantiles.
#' @param sigma Standard deviation parameter(s).
#' @return The log-density.
#' @export
centered_tri_mvnormal_lpdf <- function(x, sigma = 1) {
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
      lp <- lp - 0.5 * (K - 1) * log(2 * pi) - (K - 1) * log(sd_d) - 0.5 * sum(x_d^2) / sd_d^2
    }
  }
  return(lp)
}
