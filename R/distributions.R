#' Normal log-probability density function
#'
#' @param x Vector of quantiles.
#' @param mean Vector of means.
#' @param sd Vector of standard deviations.
#' @return The sum of the log-density.
#' @export
normal_lpdf <- function(x, mean, sd) {
  sum(dnorm(x, mean = mean, sd = sd, log = TRUE))
}

#' Lognormal log-probability density function
#'
#' @param x Vector of quantiles.
#' @param meanlog Mean of the distribution on the log scale.
#' @param sdlog Standard deviation of the distribution on the log scale.
#' @return The sum of the log-density.
#' @export
lognormal_lpdf <- function(x, meanlog, sdlog) {
  sum(suppressWarnings(dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE)))
}

#' Exponential log-probability density function
#'
#' @param x Vector of quantiles.
#' @param rate Vector of rates.
#' @return The sum of the log-density.
#' @export
exponential_lpdf <- function(x, rate) {
  sum(dexp(x, rate, log = TRUE))
}

#' Beta log-probability density function
#'
#' @param x Vector of quantiles.
#' @param a Shape parameter alpha.
#' @param b Shape parameter beta.
#' @return The sum of the log-density.
#' @export
beta_lpdf <- function(x, a, b) {
  sum(dbeta(x, a, b, log = TRUE))
}

#' Gamma log-probability density function
#'
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param rate Rate parameter.
#' @return The sum of the log-density.
#' @export
gamma_lpdf <- function(x, shape, rate) {
  sum(suppressWarnings(dgamma(x, shape = shape, rate = rate, log = TRUE)))
}

#' Inverse-gamma log-probability density function
#'
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @export
inverse_gamma_lpdf <- function(x, shape, scale) {
  sum(shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x)
}

#' Cauchy log-probability density function
#'
#' @param x Vector of quantiles.
#' @param location Location parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @export
cauchy_lpdf <- function(x, location, scale) {
  sum(suppressWarnings(dcauchy(x, location = location, scale = scale, log = TRUE)))
}

#' Student-t log-probability density function
#'
#' @param x Vector of quantiles.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#' @return The sum of the log-density.
#' @export
student_t_lpdf <- function(x, df, mu = 0, sigma = 1) {
  sum(suppressWarnings(dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma)))
}

#' Laplace log-probability density function
#'
#' @param x Vector of quantiles.
#' @param location Location parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @export
laplace_lpdf <- function(x, location = 0, scale = 1) {
  sum(-log(2 * scale) - abs(x - location) / scale)
}

#' Weibull log-probability density function
#'
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @export
weibull_lpdf <- function(x, shape, scale) {
  sum(suppressWarnings(dweibull(x, shape = shape, scale = scale, log = TRUE)))
}

#' Uniform log-probability density function
#'
#' @param x Vector of quantiles.
#' @param a Lower limit of the distribution.
#' @param b Upper limit of the distribution.
#' @return The sum of the log-density.
#' @export
uniform_lpdf <- function(x, a, b) {
  sum(dunif(x, a, b, log = TRUE))
}

#' LKJ correlation log-probability density function
#'
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

#' LKJ correlation log-probability density function for Cholesky factors
#'
#' @param CF_Omega Cholesky factor of a correlation matrix.
#' @param eta Shape parameter.
#' @return The log-density.
#' @export
lkj_CF_corr_lpdf <- function(CF_Omega, eta = 1) {
  if (eta == 1) return(0)
  return((eta - 1) * 2 * sum(log(diag(CF_Omega))))
}

#' Bernoulli log-probability mass function
#'
#' @param x Vector of binary outcomes (0 or 1).
#' @param prob Probability of success.
#' @return The sum of the log-probability.
#' @export
bernoulli_lpmf <- function(x, prob) {
  sum(dbinom(x, 1, prob, log = TRUE))
}

#' Binomial log-probability mass function
#'
#' @param x Vector of quantiles (number of successes).
#' @param size Number of trials.
#' @param prob Probability of success.
#' @return The sum of the log-probability.
#' @export
binomial_lpmf <- function(x, size, prob) {
  sum(dbinom(x, size, prob, log = TRUE))
}

#' Bernoulli log-probability mass function with logit parameterization
#'
#' @param x Vector of binary outcomes (0 or 1).
#' @param mu Linear predictor (logit probability).
#' @return The sum of the log-probability.
#' @export
bernoulli_logit_lpmf <- function(x, eta) {
  x * eta - log1p_exp(eta)
}

#' Binomial log-probability mass function with logit parameterization
#'
#' @param x Vector of quantiles (number of successes).
#' @param size Number of trials.
#' @param mu Linear predictor (logit probability).
#' @return The sum of the log-probability.
#' @export
binomial_logit_lpmf <- function(x, size, eta) {
  lgamma(size + 1) - lgamma(x + 1) - lgamma(size - x + 1) + x * eta - size * log1p_exp(eta)
}

#' Poisson log-probability mass function
#'
#' @param x Vector of quantiles.
#' @param mean Expected value.
#' @return The sum of the log-probability.
#' @export
poisson_lpmf <- function(x, mean) {
  sum(dpois(x, mean, log = TRUE))
}

#' Negative binomial log-probability mass function
#'
#' @param x Vector of quantiles.
#' @param size Target for number of successful trials.
#' @param prob Probability of success in each trial.
#' @return The sum of the log-probability.
#' @export
neg_binomial_lpmf <- function(x, size, prob) {
  sum(suppressWarnings(dnbinom(x, size = size, prob = prob, log = TRUE)))
}

#' Negative binomial log-probability mass function (alternative parameterization)
#'
#' @param x Vector of quantiles.
#' @param mu Mean parameter.
#' @param size Dispersion parameter.
#' @return The sum of the log-probability.
#' @export
neg_binomial_2_lpmf <- function(x, mu, size) {
  sum(suppressWarnings(dnbinom(x, size = size, mu = mu, log = TRUE)))
}

#' Categorical log-probability mass function
#'
#' @param x Vector of categorical outcomes.
#' @param prob Vector or matrix of probabilities.
#' @return The sum of the log-probability.
#' @export
categorical_lpmf <- function(x, prob) {
  if (is.matrix(prob)) {
    sum(log(prob[cbind(seq_along(x), x)]))
  } else {
    sum(log(prob[x]))
  }
}

#' Categorical log-probability mass function with logit parameterization
#'
#' @param x Integer or integer vector. The observed category index.
#' @param eta A numeric vector of linear predictors (unnormalized log-probabilities).
#' @return The sum of the log-probability.
#' @export
categorical_logit_lpmf <- function(x, eta) {
  sum(eta[x]) - length(x) * log_sum_exp(eta)
}

#' Multinomial log-probability mass function
#'
#' @param x Vector of counts.
#' @param size Total number of trials.
#' @param prob Vector of probabilities for each category.
#' @return The sum of the log-probability.
#' @export
multinomial_lpmf <- function(x, size, prob) {
  sum(suppressWarnings(dmultinom(x, size = size, prob = prob, log = TRUE)))
}

#' Ordered logistic log-probability mass function
#'
#' @param x Vector of ordered categorical outcomes.
#' @param eta Linear predictor.
#' @param cutpoints Vector of cutpoints.
#' @return The sum of the log-probability.
#' @export
ordered_logistic_lpmf <- function(x, eta, cutpoints) {
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

#' Dirichlet log-probability density function
#'
#' @param x Vector or matrix of simplexes.
#' @param alpha Vector of concentration parameters.
#' @return The sum of the log-density.
#' @export
dirichlet_lpdf <- function(x, alpha) {
  sum((alpha - 1) * log(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha))
}
#' Lower-triangular normal log-probability density function
#'
#' @param x A matrix of lower-triangular parameters.
#' @param mean Mean of the normal distribution.
#' @param sd Standard deviation of the normal distribution.
#' @return The log-density calculated only for the non-zero lower triangular elements.
#' @export
lower_tri_normal_lpdf <- function(x, mean = 0, sd = 1) {
  R <- nrow(x)
  C <- ncol(x)
  lp <- 0
  for (j in 1:min(R, C)) {
    val <- x[j:R, j]
    lp <- lp + sum(dnorm(val, mean, sd, log = TRUE))
  }
  return(lp)
}
#' Multivariate normal log-probability density function parameterized by Cholesky factor of correlation matrix
#'
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

#' Multivariate normal log-probability density function
#'
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

#' Sum-to-zero multivariate normal log-probability density function
#'
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

#' Centered triangular multivariate normal log-probability density function
#'
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
#' Positive centered triangular multivariate normal log-probability density function
#'
#' @param x Matrix of quantiles.
#' @param sigma Standard deviation parameter(s).
#' @return The log-density.
#' @export
positive_centered_tri_mvnormal_lpdf <- function(x, sigma = 1) {
  R <- nrow(x)
  C <- ncol(x)
  max_d <- min(C, R - 1)

  # 元の centered_tri_mvnormal_lpdf を計算
  lp <- centered_tri_mvnormal_lpdf(x, sigma)

  # 切断正規分布としての正規化定数を補正（各列で1/2になるため、密度を 2^max_d 倍する）
  lp <- lp + max_d * log(2)

  return(lp)
}

#' Best-Worst Scaling log-probability mass function
#'
#' @param x Integer. The index of the chosen best-worst pair (1 to C*(C-1)).
#' @param U A numeric vector of utilities for the items in the current choice set.
#' @param lambda A numeric scalar for the scale parameter (default is 1).
#' @return The log-probability of the observed best-worst choice.
#' @export
bw_categorical_logit_lpmf <- function(x, U, lambda = 1) {
  C <- length(U)
  ad_zero <- U[1] * 0
  U_dif <- rep(ad_zero, C * (C - 1))

  q_idx <- 1
  for (i in 1:C) {
    for (j in 1:C) {
      if (i != j) {
        U_dif[q_idx] <- U[i] - U[j]
        q_idx <- q_idx + 1
      }
    }
  }

  eta <- lambda * U_dif
  return(eta[x] - log_sum_exp(eta))
}
#' Wishart log-probability density function
#'
#' @param X Symmetric positive-definite matrix (scatter matrix).
#' @param n Degrees of freedom (sample size N).
#' @param V Scale matrix (covariance matrix Omega).
#' @return The log-density with all constant terms.
#' @export
wishart_lpdf <- function(X, n, V) {
  p <- nrow(X)
  L_V <- chol(V)
  log_det_V <- 2 * sum(log(diag(L_V)))
  trace_term <- sum(diag(solve(V, X)))
  log_det_X <- determinant(X, logarithm = TRUE)$modulus
  log_gamma_p <- (p * (p - 1) / 4) * log(pi) + sum(lgamma((n - (1:p) + 1) / 2))
  lp <- (n - p - 1) / 2 * log_det_X - 0.5 * trace_term -
    (n * p / 2) * log(2) - (n / 2) * log_det_V - log_gamma_p

  return(as.numeric(lp))
}
#' Factor analysis multivariate normal log-probability density function
#'
#' Woodbury matrix identity is used for efficient computation.
#'
#' @param x Vector or matrix of quantiles.
#' @param mu Vector of means.
#' @param Lambda Factor loading matrix (P x K).
#' @param psi Vector of unique variances (P).
#' @return The sum of the log-density.
#' @export
fa_multi_normal_lpdf <- function(x, mu, Lambda, psi) {
  P <- nrow(Lambda)
  K <- ncol(Lambda)
  inv_psi <- 1 / psi

  Lambda_scaled <- Lambda * inv_psi
  M <- diag(1, K) + (t(Lambda) %*% Lambda_scaled)
  L_M <- chol(M)
  log_det_Sigma <- sum(log(psi)) + 2 * sum(log(diag(L_M)))

  if (is.matrix(x)) {
    N <- nrow(x)
    y_c <- t(t(x) - mu)
    z_scaled <- t(t(y_c) * inv_psi)
    term1 <- sum(y_c * z_scaled)
    z_lambda <- z_scaled %*% Lambda
    theta_T <- solve(M, t(z_lambda))
    term2 <- sum(z_lambda * t(theta_T))

    lp <- -0.5 * (N * P * 1.83787706640935 + N * log_det_Sigma + term1 - term2)
    return(lp)

  } else {
    # ベクトル(1サンプル)の場合
    y_c <- x - mu
    term1 <- sum((y_c^2) * inv_psi)

    z_scaled <- y_c * inv_psi
    z_lambda <- as.vector(t(Lambda) %*% z_scaled)
    theta <- solve(M, z_lambda)
    term2 <- sum(z_lambda * theta)

    lp <- -0.5 * (P * 1.83787706640935 + log_det_Sigma + term1 - term2)
    return(lp)
  }
}
