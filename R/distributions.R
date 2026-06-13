#' Probability Distributions for RTMB Models
#'
#' @description
#' This page summarizes the probability distributions available within the `rtmb_code` block.
#' These functions are designed to be used with the Stan-like sampling syntax (`~`),
#' which internally adds the log-density to the model's total log-posterior.
#'
#' @details
#' \strong{Syntax Styles:}
#' In `rtmb_code`, you can specify distributions in two ways:
#' \itemize{
#'   \item \strong{Sampling Syntax (Recommended):} \code{y ~ normal(mu, sigma)}
#'   \item \strong{Explicit Function Call:} \code{lp <- lp + normal_lpdf(y, mu, sigma)}
#' }
#'
#' \strong{Continuous Distributions (LPDF):}
#' \itemize{
#'   \item \code{normal(mean, sd)}: Normal distribution.
#'   \item \code{lognormal(meanlog, sdlog)}: Lognormal distribution.
#'   \item \code{exponential(rate)}: Exponential distribution.
#'   \item \code{cauchy(location, scale)}: Cauchy distribution.
#'   \item \code{student_t(df, mu, sigma)}: Student's t-distribution.
#'   \item \code{gamma(shape, rate)}: Gamma distribution.
#'   \item \code{inverse_gamma(shape, scale)}: Inverse-gamma distribution.
#'   \item \code{beta(a, b)}: Beta distribution.
#' }
#'
#' \strong{Discrete Distributions (LPMF):}
#' \itemize{
#'   \item \code{bernoulli(prob)} / \code{bernoulli_logit(eta)}: Binary outcomes.
#'   \item \code{binomial(size, prob)} / \code{binomial_logit(size, eta)}: Binomial outcomes.
#'   \item \code{poisson(mean)}: Poisson count data.
#'   \item \code{neg_binomial_2(mu, size)}: Negative binomial (mean/dispersion parameterization).
#'   \item \code{ordered_logistic(eta, cutpoints)}: Ordered categorical outcomes.
#'   \item \code{sequential_logistic(eta, cutpoints)}: Sequential ordered categorical outcomes.
#' }
#'
#' \strong{Multivariate and Matrix Distributions:}
#' \itemize{
#'   \item \code{multi_normal(mean, Sigma)}: Standard multivariate normal distribution.
#'   \item \code{multi_student_t(df, mean, Sigma)}: Multivariate Student's t-distribution.
#'   \item \code{multi_cauchy(mean, Sigma)}: Multivariate Cauchy distribution (Student-t with df=1).
#'   \item \code{gaussian_process(x, mean, magnitude, smoothing, noise)}: Gaussian process prior with a squared exponential kernel.
#'   \item \code{lkj_corr(eta)}: LKJ prior for correlation matrices.
#'   \item \code{dirichlet(alpha)}: Dirichlet distribution for simplexes.
#'   \item \code{lower_tri_normal(mean, sd)}: Normal distribution for elements of a lower-triangular matrix.
#'   \item \code{centered_tri_multi_normal(sigma)}: Multivariate normal for centered triangular matrices (used in identification constraints).
#'   \item \code{sufficient_multi_normal_fa(S_mat, N, y_bar, mu, psi, Lambda)}: Factor analysis likelihood using sufficient statistics (highly efficient for large sample sizes).
#' }
#'
#' \strong{Vectorization:}
#' Most univariate distributions are vectorized. If \code{y} and \code{mu} are vectors,
#' \code{y ~ normal(mu, sigma)} will calculate the sum of log-densities for all elements
#' efficiently.
#'
#' @name distributions
#' @family distributions
#' @import RTMB
NULL

#' Normal log-probability density function
#'
#' @param x Vector of quantiles.
#' @param mean Vector of means.
#' @param sd Vector of standard deviations.
#' @param sum Logical; if TRUE (default), returns the sum of log-densities. If FALSE, returns element-wise log-densities.
#' @return The sum of the log-density (scalar) or a vector of log-densities.
#' @keywords internal
normal_lpdf <- function(x, mean, sd, sum = TRUE) {
  res <- dnorm(x, mean = mean, sd = sd, log = TRUE)
  if(sum) sum(res) else res
}

#' Lognormal log-probability density function
#'
#' @param x Vector of quantiles.
#' @param meanlog Mean of the distribution on the log scale.
#' @param sdlog Standard deviation of the distribution on the log scale.
#' @return The sum of the log-density.
#' @keywords internal
lognormal_lpdf <- function(x, meanlog, sdlog, sum = TRUE) {
  res <- suppressWarnings(dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE))
  if(sum) sum(res) else res
}

#' Exponential log-probability density function
#'
#' @param x Vector of quantiles.
#' @param rate Vector of rates.
#' @return The sum of the log-density.
#' @keywords internal
exponential_lpdf <- function(x, rate, sum = TRUE) {
  res <- dexp(x, rate, log = TRUE)
  if(sum) sum(res) else res
}

#' Half-Normal log-probability density function
#'
#' @param x Vector of quantiles (must be non-negative).
#' @param sd Vector of standard deviations (scale parameter).
#' @return The sum of the log-density.
#' @keywords internal
half_normal_lpdf <- function(x, sd, sum = TRUE) {
  res <- dnorm(x, mean = 0, sd = sd, log = TRUE) + log(2)
  if(sum) sum(res) else res
}

#' Beta log-probability density function
#'
#' @param x Vector of quantiles.
#' @param a Shape parameter alpha.
#' @param b Shape parameter beta.
#' @return The sum of the log-density.
#' @keywords internal
beta_lpdf <- function(x, a, b, sum = TRUE) {
  res <- dbeta(x, a, b, log = TRUE)
  if(sum) sum(res) else res
}

#' Gamma log-probability density function
#'
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param rate Rate parameter.
#' @return The sum of the log-density.
#' @keywords internal
gamma_lpdf <- function(x, shape, rate, sum = TRUE) {
  res <- suppressWarnings(dgamma(x, shape = shape, rate = rate, log = TRUE))
  if(sum) sum(res) else res
}

#' Inverse-gamma log-probability density function
#'
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @keywords internal
inverse_gamma_lpdf <- function(x, shape, scale, sum = TRUE) {
  res <- shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
  if(sum) sum(res) else res
}

#' Cauchy log-probability density function
#'
#' @param x Vector of quantiles.
#' @param location Location parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @keywords internal
cauchy_lpdf <- function(x, location, scale, sum = TRUE) {
  res <- suppressWarnings(dcauchy(x, location = location, scale = scale, log = TRUE))
  if(sum) sum(res) else res
}

#' Student-t log-probability density function
#'
#' @param x Vector of quantiles.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#' @return The sum of the log-density.
#' @keywords internal
student_t_lpdf <- function(x, df, mu = 0, sigma = 1, sum = TRUE) {
  res <- suppressWarnings(dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma))
  if(sum) sum(res) else res
}

#' Laplace log-probability density function
#'
#' @param x Vector of quantiles.
#' @param location Location parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @keywords internal
laplace_lpdf <- function(x, location = 0, scale = 1, sum = TRUE) {
  res <- -log(2 * scale) - sqrt((x - location)^2 + 1e-8) / scale
  if(sum) sum(res) else res
}

#' Logit-Normal log-probability density function
#'
#' @param x Vector of quantiles (must be strictly between 0 and 1).
#' @param mu Vector of means on the logit scale.
#' @param sd Vector of standard deviations on the logit scale.
#' @return The sum of the log-density.
#' @keywords internal
logit_normal_lpdf <- function(x, mu, sd, sum = TRUE) {
  logit_x <- qlogis(x)
  log_jacobian <- -log(x) - log1p(-x)
  res <- dnorm(logit_x, mean = mu, sd = sd, log = TRUE) + log_jacobian
  if(sum) sum(res) else res
}

#' Weibull log-probability density function
#'
#' @param x Vector of quantiles.
#' @param shape Shape parameter.
#' @param scale Scale parameter.
#' @return The sum of the log-density.
#' @keywords internal
weibull_lpdf <- function(x, shape, scale, sum = TRUE) {
  res <- suppressWarnings(dweibull(x, shape = shape, scale = scale, log = TRUE))
  if(sum) sum(res) else res
}

#' Uniform log-probability density function
#'
#' @param x Vector of quantiles.
#' @param a Lower limit of the distribution.
#' @param b Upper limit of the distribution.
#' @return The sum of the log-density.
#' @keywords internal
uniform_lpdf <- function(x, a, b, sum = TRUE) {
  res <- dunif(x, a, b, log = TRUE)
  if(sum) sum(res) else res
}

#' LKJ correlation log-probability density function
#'
#' @param Omega Correlation matrix.
#' @param eta Shape parameter.
#' @return The log-density.
#' @keywords internal
lkj_corr_lpdf <- function(Omega, eta = 1) {
  if (length(Omega) == 1) {
    rho <- Omega
    x <- (rho + 1) / 2
    return(dbeta(x, shape1 = eta, shape2 = eta, log = TRUE) - log(2))
  }
  # Matrix case (3 or more variables)
  K <- nrow(Omega)
  log_C <- 0
  for (k in 1:(K - 1)) {
    log_C <- log_C - (lgamma(eta + (K - 1 - k) / 2) +
                        lgamma(k / 2) -
                        lgamma(eta + (K - 1) / 2))
  }
  if (eta == 1) {
    return(log_C)
  }
  safe_Omega <- Omega + diag(1e-11, K)
  U <- chol(safe_Omega)
  log_det <- 2 * sum(log(diag(U)))
  return(log_C + (eta - 1) * log_det)
}

#' LKJ correlation log-probability density function for Cholesky factors
#'
#' @param CF_Omega Cholesky factor of a correlation matrix.
#' @param eta Shape parameter.
#' @return The log-density.
#' @keywords internal
lkj_CF_corr_lpdf <- function(CF_Omega, eta = 1) {
  K <- nrow(CF_Omega)
  log_C <- 0
  if (K > 1) {
    for (k in 1:(K - 1)) {
      log_C <- log_C - (lgamma(eta + (K - 1 - k) / 2) +
                          lgamma(k / 2) -
                          lgamma(eta + (K - 1) / 2))
    }
  }
  if (eta == 1) {
    return(log_C)
  }
  log_det <- 2 * sum(log(diag(CF_Omega)))
  return(log_C + (eta - 1) * log_det)
}

#' Bernoulli log-probability mass function
#'
#' @param x Vector of binary outcomes (0 or 1).
#' @param prob Probability of success.
#' @return The sum of the log-probability.
#' @keywords internal
bernoulli_lpmf <- function(x, prob, sum = TRUE) {
  res <- dbinom(x, 1, prob, log = TRUE)
  if(sum) sum(res) else res
}

#' Binomial log-probability mass function
#'
#' @param x Vector of quantiles (number of successes).
#' @param size Number of trials.
#' @param prob Probability of success.
#' @return The sum of the log-probability.
#' @keywords internal
binomial_lpmf <- function(x, size, prob, sum = TRUE) {
  res <- dbinom(x, size, prob, log = TRUE)
  if(sum) sum(res) else res
}

#' Bernoulli log-probability mass function with logit parameterization
#'
#' @param x Vector of binary outcomes (0 or 1).
#' @param eta Linear predictor (logit probability).
#' @return The sum of the log-probability.
#' @keywords internal
bernoulli_logit_lpmf <- function(x, eta, sum = TRUE) {
  res <- x * eta - log1p_exp(eta)
  if(sum) sum(res) else res
}

#' Binomial log-probability mass function with logit parameterization
#'
#' @param x Vector of quantiles (number of successes).
#' @param size Number of trials.
#' @param eta Linear predictor (logit probability).
#' @return The sum of the log-probability.
#' @keywords internal
binomial_logit_lpmf <- function(x, size, eta, sum = TRUE) {
  res <- lgamma(size + 1) - lgamma(x + 1) - lgamma(size - x + 1) + x * eta - size * log1p_exp(eta)
  if(sum) sum(res) else res
}

#' Poisson log-probability mass function
#'
#' @param x Vector of quantiles.
#' @param mean Expected value.
#' @return The sum of the log-probability.
#' @keywords internal
poisson_lpmf <- function(x, mean, sum = TRUE) {
  res <- dpois(x, mean, log = TRUE)
  if(sum) sum(res) else res
}

#' Negative binomial log-probability mass function
#'
#' @param x Vector of quantiles.
#' @param size Target for number of successful trials.
#' @param prob Probability of success in each trial.
#' @return The sum of the log-probability.
#' @keywords internal
neg_binomial_lpmf <- function(x, size, prob, sum = TRUE) {
  res <- suppressWarnings(dnbinom(x, size = size, prob = prob, log = TRUE))
  if(sum) sum(res) else res
}

#' Negative binomial log-probability mass function (alternative parameterization)
#'
#' @param x Vector of quantiles.
#' @param mu Mean parameter.
#' @param size Dispersion parameter.
#' @return The sum of the log-probability.
#' @keywords internal
neg_binomial_2_lpmf <- function(x, mu, size, sum = TRUE) {
  res <- lgamma(x + size) - lgamma(size) - lgamma(x + 1) +
    size * log(size) + x * log(mu) - (x + size) * log(size + mu)
  if(sum) sum(res) else res
}

#' Categorical log-probability mass function
#'
#' @param x Vector of categorical outcomes.
#' @param prob Vector or matrix of probabilities.
#' @return The sum of the log-probability.
#' @keywords internal
categorical_lpmf <- function(x, prob, sum = TRUE) {
  if (is.matrix(prob)) {
    res <- log(prob[cbind(seq_along(x), x)])
  } else {
    res <- log(prob[x])
  }
  if(sum) sum(res) else res
}

#' Categorical log-probability mass function with logit parameterization
#'
#' @param x Integer or integer vector. The observed category index.
#' @param eta A numeric vector of linear predictors (unnormalized log-probabilities).
#' @return The sum of the log-probability.
#' @keywords internal
categorical_logit_lpmf <- function(x, eta, sum = TRUE) {
  if (any(is.na(x))) {
    stop("categorical_logit() received NA category index.", call. = FALSE)
  }
  n_cat <- length(eta)
  if (any(x < 1 | x > n_cat)) {
    stop(
      sprintf(
        "categorical_logit() category index is outside the logit vector length: valid categories are 1..%d, but observed range is %d..%d. If this comes from matrix data, check whether you used x[t] instead of x[t, ].",
        n_cat,
        min(x),
        max(x)
      ),
      call. = FALSE
    )
  }
  res <- eta[x] - log_sum_exp(eta)
  if(sum) sum(res) else res
}

#' Multinomial log-probability mass function
#'
#' @param x Vector of counts.
#' @param size Total number of trials.
#' @param prob Vector of probabilities for each category.
#' @return The sum of the log-probability.
#' @keywords internal
multinomial_lpmf <- function(x, size, prob) {
  sum(suppressWarnings(dmultinom(x, size = size, prob = prob, log = TRUE)))
}

#' Ordered logistic log-probability mass function
#'
#' @param x Vector of ordered categorical outcomes.
#' @param eta Linear predictor.
#' @param cutpoints Vector of cutpoints.
#' @return The sum of the log-probability.
#' @keywords internal
ordered_logistic_lpmf <- function(x, eta, cutpoints, sum = TRUE) {
  N <- length(x)
  K <- length(cutpoints) + 1

  log1p_exp <- function(v) {
    max_val <- (v + sqrt(v^2 + 1e-11)) / 2
    return(max_val + log(exp(v - max_val) + exp(-max_val)))
  }
  log1m_exp <- function(v) {
    # In ordered-logistic middle categories, v = log(F_lower) - log(F_upper)
    # is non-positive by construction. Avoid pmin()/ifelse() here because
    # comparisons on RTMB AD types are unsafe.
    return(log(-expm1(v)))
  }

  if (length(eta) == 1 && N > 1) {
    eta_vec <- rep(eta, N)
  } else {
    eta_vec <- eta
  }

  # Perform vector operations using indices
  idx_1 <- which(x == 1)
  idx_K <- which(x == K)
  idx_mid <- which(x > 1 & x < K)

  if (sum) {
    lp <- 0
    if (length(idx_1) > 0) {
      lp <- lp - sum(log1p_exp(-(cutpoints[1] - eta_vec[idx_1])))
    }
    if (length(idx_K) > 0) {
      lp <- lp - sum(log1p_exp(cutpoints[K - 1] - eta_vec[idx_K]))
    }
    if (length(idx_mid) > 0) {
      for (ii in idx_mid) {
        y <- as.integer(x[ii])
        A <- -log1p_exp(-(cutpoints[y] - eta_vec[ii]))
        B <- -log1p_exp(-(cutpoints[y - 1L] - eta_vec[ii]))
        lp <- lp + A + log1m_exp(B - A)
      }
    }
    return(lp)
  } else {
    # Initialize with AD type safety
    res <- rep(cutpoints[1] * 0, N)
    if (length(idx_1) > 0) {
      res[idx_1] <- -log1p_exp(-(cutpoints[1] - eta_vec[idx_1]))
    }
    if (length(idx_K) > 0) {
      res[idx_K] <- -log1p_exp(cutpoints[K - 1] - eta_vec[idx_K])
    }
    if (length(idx_mid) > 0) {
      for (ii in idx_mid) {
        y <- as.integer(x[ii])
        A <- -log1p_exp(-(cutpoints[y] - eta_vec[ii]))
        B <- -log1p_exp(-(cutpoints[y - 1L] - eta_vec[ii]))
        res[ii] <- A + log1m_exp(B - A)
      }
    }
    return(res)
  }
}

#' Sequential logistic log-probability mass function
#'
#' @param x Vector of ordered categorical outcomes.
#' @param eta Linear predictor.
#' @param cutpoints Vector of sequential cutpoints.
#' @return The sum of the log-probability.
#' @keywords internal
sequential_logistic_lpmf <- function(x, eta, cutpoints, sum = TRUE) {
  N <- length(x)
  K <- length(cutpoints) + 1

  log1p_exp <- function(v) {
    max_val <- (v + sqrt(v^2 + 1e-11)) / 2
    max_val + log(exp(v - max_val) + exp(-max_val))
  }

  eta_is_matrix <- is.matrix(eta)
  if (!eta_is_matrix) {
    if (length(eta) == 1 && N > 1) {
      eta_vec <- rep(eta, N)
    } else {
      eta_vec <- eta
    }
  }

  log_prob_one <- function(i, y) {
    lp_i <- cutpoints[1] * 0
    if (y > 1L) {
      for (kk in 1:(y - 1L)) {
        eta_ik <- if (eta_is_matrix) eta[i, kk] else eta_vec[i]
        v <- eta_ik - cutpoints[kk]
        lp_i <- lp_i - log1p_exp(-v)
      }
    }
    if (y < K) {
      eta_iy <- if (eta_is_matrix) eta[i, y] else eta_vec[i]
      v <- eta_iy - cutpoints[y]
      lp_i <- lp_i - log1p_exp(v)
    }
    lp_i
  }

  if (sum) {
    lp <- cutpoints[1] * 0
    for (ii in seq_len(N)) {
      lp <- lp + log_prob_one(ii, as.integer(x[ii]))
    }
    lp
  } else {
    res <- rep(cutpoints[1] * 0, N)
    for (ii in seq_len(N)) {
      res[ii] <- log_prob_one(ii, as.integer(x[ii]))
    }
    res
  }
}

#' Dirichlet log-probability density function
#'
#' @param x Vector or matrix of simplexes.
#' @param alpha Vector of concentration parameters.
#' @return The sum of the log-density.
#' @keywords internal
dirichlet_lpdf <- function(x, alpha) {
  sum((alpha - 1) * log(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha))
}
#' Lower-triangular normal log-probability density function
#'
#' @param x A matrix of lower-triangular parameters.
#' @param mean Mean of the normal distribution.
#' @param sd Standard deviation of the normal distribution.
#' @return The log-density calculated only for the non-zero lower triangular elements.
#' @keywords internal
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
#' Positive lower-triangular normal log-probability density function
#'
#' @param x A matrix of lower-triangular parameters (Cholesky factor with positive diagonals).
#' @param mean Mean of the normal distribution (assumed to be 0 for the half-normal correction).
#' @param sd Standard deviation of the normal distribution.
#' @return The log-density calculated for the non-zero lower triangular elements.
#' @keywords internal
positive_lower_tri_normal_lpdf <- function(x, mean = 0, sd = 1) {
  R <- nrow(x)
  C <- ncol(x)
  lp <- 0

  for (j in 1:min(R, C)) {
    val <- x[j:R, j]
    lp <- lp + sum(dnorm(val, mean, sd, log = TRUE))
  }

  num_diag <- min(R, C)
  lp <- lp + num_diag * log(2)

  return(lp)
}
#' Multivariate normal log-probability density function parameterized by Cholesky factor of correlation matrix
#'
#' @param x Vector or matrix of quantiles.
#' @param mean Vector or matrix of means.
#' @param sd Vector of standard deviations.
#' @param CF_Omega Cholesky factor of the correlation matrix.
#' @return The sum of the log-density.
#' @keywords internal
multi_normal_CF_lpdf <- function(x, mean, sd, CF_Omega, sum = TRUE) {
  if (length(sd) == 1) {
    L_Sigma <- matrix(sd, 1, 1) %*% CF_Omega
  } else {
    L_Sigma <- diag(sd) %*% CF_Omega
  }

  N <- if (is.matrix(x)) nrow(x) else 1
  K <- ncol(L_Sigma)

  log_det <- 2 * sum(log(diag(L_Sigma)))
  const <- -0.5 * (K * log(2 * pi) + log_det)

  if (is.matrix(x)) {
    resid_t <- if (is.matrix(mean)) t(x - mean) else t(x) - mean
    z <- solve(L_Sigma, resid_t)
    if (sum) {
      quad_form <- sum(z^2)
      return(N * const - 0.5 * quad_form)
    } else {
      quad_per_obs <- colSums(z^2)
      return(const - 0.5 * quad_per_obs)
    }
  } else {
    z <- solve(L_Sigma, x - mean)
    quad_form <- sum(z^2)
    return(const - 0.5 * quad_form)
  }
}

#' Multivariate normal log-probability density function parameterized by Cholesky factor of correlation matrix (FIML)
#'
#' @param x Vector or matrix of quantiles containing NAs.
#' @param mean Vector of means.
#' @param sd Vector of standard deviations.
#' @param CF_Omega Cholesky factor of the correlation matrix.
#' @return The sum of the log-density using Full Information Maximum Likelihood.
#' @keywords internal
fiml_multi_normal_CF_lpdf <- function(x, mean, sd, CF_Omega) {
  if (length(sd) == 1) {
    L_Sigma <- matrix(sd, 1, 1) %*% CF_Omega
  } else {
    L_Sigma <- diag(sd) %*% CF_Omega
  }
  Sigma <- L_Sigma %*% t(L_Sigma)

  N <- nrow(x)
  lp <- 0
  
  # Ensure positive definiteness with a scale-aware jitter
  diag_Sigma <- diag(Sigma)
  safe_Sigma <- Sigma + diag(diag_Sigma * 1e-6 + 1e-8)

  for (i in 1:N) {
    y_i <- x[i, ]
    keep <- !is.na(y_i)
    if (any(keep)) {
      y_obs <- y_i[keep]
      mean_obs <- mean[keep]
      Sigma_obs <- safe_Sigma[keep, keep, drop = FALSE]

      log_det_obj <- determinant(Sigma_obs, logarithm = TRUE)
      log_det <- log_det_obj$modulus
      k <- sum(keep)
      const <- -0.5 * (k * log(2 * pi) + log_det)

      resid <- y_obs - mean_obs
      quad_form <- sum(resid * solve(Sigma_obs, resid))
      lp <- lp + const - 0.5 * quad_form
    }
  }
  return(lp)
}

#' Multivariate normal log-probability density function
#'
#' @param x Vector or matrix of quantiles.
#' @param mean Vector or matrix of means.
#' @param Sigma Covariance matrix.
#' @return The sum of the log-density.
#' @keywords internal
multi_normal_lpdf <- function(x, mean, Sigma, sum = TRUE) {

  log_det_chol <- function(L) {
    return(2 * sum(log(diag(L))))
  }
  quad_form_chol <- function(x, L) {
    z <- solve(L, x)
    return(sum(z^2))
  }

  K <- nrow(Sigma)
  if (is.null(K) || length(Sigma) == 1) {
    res <- dnorm(x, mean = mean, sd = sqrt(Sigma + 1e-11), log = TRUE)
    return(if (sum) sum(res) else res)
  }
  # Ensure positive definiteness with a scale-aware jitter
  # This prevents singularity even if the elements of Sigma are very large or very small.
  diag_Sigma <- diag(Sigma)
  safe_Sigma <- Sigma + diag(diag_Sigma * 1e-6 + 1e-8)

  # Calculate log-determinant and quadratic form more robustly
  # For RTMB, using determinant() and solve() can sometimes be more stable than chol()
  # if the matrix is very near-singular.
  
  log_det_obj <- determinant(safe_Sigma, logarithm = TRUE)
  log_det <- log_det_obj$modulus
  
  const <- -0.5 * (K * log(2 * pi) + log_det)

  if (is.matrix(x) && ncol(x) > 1) {
    n <- nrow(x)
    resid_t <- if (is.matrix(mean)) t(x - mean) else t(x) - mean
    quad_form <- colSums(resid_t * solve(safe_Sigma, resid_t))
    res <- const - 0.5 * quad_form
    return(if (sum) sum(res) else res)
  } else {
    # Single observation
    x_vec <- if(is.matrix(x)) drop(x) else x
    mean_vec <- if(is.matrix(mean)) drop(mean) else mean
    resid <- x_vec - mean_vec
    quad_form <- sum(resid * solve(safe_Sigma, resid))
    return(const - 0.5 * quad_form)
  }
}

#' Multivariate Student-t log-probability density function
#'
#' @param x Vector or matrix of quantiles.
#' @param df Degrees of freedom.
#' @param mean Vector or matrix of means.
#' @param Sigma Scale matrix.
#' @return The sum of the log-density.
#' @keywords internal
multi_student_t_lpdf <- function(x, df, mean, Sigma, sum = TRUE) {
  # Robust dimension check for RTMB
  Sigma_val <- drop(Sigma)
  K <- if (is.matrix(Sigma)) nrow(Sigma) else length(Sigma_val)
  
  if (length(Sigma_val) == 1) {
    res <- student_t_lpdf(x, df = df, mu = mean, sigma = sqrt(Sigma_val + 1e-11), sum = FALSE)
    return(if (sum) sum(res) else res)
  }
  
  # Ensure positive definiteness with a scale-aware jitter
  diag_Sigma <- diag(Sigma)
  safe_Sigma <- Sigma + diag(diag_Sigma * 1e-6 + 1e-8)
  
  # Calculate log-determinant and quadratic form more robustly
  log_det_obj <- determinant(safe_Sigma, logarithm = TRUE)
  log_det <- log_det_obj$modulus
  
  const <- lgamma((df + K) / 2) - lgamma(df / 2) - 0.5 * (K * log(df * pi) + log_det)
  
  if (is.matrix(x) && ncol(x) > 1) {
    n <- nrow(x)
    resid_t <- if (is.matrix(mean)) t(x - mean) else t(x) - mean
    quad_form <- colSums(resid_t * solve(safe_Sigma, resid_t))
    res <- const - 0.5 * (df + K) * log1p(quad_form / df)
    return(if (sum) sum(res) else res)
  } else {
    # Single observation (vector or column-matrix)
    x_vec <- if(is.matrix(x)) drop(x) else x
    mean_vec <- if(is.matrix(mean)) drop(mean) else mean
    resid <- x_vec - mean_vec
    quad_form <- sum(resid * solve(safe_Sigma, resid))
    return(const - 0.5 * (df + K) * log1p(quad_form / df))
  }
}

#' Multivariate Cauchy log-probability density function
#'
#' @param x Vector or matrix of quantiles.
#' @param mean Vector or matrix of means.
#' @param Sigma Scale matrix.
#' @return The sum of the log-density.
#' @keywords internal
multi_cauchy_lpdf <- function(x, mean, Sigma, sum = TRUE) {
  multi_student_t_lpdf(x, df = 1, mean = mean, Sigma = Sigma, sum = sum)
}

#' Normal Mixture log-probability density function
#'
#' @param x Vector of quantiles (observations).
#' @param pi_w Vector of mixing proportions for each cluster.
#' @param mean Vector of means for each cluster.
#' @param sd Vector of standard deviations for each cluster.
#' @return The sum of the log-density.
#' @keywords internal
normal_mixture_lpdf <- function(x, pi_w, mean, sd, sum = TRUE) {
  log_pi <- log(pi_w)
  K <- length(pi_w)
  N <- length(x)

  # Fully vectorized version using matrix operations
  # log_dens: N x K matrix
  # Initialize with AD-safe values (using 0 based on log_pi type)
  log_dens <- matrix(log_pi[1] * 0, N, K)
  for (k in 1:K) {
    log_dens[, k] <- dnorm(x, mean = mean[k], sd = sd[k], log = TRUE) + log_pi[k]
  }

  # Log-Sum-Exp row-wise
  max_lp <- apply(log_dens, 1, max)
  res <- max_lp + log(rowSums(exp(log_dens - max_lp)))

  if(sum) sum(res) else res
}

#' Generic Mixture log-probability density function
#'
#' @param x Vector of quantiles.
#' @param pi_w Vector or matrix of mixing proportions.
#' @param lpdf_list A list of log-density vectors (one for each cluster).
#' @return The sum of the log-density.
#' @keywords internal
mixture_lpdf <- function(x, pi_w, lpdf_list, sum = TRUE) {
  N <- if (is.matrix(x)) nrow(x) else length(x)
  K <- length(lpdf_list)

  # log_pi: N x K matrix
  if (is.matrix(pi_w)) {
    log_pi <- log(pi_w)
  } else {
    log_pi <- matrix(log(pi_w), N, K, byrow = TRUE)
  }

  # log_dens: N x K matrix
  log_dens <- matrix(log_pi[1] * 0, N, K)
  for (k in 1:K) {
    log_dens[, k] <- lpdf_list[[k]] + log_pi[, k]
  }

  # Log-Sum-Exp row-wise
  max_lp <- apply(log_dens, 1, max)
  res <- max_lp + log(rowSums(exp(log_dens - max_lp)))

  if(sum) sum(res) else res
}


#' Centered / Centered matrix multivariate normal log-probability density function
#'
#' @param x Vector or matrix of quantiles. If a matrix, it assumes each column sums to zero.
#' @param sigma Standard deviation parameter. Can be a scalar, or a vector of length ncol(x) if x is a matrix.
#' @param K Dimension. If NULL, inferred from the length or nrow of x.
#' @return The log-density.
#' @keywords internal
centered_multi_normal_lpdf <- function(x, sigma, K = NULL) {
  if (is.matrix(x)) {
    R <- nrow(x)
    C <- ncol(x)
    if (length(sigma) == 1) {
      df <- (R - 1) * C
      lp <- -0.5 * df * log(2 * pi) - df * log(sigma) - 0.5 * sum(x^2) / sigma^2
    } else {
      lp <- 0
      for (c in 1:C) {
        x_c <- x[, c]
        sd_c <- sigma[c]
        lp <- lp - 0.5 * (R - 1) * log(2 * pi) - (R - 1) * log(sd_c) - 0.5 * sum(x_c^2) / sd_c^2
      }
    }
    return(lp)
  } else {
    if (is.null(K)) K <- length(x)
    lp <- -0.5 * (K - 1) * log(2 * pi) - (K - 1) * log(sigma) - 0.5 * sum(x^2) / sigma^2
    return(lp)
  }
}

#' Centered triangular multivariate normal log-probability density function
#'
#' @param x Matrix of quantiles.
#' @param sigma Standard deviation parameter(s).
#' @return The log-density.
#' @keywords internal
centered_tri_multi_normal_lpdf <- function(x, sigma) {
  R <- nrow(x)
  C <- ncol(x)
  max_d <- min(C, R - 1)

  lp <- 0
  if (length(sigma) == 1) {
    df <- sum(R - (1:max_d))
    lp <- -0.5 * df * log(2 * pi) - df * log(sigma) - 0.5 * sum(x^2) / sigma^2
  } else {
    # When specifying a different sigma for each dimension
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
#' @param x Matrix of quantiles (centered triangular matrix with positive diagonals).
#' @param sigma Standard deviation parameter(s).
#' @return The log-density.
#' @keywords internal
positive_centered_tri_multi_normal_lpdf <- function(x, sigma = 1) {
  R <- nrow(x)
  C <- ncol(x)
  max_d <- min(C, R - 1)

  # Calculate the original centered_tri_multi_normal_lpdf
  lp <- centered_tri_multi_normal_lpdf(x, sigma)

  # Correct the normalizing constant for the truncated normal distribution (multiply density by 2^max_d since it's halved for each column)
  lp <- lp + max_d * log(2)

  return(lp)
}

#' Best-Worst Categorical Logit Log-PMF (RTMB optimized)
#'
#' @param x Vector of indices of the selected pairs.
#' @param U Vector of utilities for each item.
#' @param lambda Scaling parameter.
#' @return Log-probability (scalar advector).
#' @keywords internal
bw_categorical_logit_lpmf <- function(x, U, lambda = 1) {
  # 1. Handle missing values
  if (any(is.na(as.numeric(x)))) return(0)

  C <- length(U)

  # 2. Generate all pairs (i, j) where i != j using vectorized operations
  # Instead of loops, we use rep() to create index vectors
  idx_i <- rep(1:C, each = C)
  idx_j <- rep(1:C, times = C)
  mask <- idx_i != idx_j # Logical vector to exclude i == j

  # 3. Compute utility differences for all pairs at once
  # This avoids partial assignment and keeps the AD tape clean
  U_dif <- U[idx_i[mask]] - U[idx_j[mask]]

  # 4. Compute log-probabilities
  # eta is a vector of length C*(C-1)
  eta <- lambda * U_dif

  # Log-sum-exp is computed once per set
  # sum(eta[x]) handles multiple observations if x is a vector
  log_p_obs <- sum(eta[x])
  log_p_denom <- length(x) * log_sum_exp(eta)

  return(log_p_obs - log_p_denom)
}
#' Wishart log-probability density function
#'
#' @param X Symmetric positive-definite matrix (scatter matrix).
#' @param n Degrees of freedom (sample size N).
#' @param V Scale matrix (covariance matrix Omega).
#' @return The log-density with all constant terms.
#' @keywords internal
wishart_lpdf <- function(X, n, V) {
  p <- nrow(X)
  L_V <- chol(V)
  log_det_V <- 2 * sum(log(diag(L_V)))
  trace_term <- sum(diag(solve(V, X)))
  log_det_X <- determinant(X, logarithm = TRUE)$modulus
  log_gamma_p <- (p * (p - 1) / 4) * log(pi) + sum(lgamma((n - (1:p) + 1) / 2))
  lp <- (n - p - 1) / 2 * log_det_X - 0.5 * trace_term -
    (n * p / 2) * log(2) - (n / 2) * log_det_V - log_gamma_p

  return(lp)
}
#' Wishart log-probability density function parameterized by Cholesky factor of correlation matrix
#'
#' @param X Symmetric positive-definite matrix (scatter matrix).
#' @param n Degrees of freedom (sample size N).
#' @param sd Vector of standard deviations.
#' @param CF_Omega Cholesky factor of the correlation matrix.
#' @return The log-density with all constant terms.
#' @keywords internal
wishart_CF_lpdf <- function(X, n, sd, CF_Omega) {
  p <- nrow(X)

  if (length(sd) == 1) {
    L_Sigma <- matrix(sd, 1, 1) %*% CF_Omega
  } else {
    L_Sigma <- diag(sd) %*% CF_Omega
  }

  log_det_Sigma <- 2 * sum(log(diag(L_Sigma)))
  L_inv <- solve(L_Sigma)
  Sigma_inv <- t(L_inv) %*% L_inv
  trace_term <- sum(Sigma_inv * X)

  log_det_X <- determinant(X, logarithm = TRUE)$modulus
  log_gamma_p <- (p * (p - 1) / 4) * log(pi) + sum(lgamma((n - (1:p) + 1) / 2))

  lp <- (n - p - 1) / 2 * log_det_X - 0.5 * trace_term -
    (n * p / 2) * log(2) - (n / 2) * log_det_Sigma - log_gamma_p

  return(lp)
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
#' @keywords internal
fa_multi_normal_lpdf <- function(x, mu, Lambda, psi, sum = TRUE) {
  P <- nrow(Lambda)
  K <- ncol(Lambda)
  inv_psi <- 1 / psi

  Lambda_scaled <- Lambda * inv_psi
  M <- diag(1, K) + (t(Lambda) %*% Lambda_scaled)
  L_M <- chol(M)
  log_det_Sigma <- sum(log(psi)) + 2 * sum(log(diag(L_M)))
  const <- -0.5 * (P * 1.83787706640935 + log_det_Sigma)

  if (is.matrix(x)) {
    N <- nrow(x)
    y_c <- t(t(x) - mu)
    z_scaled <- t(t(y_c) * inv_psi)
    z_lambda <- z_scaled %*% Lambda
    theta_T <- solve(M, t(z_lambda))

    if (sum) {
      term1 <- sum(y_c * z_scaled)
      term2 <- sum(z_lambda * t(theta_T))
      return(N * const - 0.5 * (term1 - term2))
    } else {
      term1_per_obs <- rowSums(y_c * z_scaled)
      term2_per_obs <- rowSums(z_lambda * t(theta_T))
      return(const - 0.5 * (term1_per_obs - term2_per_obs))
    }

  } else {
    # Vector case (single sample)
    y_c <- x - mu
    term1 <- sum((y_c^2) * inv_psi)

    z_scaled <- y_c * inv_psi
    z_lambda <- as.vector(t(Lambda) %*% z_scaled)
    theta <- solve(M, z_lambda)
    term2 <- sum(z_lambda * theta)

    return(const - 0.5 * (term1 - term2))
  }
}
#' Factor analysis multivariate normal log-probability density function (FIML)
#'
#' @param x Matrix of quantiles containing NAs.
#' @param mu Vector of means.
#' @param Lambda Factor loading matrix (P x K).
#' @param psi Vector of unique variances (P).
#' @return The sum of the log-density.
#' @keywords internal
fiml_multi_normal_fa_lpdf <- function(x, mu, Lambda, psi) {
  P <- nrow(Lambda)
  K <- ncol(Lambda)
  
  # Implied Covariance: Sigma = Lambda %*% t(Lambda) + diag(psi^2)
  # But for FIML, we can either construct Sigma and invert subsets, or use Woodbury on subsets.
  # Woodbury on subsets is complex because M = I + Lambda_obs' * inv_psi_obs * Lambda_obs changes per missing pattern.
  # Given typical RTMB usage, constructing Sigma and subsetting is simplest.
  Sigma <- Lambda %*% t(Lambda)
  diag(Sigma) <- diag(Sigma) + psi^2
  
  N <- nrow(x)
  lp <- 0
  
  diag_Sigma <- diag(Sigma)
  safe_Sigma <- Sigma + diag(diag_Sigma * 1e-6 + 1e-8)

  for (i in 1:N) {
    y_i <- x[i, ]
    keep <- !is.na(y_i)
    if (any(keep)) {
      y_obs <- y_i[keep]
      mean_obs <- mu[keep]
      Sigma_obs <- safe_Sigma[keep, keep, drop = FALSE]

      log_det_obj <- determinant(Sigma_obs, logarithm = TRUE)
      log_det <- log_det_obj$modulus
      k <- sum(keep)
      const <- -0.5 * (k * log(2 * pi) + log_det)

      resid <- y_obs - mean_obs
      quad_form <- sum(resid * solve(Sigma_obs, resid))
      lp <- lp + const - 0.5 * quad_form
    }
  }
  return(lp)
}

#' Sufficient statistics multivariate normal log-probability density function
#'
#' @param S_mat Deviation sum of squares matrix.
#' @param N Sample size.
#' @param y_bar Sample mean vector.
#' @param mean Mean parameter vector.
#' @param sd Standard deviation parameter vector.
#' @param CF_Omega Cholesky factor of correlation matrix.
#' @return The exact log-likelihood of the N raw observations.
#' @keywords internal
sufficient_multi_normal_CF_lpdf <- function(S_mat, N, y_bar, mean, sd, CF_Omega) {
  p <- length(y_bar)

  # 1. Calculate L_Sigma and log|Sigma|
  if (length(sd) == 1) {
    L_Sigma <- matrix(sd, 1, 1) %*% CF_Omega
  } else {
    L_Sigma <- diag(sd) %*% CF_Omega
  }
  log_det_Sigma <- 2 * sum(log(diag(L_Sigma)))

  # 2. Mahalanobis distance of the mean vector (N * (y_bar - mu)^T Sigma^-1 (y_bar - mu))
  z <- solve(L_Sigma, y_bar - mean)
  quad_mean <- sum(z^2)

  # 3. Trace term of the deviation sum of squares matrix (tr(Sigma^-1 S_mat))
  L_inv <- solve(L_Sigma)
  Sigma_inv <- t(L_inv) %*% L_inv
  trace_term <- sum(Sigma_inv * S_mat)

  # 4. Full log-likelihood of the original multivariate normal distribution
  lp <- -0.5 * N * p * log(2 * pi) -
    0.5 * N * log_det_Sigma -
    0.5 * trace_term -
    0.5 * N * quad_mean

  return(lp)
}

#' Sufficient statistics multivariate normal log-probability density function (covariance parameterization)
#'
#' @param S_mat Deviation sum of squares matrix.
#' @param N Sample size.
#' @param y_bar Sample mean vector.
#' @param mean Mean parameter vector.
#' @param Sigma Covariance matrix.
#' @return The exact log-likelihood of the N raw observations.
#' @keywords internal
sufficient_multi_normal_lpdf <- function(S_mat, N, y_bar, mean, Sigma) {
  p <- length(y_bar)

  # 1. Calculate log|Sigma|
  # Ensure positive definiteness with a scale-aware jitter
  diag_Sigma <- diag(Sigma)
  safe_Sigma <- Sigma + diag(diag_Sigma * 1e-6 + 1e-8)

  log_det_obj <- determinant(safe_Sigma, logarithm = TRUE)
  log_det_Sigma <- log_det_obj$modulus

  # 2. Mahalanobis distance of the mean vector
  resid <- y_bar - mean
  Sigma_inv <- solve(safe_Sigma)
  quad_mean <- sum(resid * (Sigma_inv %*% resid))

  # 3. Trace term of the deviation sum of squares matrix (tr(Sigma^-1 S_mat))
  trace_term <- sum(Sigma_inv * S_mat)

  # 4. Full log-likelihood
  lp <- -0.5 * N * p * log(2 * pi) -
    0.5 * N * log_det_Sigma -
    0.5 * trace_term -
    0.5 * N * quad_mean

  return(lp)
}

#' Sufficient statistics factor analysis multivariate normal log-probability density function
#'
#' @param S_mat Deviation sum of squares matrix.
#' @param N Sample size.
#' @param y_bar Sample mean vector.
#' @param mu Mean parameter vector.
#' @param Lambda Factor loading matrix (P x K).
#' @param psi Vector of unique variances (P).
#' @return The exact log-likelihood of the N raw observations.
#' @keywords internal
sufficient_multi_normal_fa_lpdf <- function(S_mat, N, y_bar, mu, psi,Lambda) {
  P <- nrow(Lambda)
  K <- ncol(Lambda)

  psi_safe <- psi^2 + 1e-8
  inv_psi <- 1 / psi_safe

  # Modification 1: Avoid recycling multiplication with vectors and use explicit matrix multiplication
  Lambda_scaled <- diag(inv_psi) %*% Lambda

  M <- diag(1, K) + (t(Lambda) %*% Lambda_scaled)
  L_M <- chol(M)

  log_det_Sigma <- sum(log(psi_safe)) + 2 * sum(log(diag(L_M)))

  # --- 1. Calculate trace term ---
  term1_trace <- sum(diag(S_mat) * inv_psi)

  Q <- t(Lambda_scaled) %*% S_mat %*% Lambda_scaled
  term2_trace <- sum(diag(solve(M, Q)))
  trace_term <- term1_trace - term2_trace

  # --- 2. Calculate mean term ---
  d <- y_bar - mu
  term1_mean <- sum((d^2) * inv_psi)

  # Modification 2: Avoid as.vector() and treat strictly as a column vector (matrix)
  d_mat <- matrix(d, ncol = 1)
  z_lambda <- t(Lambda_scaled) %*% d_mat
  theta <- solve(M, z_lambda)

  # Modification 3: Extract scalar from inner product of matrices
  term2_mean <- sum(t(z_lambda) %*% theta)

  quad_mean <- term1_mean - term2_mean

  # --- 3. Calculate full log-likelihood ---
  lp <- -0.5 * (N * P * 1.83787706640935 +
                  N * log_det_Sigma +
                  trace_term +
                  N * quad_mean)

  return(sum(lp))
}

#' Gaussian Process Log-Density (Squared Exponential Kernel)
#'
#' @description
#' Calculates the log-density of a Gaussian Process with a Squared Exponential (RBF) kernel.
#'
#' @param y Observation vector (N), or matrix (M x N) whose rows are independent GP realizations on the same coordinates.
#' @param x Coordinate vector or matrix (N x D), where N matches \code{length(y)} for vector \code{y} or \code{ncol(y)} for matrix \code{y}.
#' @param mean Mean vector (scalar, length N, or M x N matrix for matrix \code{y}).
#' @param magnitude Signal standard deviation (alpha).
#' @param smoothing Length-scale (rho).
#' @param noise Measurement noise standard deviation (sigma).
#' @param sum Logical; whether to return the sum of log-densities. If \code{FALSE} and \code{y} is a matrix, returns one log-density per row.
#' @return Log-density value.
#' @export
gaussian_process_lpdf <- function(y, x, mean = 0, magnitude = 1, smoothing = 1, noise = 0.01, sum = TRUE) {
  y_is_matrix <- is.matrix(y)
  N <- if (y_is_matrix) ncol(y) else length(y)
  X <- if (is.matrix(x)) x else as.matrix(x)
  if (nrow(X) != N) {
    stop(
      "gaussian_process_lpdf(): number of coordinate rows must match ",
      if (y_is_matrix) "ncol(y)" else "length(y)",
      ".",
      call. = FALSE
    )
  }

  # Construct covariance matrix K
  K <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in i:N) {
      # Squared Euclidean distance
      d2 <- sum((X[i, ] - X[j, ])^2)

      # Squared Exponential kernel
      val <- magnitude^2 * exp(-d2 / (2 * smoothing^2))
      K[i, j] <- val
      if (i != j) K[j, i] <- val # Symmetry
    }
  }

  # Add measurement noise and a small jitter for numerical stability
  diag(K) <- diag(K) + noise^2 + 1e-8

  # Evaluate as Multivariate Normal
  if (y_is_matrix && N == 1L) {
    mean_vec <- if (is.matrix(mean)) mean[, 1] else mean
    return(normal_lpdf(y[, 1], mean = mean_vec, sd = sqrt(K[1, 1]), sum = sum))
  }
  return(multi_normal_lpdf(y, mean = mean, Sigma = K, sum = sum))
}
