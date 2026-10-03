#' Log-CDF and Log-CCDF Functions
#'
#' @description
#' Log cumulative distribution functions and log complementary cumulative
#' distribution functions for univariate distributions used in
#' `rtmb_code()`. These functions are useful for censored and interval-censored
#' likelihoods and are compatible with RTMB automatic differentiation.
#'
#' @details
#' Each function returns the sum of the element-wise log probabilities by
#' default. Set `sum = FALSE` to return one value per observation. The naming
#' follows Stan: `_lcdf` computes `log(P[X <= x])`, and `_lccdf` computes
#' `log(P[X > x])`. For discrete distributions this distinction is important.
#'
#' Supported continuous distributions are normal, lognormal, exponentially
#' modified normal, exponential, half-normal, beta, gamma, inverse-gamma,
#' Cauchy, Student-t, Laplace, logit-normal, Weibull, and uniform, as well as
#' normal and generic mixtures. Supported discrete distributions are Bernoulli,
#' binomial,
#' Bernoulli-logit, binomial-logit, Poisson, negative binomial,
#' negative-binomial-2, ordered logistic, and sequential logistic.
#'
#' Multivariate and matrix distributions, nominal categorical distributions,
#' and diffusion models do not have a single generally applicable scalar CDF
#' and are therefore not included.
#'
#' @param x Quantile or vector of quantiles.
#' @param mean Mean parameter. For Poisson functions, the rate parameter.
#' @param sd Positive standard deviation.
#' @param meanlog,sdlog Mean and standard deviation on the log scale.
#' @param mu Location or mean parameter.
#' @param sigma Positive scale parameter.
#' @param lambda Positive exponential rate.
#' @param rate Positive rate parameter.
#' @param a,b Positive shape parameters, or lower and upper bounds for the
#'   uniform distribution.
#' @param shape Positive shape parameter.
#' @param scale Positive scale parameter.
#' @param location Location parameter.
#' @param df Positive degrees of freedom.
#' @param pi_w Mixture weights that sum to one.
#' @param lcdf_list,lccdf_list Lists of component-wise log-CDF or log-CCDF
#'   vectors, respectively.
#' @param size Number of trials or positive dispersion parameter.
#' @param prob Probability parameter.
#' @param eta Linear predictor.
#' @param cutpoints Ordered cutpoints.
#' @param sum Logical; return the sum when `TRUE`, otherwise return
#'   element-wise values.
#' @return A scalar sum of log probabilities when `sum = TRUE`, or a vector
#'   of element-wise log probabilities when `sum = FALSE`.
#' @name distribution_cdf
#' @family distributions
#' @keywords internal
NULL

.rtmb_finish_log_cdf <- function(res, sum = TRUE) {
  if (sum) base::sum(res) else res
}

#' @rdname distribution_cdf
normal_lcdf <- function(x, mean, sd, sum = TRUE) {
  .rtmb_finish_log_cdf(
    RTMB::pnorm(x, mean = mean, sd = sd, lower.tail = TRUE, log.p = TRUE),
    sum
  )
}

#' @rdname distribution_cdf
normal_lccdf <- function(x, mean, sd, sum = TRUE) {
  .rtmb_finish_log_cdf(
    RTMB::pnorm(x, mean = mean, sd = sd, lower.tail = FALSE, log.p = TRUE),
    sum
  )
}

#' @rdname distribution_cdf
lognormal_lcdf <- function(x, meanlog, sdlog, sum = TRUE) {
  x <- 0.5 * (x + abs(x))
  .rtmb_finish_log_cdf(
    RTMB::pnorm(log(x), mean = meanlog, sd = sdlog,
                lower.tail = TRUE, log.p = TRUE),
    sum
  )
}

#' @rdname distribution_cdf
lognormal_lccdf <- function(x, meanlog, sdlog, sum = TRUE) {
  x <- 0.5 * (x + abs(x))
  .rtmb_finish_log_cdf(
    RTMB::pnorm(log(x), mean = meanlog, sd = sdlog,
                lower.tail = FALSE, log.p = TRUE),
    sum
  )
}

#' @rdname distribution_cdf
exp_mod_normal_lcdf <- function(x, mu, sigma, lambda, sum = TRUE) {
  z <- (x - mu) / sigma
  log_weight <- -lambda * (x - mu) + 0.5 * (lambda * sigma)^2
  log_first <- RTMB::pnorm(z, log.p = TRUE)
  log_second <- log_weight + RTMB::pnorm(z - lambda * sigma, log.p = TRUE)
  res <- RTMB::logspace_sub(log_first, log_second)
  .rtmb_finish_log_cdf(res, sum)
}

#' @rdname distribution_cdf
exp_mod_normal_lccdf <- function(x, mu, sigma, lambda, sum = TRUE) {
  z <- (x - mu) / sigma
  log_weight <- -lambda * (x - mu) + 0.5 * (lambda * sigma)^2
  log_first <- RTMB::pnorm(z, lower.tail = FALSE, log.p = TRUE)
  log_second <- log_weight + RTMB::pnorm(z - lambda * sigma, log.p = TRUE)
  res <- RTMB::logspace_add(log_first, log_second)
  .rtmb_finish_log_cdf(res, sum)
}

#' @rdname distribution_cdf
exponential_lcdf <- function(x, rate, sum = TRUE) {
  x <- 0.5 * (x + abs(x))
  .rtmb_finish_log_cdf(log1m_exp(-rate * x), sum)
}

#' @rdname distribution_cdf
exponential_lccdf <- function(x, rate, sum = TRUE) {
  x <- 0.5 * (x + abs(x))
  .rtmb_finish_log_cdf(-rate * x, sum)
}

#' @rdname distribution_cdf
half_normal_lccdf <- function(x, sd, sum = TRUE) {
  x <- 0.5 * (x + abs(x))
  res <- log(2) + RTMB::pnorm(x / sd, lower.tail = FALSE, log.p = TRUE)
  .rtmb_finish_log_cdf(res, sum)
}

#' @rdname distribution_cdf
half_normal_lcdf <- function(x, sd, sum = TRUE) {
  log_survival <- half_normal_lccdf(x, sd, sum = FALSE)
  .rtmb_finish_log_cdf(log1m_exp(log_survival), sum)
}

#' @rdname distribution_cdf
beta_lcdf <- function(x, a, b, sum = TRUE) {
  .rtmb_finish_log_cdf(log(RTMB::pbeta(x, a, b)), sum)
}

#' @rdname distribution_cdf
beta_lccdf <- function(x, a, b, sum = TRUE) {
  .rtmb_finish_log_cdf(log(RTMB::pbeta(1 - x, b, a)), sum)
}

#' @rdname distribution_cdf
gamma_lcdf <- function(x, shape, rate, sum = TRUE) {
  .rtmb_finish_log_cdf(log(RTMB::pgamma(x, shape, rate = rate)), sum)
}

#' @rdname distribution_cdf
gamma_lccdf <- function(x, shape, rate, sum = TRUE) {
  log_cdf <- log(RTMB::pgamma(x, shape, rate = rate))
  .rtmb_finish_log_cdf(log1m_exp(log_cdf), sum)
}

#' @rdname distribution_cdf
inverse_gamma_lcdf <- function(x, shape, scale, sum = TRUE) {
  log_gamma_cdf <- log(RTMB::pgamma(scale / x, shape, rate = 1))
  .rtmb_finish_log_cdf(log1m_exp(log_gamma_cdf), sum)
}

#' @rdname distribution_cdf
inverse_gamma_lccdf <- function(x, shape, scale, sum = TRUE) {
  .rtmb_finish_log_cdf(
    log(RTMB::pgamma(scale / x, shape, rate = 1)),
    sum
  )
}

#' @rdname distribution_cdf
cauchy_lcdf <- function(x, location, scale, sum = TRUE) {
  z <- (x - location) / scale
  .rtmb_finish_log_cdf(log(0.5 + atan(z) / pi), sum)
}

#' @rdname distribution_cdf
cauchy_lccdf <- function(x, location, scale, sum = TRUE) {
  z <- (x - location) / scale
  .rtmb_finish_log_cdf(log(0.5 - atan(z) / pi), sum)
}

#' @rdname distribution_cdf
student_t_lcdf <- function(x, df, mu = 0, sigma = 1, sum = TRUE) {
  z <- (x - mu) / sigma
  u <- 0.5 * (1 + z / sqrt(df + z^2))
  .rtmb_finish_log_cdf(log(RTMB::pbeta(u, df / 2, df / 2)), sum)
}

#' @rdname distribution_cdf
student_t_lccdf <- function(x, df, mu = 0, sigma = 1, sum = TRUE) {
  z <- (x - mu) / sigma
  u <- 0.5 * (1 + z / sqrt(df + z^2))
  .rtmb_finish_log_cdf(log(RTMB::pbeta(1 - u, df / 2, df / 2)), sum)
}

.rtmb_laplace_prob <- function(x, location, scale, lower.tail) {
  z <- (x - location) / scale
  positive <- 0.5 * (z + abs(z))
  signed_mass <- exp(z - positive) - exp(-positive)
  if (lower.tail) {
    0.5 * (1 + signed_mass)
  } else {
    0.5 * (1 - signed_mass)
  }
}

#' @rdname distribution_cdf
laplace_lcdf <- function(x, location = 0, scale = 1, sum = TRUE) {
  .rtmb_finish_log_cdf(
    log(.rtmb_laplace_prob(x, location, scale, lower.tail = TRUE)),
    sum
  )
}

#' @rdname distribution_cdf
laplace_lccdf <- function(x, location = 0, scale = 1, sum = TRUE) {
  .rtmb_finish_log_cdf(
    log(.rtmb_laplace_prob(x, location, scale, lower.tail = FALSE)),
    sum
  )
}

#' @rdname distribution_cdf
logit_normal_lcdf <- function(x, mu, sd, sum = TRUE) {
  x <- .rtmb_unit_interval(x)
  normal_lcdf(log(x) - log1p(-x), mean = mu, sd = sd, sum = sum)
}

#' @rdname distribution_cdf
logit_normal_lccdf <- function(x, mu, sd, sum = TRUE) {
  x <- .rtmb_unit_interval(x)
  normal_lccdf(log(x) - log1p(-x), mean = mu, sd = sd, sum = sum)
}

#' @rdname distribution_cdf
weibull_lcdf <- function(x, shape, scale, sum = TRUE) {
  x <- 0.5 * (x + abs(x))
  .rtmb_finish_log_cdf(log1m_exp(-(x / scale)^shape), sum)
}

#' @rdname distribution_cdf
weibull_lccdf <- function(x, shape, scale, sum = TRUE) {
  x <- 0.5 * (x + abs(x))
  .rtmb_finish_log_cdf(-(x / scale)^shape, sum)
}

.rtmb_unit_interval <- function(x) {
  positive <- 0.5 * (x + abs(x))
  0.5 * (positive + 1 - abs(positive - 1))
}

#' @rdname distribution_cdf
uniform_lcdf <- function(x, a, b, sum = TRUE) {
  prob <- .rtmb_unit_interval((x - a) / (b - a))
  .rtmb_finish_log_cdf(log(prob), sum)
}

#' @rdname distribution_cdf
uniform_lccdf <- function(x, a, b, sum = TRUE) {
  prob <- .rtmb_unit_interval((b - x) / (b - a))
  .rtmb_finish_log_cdf(log(prob), sum)
}

.rtmb_normal_mixture_logcdf <- function(x, pi_w, mean, sd, lower.tail) {
  log_pi <- log(pi_w)
  N <- length(x)
  K <- length(pi_w)
  log_terms <- rtmb_array(log_pi[1] * 0, c(N, K), seed = log_pi[1])
  for (k in seq_len(K)) {
    log_terms[, k] <- log_pi[k] + RTMB::pnorm(
      x,
      mean = mean[k],
      sd = sd[k],
      lower.tail = lower.tail,
      log.p = TRUE
    )
  }
  log_sum_exp_matrix(log_terms)
}

#' @rdname distribution_cdf
normal_mixture_lcdf <- function(x, pi_w, mean, sd, sum = TRUE) {
  .rtmb_finish_log_cdf(
    .rtmb_normal_mixture_logcdf(x, pi_w, mean, sd, lower.tail = TRUE),
    sum
  )
}

#' @rdname distribution_cdf
normal_mixture_lccdf <- function(x, pi_w, mean, sd, sum = TRUE) {
  .rtmb_finish_log_cdf(
    .rtmb_normal_mixture_logcdf(x, pi_w, mean, sd, lower.tail = FALSE),
    sum
  )
}

.rtmb_mixture_logcdf <- function(x, pi_w, log_prob_list) {
  N <- if (is.matrix(x)) nrow(x) else length(x)
  K <- length(log_prob_list)
  seed <- .rtmb_ad_seed(c(list(pi_w), log_prob_list))
  log_terms <- if (is.null(seed)) {
    matrix(0, N, K)
  } else {
    rtmb_array(0, c(N, K), seed = seed)
  }

  for (k in seq_len(K)) {
    log_weight <- if (is.matrix(pi_w)) log(pi_w[, k]) else log(pi_w[k])
    log_terms[, k] <- log_weight + log_prob_list[[k]]
  }
  log_sum_exp_matrix(log_terms)
}

#' @rdname distribution_cdf
mixture_lcdf <- function(x, pi_w, lcdf_list, sum = TRUE) {
  .rtmb_finish_log_cdf(
    .rtmb_mixture_logcdf(x, pi_w, lcdf_list),
    sum
  )
}

#' @rdname distribution_cdf
mixture_lccdf <- function(x, pi_w, lccdf_list, sum = TRUE) {
  .rtmb_finish_log_cdf(
    .rtmb_mixture_logcdf(x, pi_w, lccdf_list),
    sum
  )
}

#' @rdname distribution_cdf
bernoulli_lcdf <- function(x, prob, sum = TRUE) {
  .rtmb_finish_log_cdf(log(RTMB::pbinom(x, 1, prob)), sum)
}

#' @rdname distribution_cdf
bernoulli_lccdf <- function(x, prob, sum = TRUE) {
  log_cdf <- log(RTMB::pbinom(x, 1, prob))
  .rtmb_finish_log_cdf(log1m_exp(log_cdf), sum)
}

#' @rdname distribution_cdf
binomial_lcdf <- function(x, size, prob, sum = TRUE) {
  .rtmb_finish_log_cdf(log(RTMB::pbinom(x, size, prob)), sum)
}

#' @rdname distribution_cdf
binomial_lccdf <- function(x, size, prob, sum = TRUE) {
  log_cdf <- log(RTMB::pbinom(x, size, prob))
  .rtmb_finish_log_cdf(log1m_exp(log_cdf), sum)
}

#' @rdname distribution_cdf
bernoulli_logit_lcdf <- function(x, eta, sum = TRUE) {
  bernoulli_lcdf(x, inv_logit(eta), sum = sum)
}

#' @rdname distribution_cdf
bernoulli_logit_lccdf <- function(x, eta, sum = TRUE) {
  bernoulli_lccdf(x, inv_logit(eta), sum = sum)
}

#' @rdname distribution_cdf
binomial_logit_lcdf <- function(x, size, eta, sum = TRUE) {
  binomial_lcdf(x, size, inv_logit(eta), sum = sum)
}

#' @rdname distribution_cdf
binomial_logit_lccdf <- function(x, size, eta, sum = TRUE) {
  binomial_lccdf(x, size, inv_logit(eta), sum = sum)
}

#' @rdname distribution_cdf
poisson_lcdf <- function(x, mean, sum = TRUE) {
  .rtmb_finish_log_cdf(log(RTMB::ppois(x, mean)), sum)
}

#' @rdname distribution_cdf
poisson_lccdf <- function(x, mean, sum = TRUE) {
  .rtmb_finish_log_cdf(
    log(RTMB::pgamma(mean, floor(x) + 1, rate = 1)),
    sum
  )
}

#' @rdname distribution_cdf
neg_binomial_lcdf <- function(x, size, prob, sum = TRUE) {
  .rtmb_finish_log_cdf(log(RTMB::pbeta(prob, size, floor(x) + 1)), sum)
}

#' @rdname distribution_cdf
neg_binomial_lccdf <- function(x, size, prob, sum = TRUE) {
  .rtmb_finish_log_cdf(
    log(RTMB::pbeta(1 - prob, floor(x) + 1, size)),
    sum
  )
}

#' @rdname distribution_cdf
neg_binomial_2_lcdf <- function(x, mu, size, sum = TRUE) {
  prob <- size / (size + mu)
  neg_binomial_lcdf(x, size, prob, sum = sum)
}

#' @rdname distribution_cdf
neg_binomial_2_lccdf <- function(x, mu, size, sum = TRUE) {
  prob <- size / (size + mu)
  neg_binomial_lccdf(x, size, prob, sum = sum)
}

.rtmb_ordered_logistic_logcdf <- function(x, eta, cutpoints, lower.tail) {
  N <- length(x)
  K <- length(cutpoints) + 1L
  eta_vec <- if (length(eta) == 1L && N > 1L) rep(eta, N) else eta
  res <- rep((eta_vec[1] + cutpoints[1]) * 0, N)

  for (i in seq_len(N)) {
    y <- as.integer(x[i])
    if (y < 1L || y > K) {
      stop("ordered_logistic CDF category must be between 1 and K.", call. = FALSE)
    }
    if (y == K) {
      res[i] <- if (lower.tail) 0 else -Inf
    } else {
      v <- cutpoints[y] - eta_vec[i]
      res[i] <- if (lower.tail) -log1p_exp(-v) else -log1p_exp(v)
    }
  }
  res
}

#' @rdname distribution_cdf
ordered_logistic_lcdf <- function(x, eta, cutpoints, sum = TRUE) {
  .rtmb_finish_log_cdf(
    .rtmb_ordered_logistic_logcdf(x, eta, cutpoints, lower.tail = TRUE),
    sum
  )
}

#' @rdname distribution_cdf
ordered_logistic_lccdf <- function(x, eta, cutpoints, sum = TRUE) {
  .rtmb_finish_log_cdf(
    .rtmb_ordered_logistic_logcdf(x, eta, cutpoints, lower.tail = FALSE),
    sum
  )
}

.rtmb_sequential_logistic_lccdf <- function(x, eta, cutpoints) {
  N <- length(x)
  K <- length(cutpoints) + 1L
  eta_is_matrix <- is.matrix(eta)
  eta_vec <- if (!eta_is_matrix && length(eta) == 1L && N > 1L) {
    rep(eta, N)
  } else {
    eta
  }
  res <- rep(cutpoints[1] * 0, N)

  for (i in seq_len(N)) {
    y <- as.integer(x[i])
    if (y < 1L || y > K) {
      stop("sequential_logistic CDF category must be between 1 and K.", call. = FALSE)
    }
    if (y == K) {
      res[i] <- -Inf
    } else {
      lp <- cutpoints[1] * 0
      for (k in seq_len(y)) {
        eta_ik <- if (eta_is_matrix) eta[i, k] else eta_vec[i]
        lp <- lp - log1p_exp(-(eta_ik - cutpoints[k]))
      }
      res[i] <- lp
    }
  }
  res
}

#' @rdname distribution_cdf
sequential_logistic_lccdf <- function(x, eta, cutpoints, sum = TRUE) {
  .rtmb_finish_log_cdf(
    .rtmb_sequential_logistic_lccdf(x, eta, cutpoints),
    sum
  )
}

#' @rdname distribution_cdf
sequential_logistic_lcdf <- function(x, eta, cutpoints, sum = TRUE) {
  log_survival <- .rtmb_sequential_logistic_lccdf(x, eta, cutpoints)
  K <- length(cutpoints) + 1L
  res <- rep(cutpoints[1] * 0, length(x))
  for (i in seq_along(x)) {
    if (as.integer(x[i]) == K) {
      res[i] <- cutpoints[1] * 0
    } else {
      res[i] <- log1m_exp(log_survival[i])
    }
  }
  .rtmb_finish_log_cdf(res, sum)
}
