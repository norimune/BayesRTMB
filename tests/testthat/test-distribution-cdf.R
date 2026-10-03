test_that("continuous log-CDF functions match reference distributions", {
  x <- c(-1, 0.2, 1.5)
  expect_equal(normal_lcdf(x, 0.3, 1.2, sum = FALSE),
               pnorm(x, 0.3, 1.2, log.p = TRUE))
  expect_equal(normal_lccdf(x, 0.3, 1.2, sum = FALSE),
               pnorm(x, 0.3, 1.2, lower.tail = FALSE, log.p = TRUE))

  x_pos <- c(0.2, 0.8, 2)
  expect_equal(lognormal_lcdf(x_pos, 0.1, 0.7, sum = FALSE),
               plnorm(x_pos, 0.1, 0.7, log.p = TRUE))
  expect_equal(lognormal_lccdf(x_pos, 0.1, 0.7, sum = FALSE),
               plnorm(x_pos, 0.1, 0.7, lower.tail = FALSE, log.p = TRUE))
  expect_equal(exponential_lcdf(x_pos, 1.3, sum = FALSE),
               pexp(x_pos, 1.3, log.p = TRUE))
  expect_equal(exponential_lccdf(x_pos, 1.3, sum = FALSE),
               pexp(x_pos, 1.3, lower.tail = FALSE, log.p = TRUE))

  expect_equal(beta_lcdf(c(0.1, 0.4, 0.8), 2.2, 3.1, sum = FALSE),
               pbeta(c(0.1, 0.4, 0.8), 2.2, 3.1, log.p = TRUE))
  expect_equal(beta_lccdf(c(0.1, 0.4, 0.8), 2.2, 3.1, sum = FALSE),
               pbeta(c(0.1, 0.4, 0.8), 2.2, 3.1,
                     lower.tail = FALSE, log.p = TRUE))
  expect_equal(gamma_lcdf(x_pos, 2.2, 1.3, sum = FALSE),
               pgamma(x_pos, 2.2, rate = 1.3, log.p = TRUE))
  expect_equal(gamma_lccdf(x_pos, 2.2, 1.3, sum = FALSE),
               pgamma(x_pos, 2.2, rate = 1.3,
                      lower.tail = FALSE, log.p = TRUE))

  inv_x <- c(0.5, 1.2, 3)
  expect_equal(inverse_gamma_lcdf(inv_x, 2.5, 1.4, sum = FALSE),
               pgamma(1.4 / inv_x, 2.5, lower.tail = FALSE, log.p = TRUE))
  expect_equal(inverse_gamma_lccdf(inv_x, 2.5, 1.4, sum = FALSE),
               pgamma(1.4 / inv_x, 2.5, log.p = TRUE))

  expect_equal(cauchy_lcdf(x, 0.3, 1.2, sum = FALSE),
               pcauchy(x, 0.3, 1.2, log.p = TRUE))
  expect_equal(cauchy_lccdf(x, 0.3, 1.2, sum = FALSE),
               pcauchy(x, 0.3, 1.2, lower.tail = FALSE, log.p = TRUE))
  expect_equal(student_t_lcdf(x, 5, 0.3, 1.2, sum = FALSE),
               pt((x - 0.3) / 1.2, 5, log.p = TRUE), tolerance = 1e-12)
  expect_equal(student_t_lccdf(x, 5, 0.3, 1.2, sum = FALSE),
               pt((x - 0.3) / 1.2, 5,
                  lower.tail = FALSE, log.p = TRUE), tolerance = 1e-12)

  expect_equal(weibull_lcdf(x_pos, 1.7, 1.2, sum = FALSE),
               pweibull(x_pos, 1.7, 1.2, log.p = TRUE))
  expect_equal(weibull_lccdf(x_pos, 1.7, 1.2, sum = FALSE),
               pweibull(x_pos, 1.7, 1.2,
                        lower.tail = FALSE, log.p = TRUE))
  expect_equal(uniform_lcdf(c(-2, 0, 2), -1, 1, sum = FALSE),
               punif(c(-2, 0, 2), -1, 1, log.p = TRUE))
  expect_equal(uniform_lccdf(c(-2, 0, 2), -1, 1, sum = FALSE),
               punif(c(-2, 0, 2), -1, 1,
                     lower.tail = FALSE, log.p = TRUE))
})

test_that("special continuous log-CDF functions match direct calculations", {
  x <- c(-0.5, 0.5, 1.5)
  mu <- 0.2
  sigma <- 0.8
  lambda <- 1.1
  z <- (x - mu) / sigma
  exgauss_cdf <- pnorm(z) -
    exp(-lambda * (x - mu) + 0.5 * (lambda * sigma)^2) *
    pnorm(z - lambda * sigma)
  expect_equal(exp_mod_normal_lcdf(x, mu, sigma, lambda, sum = FALSE),
               log(exgauss_cdf), tolerance = 1e-12)
  expect_equal(exp_mod_normal_lccdf(x, mu, sigma, lambda, sum = FALSE),
               log1p(-exgauss_cdf), tolerance = 1e-12)

  half_x <- c(0, 0.5, 2)
  half_cdf <- 2 * pnorm(half_x / 1.3) - 1
  expect_equal(half_normal_lcdf(half_x, 1.3, sum = FALSE), log(half_cdf))
  expect_equal(half_normal_lccdf(half_x, 1.3, sum = FALSE),
               log1p(-half_cdf))

  laplace_cdf <- function(q, location, scale) {
    ifelse(q < location,
           0.5 * exp((q - location) / scale),
           1 - 0.5 * exp(-(q - location) / scale))
  }
  laplace_prob <- laplace_cdf(x, 0.2, 1.3)
  expect_equal(laplace_lcdf(x, 0.2, 1.3, sum = FALSE), log(laplace_prob))
  expect_equal(laplace_lccdf(x, 0.2, 1.3, sum = FALSE),
               log1p(-laplace_prob))

  logit_x <- c(0.1, 0.4, 0.8)
  expect_equal(logit_normal_lcdf(logit_x, 0.2, 1.1, sum = FALSE),
               pnorm(qlogis(logit_x), 0.2, 1.1, log.p = TRUE))
  expect_equal(logit_normal_lccdf(logit_x, 0.2, 1.1, sum = FALSE),
               pnorm(qlogis(logit_x), 0.2, 1.1,
                     lower.tail = FALSE, log.p = TRUE))

  pi_w <- c(0.3, 0.7)
  mix_cdf <- pi_w[1] * pnorm(x, -0.5, 1) +
    pi_w[2] * pnorm(x, 1, 0.6)
  expect_equal(normal_mixture_lcdf(
    x, pi_w, c(-0.5, 1), c(1, 0.6), sum = FALSE
  ), log(mix_cdf))
  expect_equal(normal_mixture_lccdf(
    x, pi_w, c(-0.5, 1), c(1, 0.6), sum = FALSE
  ), log1p(-mix_cdf))

  component_lcdf <- list(
    normal_lcdf(x, -0.5, 1, sum = FALSE),
    normal_lcdf(x, 1, 0.6, sum = FALSE)
  )
  component_lccdf <- list(
    normal_lccdf(x, -0.5, 1, sum = FALSE),
    normal_lccdf(x, 1, 0.6, sum = FALSE)
  )
  expect_equal(mixture_lcdf(x, pi_w, component_lcdf, sum = FALSE),
               log(mix_cdf))
  expect_equal(mixture_lccdf(x, pi_w, component_lccdf, sum = FALSE),
               log1p(-mix_cdf))
})

test_that("discrete log-CDF functions use P(X <= x) and P(X > x)", {
  expect_equal(bernoulli_lcdf(0:1, 0.3, sum = FALSE),
               pbinom(0:1, 1, 0.3, log.p = TRUE))
  expect_equal(bernoulli_lccdf(0:1, 0.3, sum = FALSE),
               pbinom(0:1, 1, 0.3, lower.tail = FALSE, log.p = TRUE))

  x <- c(0, 2, 4)
  expect_equal(binomial_lcdf(x, 5, 0.4, sum = FALSE),
               pbinom(x, 5, 0.4, log.p = TRUE))
  expect_equal(binomial_lccdf(x, 5, 0.4, sum = FALSE),
               pbinom(x, 5, 0.4, lower.tail = FALSE, log.p = TRUE))
  expect_equal(binomial_logit_lcdf(x, 5, qlogis(0.4), sum = FALSE),
               pbinom(x, 5, 0.4, log.p = TRUE))

  expect_equal(poisson_lcdf(x, 2.3, sum = FALSE),
               ppois(x, 2.3, log.p = TRUE))
  expect_equal(poisson_lccdf(x, 2.3, sum = FALSE),
               ppois(x, 2.3, lower.tail = FALSE, log.p = TRUE))
  expect_equal(neg_binomial_lcdf(x, 2.2, 0.4, sum = FALSE),
               pnbinom(x, 2.2, 0.4, log.p = TRUE))
  expect_equal(neg_binomial_lccdf(x, 2.2, 0.4, sum = FALSE),
               pnbinom(x, 2.2, 0.4, lower.tail = FALSE, log.p = TRUE))
  expect_equal(neg_binomial_2_lcdf(x, 2.5, 2.2, sum = FALSE),
               pnbinom(x, 2.2, mu = 2.5, log.p = TRUE))
  expect_equal(neg_binomial_2_lccdf(x, 2.5, 2.2, sum = FALSE),
               pnbinom(x, 2.2, mu = 2.5,
                       lower.tail = FALSE, log.p = TRUE))
})

test_that("ordered and sequential log-CDF functions match category sums", {
  cutpoints <- c(-1, 0.8)
  eta <- c(-0.2, 0.3, 1)
  ordered_cdf <- c(
    plogis(cutpoints[1] - eta[1]),
    plogis(cutpoints[2] - eta[2]),
    1
  )
  expect_equal(ordered_logistic_lcdf(1:3, eta, cutpoints, sum = FALSE),
               log(ordered_cdf))
  expect_equal(ordered_logistic_lccdf(1:3, eta, cutpoints, sum = FALSE),
               log1p(-ordered_cdf))

  sequential_prob <- exp(vapply(
    1:3,
    function(y) sequential_logistic_lpmf(
      y, eta = 0.3, cutpoints = cutpoints, sum = FALSE
    ),
    numeric(1)
  ))
  sequential_cdf <- cumsum(sequential_prob)
  expect_equal(sequential_logistic_lcdf(
    1:3, eta = 0.3, cutpoints = cutpoints, sum = FALSE
  ), log(sequential_cdf), tolerance = 1e-12)
  expect_equal(sequential_logistic_lccdf(
    1:3, eta = 0.3, cutpoints = cutpoints, sum = FALSE
  ), log1p(-sequential_cdf), tolerance = 1e-12)
})

test_that("log-CDF functions are compatible with RTMB tapes", {
  functions <- list(
    normal = function(p) normal_lcdf(0.2, p[1], exp(p[2])),
    lognormal = function(p) lognormal_lccdf(1.2, p[1], exp(p[2])),
    exgauss = function(p) exp_mod_normal_lccdf(
      0.5, p[1], exp(p[2]), exp(p[3])
    ),
    exponential = function(p) exponential_lcdf(0.8, exp(p[3])),
    half_normal = function(p) half_normal_lcdf(0.8, exp(p[2])),
    beta = function(p) beta_lccdf(0.4, exp(p[3]), 2),
    gamma = function(p) gamma_lccdf(1.1, exp(p[3]), 2),
    inverse_gamma = function(p) inverse_gamma_lccdf(1.3, exp(p[3]), 2),
    cauchy = function(p) cauchy_lcdf(0.1, p[1], exp(p[2])),
    student = function(p) student_t_lccdf(
      0.3, exp(p[3]) + 1, p[1], exp(p[2])
    ),
    laplace = function(p) laplace_lcdf(0.2, p[1], exp(p[2])),
    logit_normal = function(p) logit_normal_lccdf(
      0.7, p[1], exp(p[2])
    ),
    weibull = function(p) weibull_lcdf(1.1, exp(p[3]), exp(p[2])),
    uniform = function(p) uniform_lccdf(0.4, p[1] - 2, p[1] + 2),
    binomial = function(p) binomial_lcdf(2, 5, inv_logit(p[4])),
    poisson = function(p) poisson_lccdf(2, exp(p[3])),
    neg_binomial = function(p) neg_binomial_2_lcdf(2, exp(p[3]), 2),
    mixture = function(p) {
      pi_w <- rtmb_vector(0, 2, seed = p[1])
      means <- rtmb_vector(0, 2, seed = p[1])
      sds <- rtmb_vector(0, 2, seed = p[1])
      pi_w[1] <- inv_logit(p[4])
      pi_w[2] <- 1 - pi_w[1]
      means[1] <- p[1]
      means[2] <- p[1] + 1
      sds[1] <- exp(p[2])
      sds[2] <- exp(p[2]) + 0.2
      normal_mixture_lcdf(0.1, pi_w, means, sds)
    },
    generic_mixture = function(p) {
      pi_w <- rtmb_vector(0, 2, seed = p[1])
      pi_w[1] <- inv_logit(p[4])
      pi_w[2] <- 1 - pi_w[1]
      component_lcdf <- list(
        normal_lcdf(0.1, p[1], exp(p[2]), sum = FALSE),
        exponential_lcdf(0.1, exp(p[3]), sum = FALSE)
      )
      mixture_lcdf(0.1, pi_w, component_lcdf)
    },
    ordered = function(p) {
      cuts <- rtmb_vector(0, 2, seed = p[1])
      cuts[1] <- p[3] - 1
      cuts[2] <- p[3] + 1
      ordered_logistic_lcdf(2, p[1], cuts)
    },
    sequential = function(p) {
      cuts <- rtmb_vector(0, 2, seed = p[1])
      cuts[1] <- p[3] - 1
      cuts[2] <- p[3] + 1
      sequential_logistic_lccdf(1, p[1], cuts)
    },
    sequential_last = function(p) {
      cuts <- rtmb_vector(0, 2, seed = p[1])
      cuts[1] <- p[3] - 1
      cuts[2] <- p[3] + 1
      sequential_logistic_lcdf(3, p[1], cuts)
    }
  )
  at <- c(0.1, 0.2, 0.3, 0.4)

  for (f in functions) {
    tape <- RTMB::MakeTape(f, rep(0, 4))
    value <- tape(at)
    jacobian <- attr(tape, "methods")$jacobian(at)
    expect_true(all(is.finite(value)))
    expect_true(all(is.finite(jacobian)))
  }
})

test_that("log-CDF functions resolve inside rtmb_code", {
  code <- rtmb_code(
    parameters = {
      mu <- Dim()
      sigma <- Dim(lower = 0)
    },
    model = {
      lp <- lp + normal_lcdf(upper, mu, sigma)
      lp <- lp + normal_lccdf(lower, mu, sigma)
    }
  )
  model <- rtmb_model(list(lower = -1, upper = 1), code)
  obj <- model$get_ad_obj()

  expect_true(all(is.finite(obj$fn(obj$par))))
  expect_true(all(is.finite(obj$gr(obj$par))))
})
