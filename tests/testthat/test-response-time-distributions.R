test_that("exp_mod_normal_lpdf matches the closed-form density", {
  x <- c(0.35, 0.48, 0.62)
  mu <- 0.4
  sigma <- 0.08
  lambda <- 7

  z <- (x - mu) / sigma - lambda * sigma
  expected <- log(lambda) +
    lambda * (mu - x) +
    0.5 * (lambda * sigma)^2 +
    stats::pnorm(z, log.p = TRUE)

  expect_equal(exp_mod_normal_lpdf(x, mu, sigma, lambda, sum = FALSE), expected)
  expect_equal(exp_mod_normal_lpdf(x, mu, sigma, lambda), sum(expected))
})

test_that("exp_mod_normal_lpdf works with sampling syntax", {
  dat <- list(y = c(0.35, 0.48, 0.62))
  par <- list(mu = 0.4, sigma = 0.08, lambda = 7)

  fn <- model_code({
    y ~ exp_mod_normal(mu, sigma, lambda)
  })

  expect_equal(
    fn(dat, par),
    exp_mod_normal_lpdf(dat$y, par$mu, par$sigma, par$lambda),
    tolerance = 1e-12
  )
})

test_that("diffusion_lpdf returns row-wise contributions for matrix input", {
  rt <- matrix(
    c(0.45, 0.62, 0.51,
      0.77, 0.38, 0.56),
    nrow = 2,
    byrow = TRUE
  )
  response <- matrix(
    c(1, 0, 1,
      0, 1, 1),
    nrow = 2,
    byrow = TRUE
  )
  alpha <- c(1.1, 1.3)
  tau <- 0.2
  beta <- 0.55
  delta <- c(0.8, 1.2)

  row_lp <- diffusion_lpdf(rt, response, alpha, tau, beta, delta, K_diff = 20, sum = FALSE)
  expected <- c(
    diffusion_lpdf(rt[1, , drop = FALSE], response[1, , drop = FALSE],
                   alpha[1], tau, beta, delta[1], K_diff = 20),
    diffusion_lpdf(rt[2, , drop = FALSE], response[2, , drop = FALSE],
                   alpha[2], tau, beta, delta[2], K_diff = 20)
  )

  expect_equal(row_lp, expected, tolerance = 1e-12)
  expect_equal(diffusion_lpdf(rt, response, alpha, tau, beta, delta, K_diff = 20), sum(row_lp))
})

test_that("diffusion_lpdf works with obs() sampling syntax", {
  dat <- list(
    rt = matrix(c(0.45, 0.62, 0.51, 0.77), nrow = 2, byrow = TRUE),
    choice = matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  )
  par <- list(
    alpha = c(1.1, 1.3),
    tau = 0.2,
    beta = 0.55,
    delta = c(0.8, 1.2)
  )

  fn <- model_code({
    obs(rt, choice) ~ diffusion(alpha, tau, beta, delta, K_diff = 20)
  })

  expect_equal(
    fn(dat, par),
    diffusion_lpdf(dat$rt, dat$choice, par$alpha, par$tau, par$beta, par$delta, K_diff = 20),
    tolerance = 1e-12
  )
})
