mvnormal_reference <- function(x, mean, Sigma) {
  residual <- x - mean
  -0.5 * (
    length(x) * log(2 * pi) +
      as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) +
      sum(residual * solve(Sigma, residual))
  )
}

test_that("multi_normal_lpdf matches the previous stabilized density", {
  Sigma <- matrix(
    c(1.4, 0.3, -0.1,
      0.3, 0.9, 0.2,
      -0.1, 0.2, 1.2),
    nrow = 3,
    byrow = TRUE
  )
  safe_Sigma <- Sigma + diag(diag(Sigma) * 1e-6 + 1e-8)
  mean <- c(0.2, -0.4, 0.7)
  x <- matrix(
    c(0.1, -0.2, 1.0,
      0.8, -0.5, 0.3,
      -0.4, 0.2, 0.9),
    nrow = 3,
    byrow = TRUE
  )

  expected <- apply(x, 1, mvnormal_reference, mean = mean, Sigma = safe_Sigma)

  expect_equal(
    multi_normal_lpdf(x, mean, Sigma, sum = FALSE),
    expected,
    tolerance = 1e-12
  )
  expect_equal(
    multi_normal_lpdf(x, mean, Sigma),
    sum(expected),
    tolerance = 1e-12
  )
  expect_equal(
    multi_normal_lpdf(x[1, ], mean, Sigma),
    expected[1],
    tolerance = 1e-12
  )
})

test_that("multi_normal_lpdf supports observation-specific means", {
  Sigma <- matrix(c(1.2, 0.25, 0.25, 0.8), 2, 2)
  safe_Sigma <- Sigma + diag(diag(Sigma) * 1e-6 + 1e-8)
  x <- matrix(c(0.1, 0.8, -0.3, 0.4), nrow = 2, byrow = TRUE)
  mean <- matrix(c(0.0, 0.5, -0.1, 0.2), nrow = 2, byrow = TRUE)
  expected <- vapply(seq_len(nrow(x)), function(i) {
    mvnormal_reference(x[i, ], mean[i, ], safe_Sigma)
  }, numeric(1))

  expect_equal(
    multi_normal_lpdf(x, mean, Sigma, sum = FALSE),
    expected,
    tolerance = 1e-12
  )
})

test_that("fa_multi_normal_lpdf rank-one path matches explicit covariance", {
  Lambda <- matrix(c(0.8, -0.3, 1.1), ncol = 1)
  psi <- c(0.7, 1.2, 0.9)
  mu <- c(0.2, -0.1, 0.5)
  Sigma <- tcrossprod(Lambda) + diag(psi)
  x <- matrix(
    c(0.1, 0.4, 0.8,
      -0.2, -0.3, 0.6,
      1.0, 0.2, -0.1),
    nrow = 3,
    byrow = TRUE
  )
  expected <- apply(x, 1, mvnormal_reference, mean = mu, Sigma = Sigma)

  expect_equal(
    fa_multi_normal_lpdf(x, mu, Lambda, psi, sum = FALSE),
    expected,
    tolerance = 1e-12
  )
  expect_equal(
    fa_multi_normal_lpdf(x, mu, Lambda, psi),
    sum(expected),
    tolerance = 1e-12
  )
  expect_equal(
    fa_multi_normal_lpdf(x[1, ], mu, Lambda, psi),
    expected[1],
    tolerance = 1e-12
  )
})

test_that("fa_multi_normal_lpdf retains the multi-factor path", {
  Lambda <- matrix(
    c(0.8, 0.1,
      -0.3, 0.6,
      1.1, -0.2),
    nrow = 3,
    byrow = TRUE
  )
  psi <- c(0.7, 1.2, 0.9)
  mu <- c(0.2, -0.1, 0.5)
  Sigma <- tcrossprod(Lambda) + diag(psi)
  x <- c(0.1, 0.4, 0.8)

  expect_equal(
    fa_multi_normal_lpdf(x, mu, Lambda, psi),
    mvnormal_reference(x, mu, Sigma),
    tolerance = 1e-12
  )
})

test_that("fa_multi_normal_lpdf rank-one gradient matches explicit covariance", {
  x <- matrix(c(0.1, 0.4, 0.8, -0.2, -0.3, 0.6), nrow = 2, byrow = TRUE)
  mu <- c(0.2, -0.1, 0.5)
  psi <- c(0.7, 1.2, 0.9)
  parameters <- list(lambda = c(0.8, -0.3, 1.1))

  rank_one_objective <- function(par) {
    Lambda <- matrix(par$lambda, ncol = 1)
    -fa_multi_normal_lpdf(x, mu, Lambda, psi)
  }
  explicit_objective <- function(par) {
    Lambda <- matrix(par$lambda, ncol = 1)
    Sigma <- tcrossprod(Lambda) + diag(psi)
    -sum(RTMB::dmvnorm(x, mu = mu, Sigma = Sigma, log = TRUE))
  }

  rank_one <- RTMB::MakeADFun(rank_one_objective, parameters, silent = TRUE)
  explicit <- RTMB::MakeADFun(explicit_objective, parameters, silent = TRUE)

  expect_equal(rank_one$fn(), explicit$fn(), tolerance = 1e-10)
  expect_equal(rank_one$gr(), explicit$gr(), tolerance = 1e-9)
})
