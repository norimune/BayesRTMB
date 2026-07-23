library(testthat)
library(BayesRTMB)

test_that("formula wrappers can resolve data from the formula environment", {
  skip_if_not_installed("RTMB")

  set.seed(1)
  X <- seq(-1, 1, length.out = 20)
  ID <- factor(rep(seq_len(4), each = 5))
  Y <- 1 + 0.5 * X + rnorm(length(X), sd = 0.2)
  dat <- data.frame(Y = Y, X = X, ID = ID)

  mdl_lm <- suppressMessages(rtmb_lm(Y ~ X))
  mdl_lm_data <- suppressMessages(rtmb_lm(Y ~ X, data = dat))
  expect_true(inherits(mdl_lm, "RTMB_Model"))
  expect_equal(mdl_lm$raw_data$Y, mdl_lm_data$raw_data$Y)
  expect_equal(mdl_lm$raw_data$X, mdl_lm_data$raw_data$X)

  mdl_glm <- suppressMessages(rtmb_glm(Y ~ X, family = "gaussian"))
  mdl_lmer <- suppressMessages(rtmb_lmer(Y ~ X + (1 | ID), laplace = FALSE))
  mdl_glmer <- suppressMessages(rtmb_glmer(Y ~ X + (1 | ID), laplace = FALSE))

  expect_true(inherits(mdl_glm, "RTMB_Model"))
  expect_true(inherits(mdl_lmer, "RTMB_Model"))
  expect_true(inherits(mdl_glmer, "RTMB_Model"))
})

test_that("formula shorthand requires explicit data when data is omitted", {
  skip_if_not_installed("RTMB")

  Y <- c(1, 2, 3)
  expect_error(
    rtmb_lm(Y ~ .),
    "The '\\.' shorthand requires an explicit 'data' argument"
  )
})

test_that("formula auxiliary values remain in the formula environment", {
  skip_if_not_installed("RTMB")

  set.seed(2)
  X <- seq_len(10)
  Y <- rnorm(length(X))
  degree <- 2
  knots <- c(3, 7)

  mdl_poly <- suppressMessages(rtmb_lm(Y ~ poly(X, degree)))
  mdl_spline <- suppressMessages(rtmb_lm(Y ~ splines::bs(X, knots = knots)))

  expected_poly <- colnames(model.matrix(Y ~ poly(X, degree)))
  expected_spline <- colnames(model.matrix(Y ~ splines::bs(X, knots = knots)))
  expect_identical(mdl_poly$extra$X_colnames, expected_poly)
  expect_identical(mdl_spline$extra$X_colnames, expected_spline)
  expect_false("degree" %in% names(mdl_poly$raw_data))
  expect_false("knots" %in% names(mdl_spline$raw_data))
})

test_that("matrix responses remain matrix columns when data is omitted", {
  skip_if_not_installed("RTMB")

  success <- c(2, 3, 4, 5, 6)
  failure <- c(8, 7, 6, 5, 4)
  Y <- cbind(success, failure)
  X <- seq_along(success)

  mdl <- suppressMessages(rtmb_glm(Y ~ X, family = "binomial"))

  expect_true(is.matrix(mdl$raw_data$Y))
  expect_equal(unclass(mdl$raw_data$Y), unclass(Y))
  expect_equal(mdl$data$Y, success)
  expect_equal(mdl$data$trials, success + failure)
})

test_that("data access operators require explicit data", {
  skip_if_not_installed("RTMB")

  dat <- data.frame(Y = 1:5, X = 6:10)
  dollar_formula <- as.formula("dat$Y ~ dat$X", env = environment())
  bracket_formula <- as.formula("dat[['Y']] ~ dat[['X']]", env = environment())

  expect_error(rtmb_lm(dollar_formula), "must use bare names")
  expect_error(rtmb_lm(bracket_formula), "must use bare names")
})

test_that("CWC variables are resolved before formula data", {
  skip_if_not_installed("RTMB")

  set.seed(3)
  X <- rep(seq_len(3), 4)
  ID <- factor(rep(seq_len(4), each = 3))
  Y <- rnorm(length(X))

  mdl <- suppressMessages(rtmb_lmer(
    Y ~ X,
    laplace = FALSE,
    cwc = list(cluster = "ID", pars = "X")
  ))

  expect_true("ID" %in% names(mdl$raw_data))
  expect_equal(
    as.numeric(tapply(mdl$raw_data$X, mdl$raw_data$ID, mean)),
    rep(0, nlevels(ID)),
    tolerance = 1e-12
  )

  expect_error(
    rtmb_lmer(
      Y ~ X,
      laplace = FALSE,
      cwc = list(cluster = ID[-1], pars = "X")
    ),
    "same length as the model data"
  )
})

test_that("rtmb_lmer preserves unquoted CWC columns from explicit data", {
  skip_if_not_installed("RTMB")

  dat <- data.frame(
    y = rnorm(12),
    x = rep(seq_len(3), 4),
    ID = factor(rep(seq_len(4), each = 3))
  )

  mdl <- suppressMessages(rtmb_lmer(
    y ~ x + (1 | ID),
    data = dat,
    laplace = FALSE,
    cwc = list(ID, "all")
  ))

  expect_equal(
    as.numeric(tapply(mdl$raw_data$x, mdl$raw_data$ID, mean)),
    rep(0, nlevels(dat$ID)),
    tolerance = 1e-12
  )
})

test_that("unquoted CWC clusters follow wide-to-long conversion", {
  skip_if_not_installed("RTMB")

  set.seed(4)
  y1 <- rnorm(6)
  y2 <- rnorm(6)
  X <- seq_len(6)
  ID <- factor(seq_len(6))

  mdl <- suppressMessages(rtmb_lmer(
    cbind(y1, y2) ~ condition + X + (1 | ID),
    laplace = FALSE,
    within = list(condition = 2),
    cwc = list(ID, "X")
  ))

  expect_equal(nrow(mdl$raw_data), 12)
  expect_equal(
    as.numeric(tapply(mdl$raw_data$X, mdl$raw_data$ID, mean)),
    rep(0, nlevels(ID)),
    tolerance = 1e-12
  )
})
