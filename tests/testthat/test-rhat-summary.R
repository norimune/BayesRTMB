test_that("rhat_summary computes R-hat without full summary", {
  set.seed(1)
  draws <- array(
    rnorm(40),
    dim = c(20, 2, 1),
    dimnames = list(NULL, NULL, "b")
  )

  fit <- structure(
    list(
      draws = function(...) draws,
      summary = function(...) stop("summary should not be called", call. = FALSE)
    ),
    class = "mcmc_fit"
  )

  out <- rhat_summary(fit)
  expect_s3_class(out, "rhat_summary")
  expect_named(out, "b")
  expect_true(is.finite(unclass(out)))
})

test_that("rhat_summary keeps constant parameters non-finite like summary", {
  set.seed(1)
  draws <- array(
    c(rnorm(40), rep(1, 40)),
    dim = c(20, 2, 2),
    dimnames = list(NULL, NULL, c("b", "constant"))
  )

  fit <- structure(
    list(draws = function(...) draws),
    class = "mcmc_fit"
  )

  out_all <- rhat_summary(fit, finite = FALSE)
  expect_named(out_all, c("b", "constant"))
  expect_true(is.finite(unclass(out_all)[["b"]]))
  expect_true(is.na(unclass(out_all)[["constant"]]))

  out_finite <- rhat_summary(fit)
  expect_named(out_finite, "b")
})

test_that("rhat_summary excludes generated quantities by default", {
  set.seed(1)
  primary_draws <- array(
    rnorm(40),
    dim = c(20, 2, 1),
    dimnames = list(NULL, NULL, "b")
  )
  generate_draws <- array(
    rnorm(40),
    dim = c(20, 2, 1),
    dimnames = list(NULL, NULL, "log_lik[1]")
  )

  fit <- structure(
    list(
      draws = function(..., inc_generate = FALSE) {
        if (!isTRUE(inc_generate)) return(primary_draws)
        out <- array(
          NA_real_,
          dim = c(20, 2, 2),
          dimnames = list(NULL, NULL, c("b", "log_lik[1]"))
        )
        out[, , 1] <- primary_draws[, , 1]
        out[, , 2] <- generate_draws[, , 1]
        out
      }
    ),
    class = "mcmc_fit"
  )

  expect_named(rhat_summary(fit), "b")
  expect_named(rhat_summary(fit, inc_generate = TRUE), c("b", "log_lik[1]"))
})
