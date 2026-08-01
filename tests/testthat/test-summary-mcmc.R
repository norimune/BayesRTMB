test_that("summary_mcmc matches MCMC_Fit summary contents", {
  set.seed(123)
  draws <- array(
    rnorm(80 * 2 * 3),
    dim = c(80, 2, 3),
    dimnames = list(
      iteration = NULL,
      chain = c("chain1", "chain2"),
      variable = c("theta[1]", "lp", "theta[2]")
    )
  )

  fit <- MCMC_Fit$new(
    model = list(view = NULL),
    fit = draws,
    random_fit = NULL,
    eps = c(chain1 = 0.1, chain2 = 0.1),
    accept = c(chain1 = 0.9, chain2 = 0.9),
    treedepth = c(chain1 = 2, chain2 = 2),
    laplace = FALSE,
    posterior_mean = numeric()
  )

  expected <- fit$summary(max_rows = NULL)
  observed <- summary_mcmc(draws, max_rows = NULL)

  expect_s3_class(observed, "summary_BayesRTMB")
  expect_equal(as.data.frame(observed), as.data.frame(expected))
  expect_identical(attr(observed, "digits"), 2L)
})

test_that("summary_mcmc selects variables and chains", {
  draws <- array(
    seq_len(20 * 2 * 3),
    dim = c(20, 2, 3),
    dimnames = list(
      NULL,
      c("chain1", "chain2"),
      c("alpha", "beta[1]", "beta[2]")
    )
  )

  out <- summary_mcmc(
    draws,
    pars = "beta",
    chains = 1,
    max_rows = NULL,
    digits = 4
  )

  expect_identical(out$variable, c("beta[1]", "beta[2]"))
  expect_equal(out$mean, c(mean(draws[, 1, 2]), mean(draws[, 1, 3])))
  expect_identical(attr(out, "digits"), 4L)

  excluded <- summary_mcmc(draws, pars = "-beta", max_rows = NULL)
  expect_identical(excluded$variable, "alpha")
})

test_that("summary_mcmc handles unnamed and constant variables", {
  draws <- array(1, dim = c(20, 2, 1))
  out <- summary_mcmc(draws)

  expect_identical(out$variable, "V1")
  expect_equal(out$mean, 1)
  expect_equal(out$map, 1)
  expect_equal(out$q2.5, 1)
  expect_equal(out$q97.5, 1)
  expect_true(all(is.na(out[c("ess_bulk", "ess_tail", "rhat")])))
})

test_that("summary_mcmc validates its array and selectors", {
  draws <- array(rnorm(40), dim = c(10, 2, 2))

  expect_error(summary_mcmc(matrix(rnorm(20), 10, 2)),
               "three-dimensional array")
  expect_error(summary_mcmc(array(letters[1:8], dim = c(2, 2, 2))),
               "must be numeric")
  expect_error(summary_mcmc(draws, chains = 3),
               "specified chains")
  expect_error(summary_mcmc(draws, pars = "missing"),
               "variable name")
  expect_error(summary_mcmc(draws, max_rows = 0),
               "positive integer")
})
