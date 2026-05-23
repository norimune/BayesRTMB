test_that("Wrappers optimize correctly", {
  
  # rtmb_lm
  fit_lm <- rtmb_lm(mpg ~ wt, data = mtcars)
  res_lm <- fit_lm$optimize()
  expect_true(inherits(res_lm, "map_fit"))
  
  # rtmb_ttest
  fit_t <- rtmb_ttest(mpg ~ am, data = mtcars)
  res_t <- fit_t$optimize()
  expect_true(inherits(res_t, "map_fit"))

  # rtmb_corr
  fit_c <- rtmb_corr(mtcars[, c("mpg", "wt")])
  res_c <- fit_c$optimize()
  expect_true(inherits(res_c, "map_fit"))
  
  # rtmb_fa
  # Added a simple fallback to handle mtcars subset
  fit_f <- rtmb_fa(scale(mtcars[, c("mpg", "disp", "hp")]), nfactors = 1)
  res_f <- fit_f$optimize()
  expect_true(inherits(res_f, "map_fit"))
})

test_that("Wrappers sample correctly (skip on CRAN)", {
  skip_on_cran()
  
  fit_lm <- rtmb_lm(mpg ~ wt, data = mtcars)
  # Extremely short sampling just to test code path
  res_lm <- fit_lm$sample(chains = 1, sampling = 10, warmup = 10)
  expect_true(inherits(res_lm, "mcmc_fit"))

  res_lm_dense <- fit_lm$sample(
    chains = 1, sampling = 5, warmup = 25,
    metric = "dense", metric_init = "hessian",
    metric_adaptation = "cumulative",
    nuts_variant = "multinomial"
  )
  expect_true(inherits(res_lm_dense, "mcmc_fit"))
  expect_identical(res_lm_dense$metric_type, "dense")
  expect_identical(res_lm_dense$metric_init, "hessian")
  expect_identical(res_lm_dense$metric_adaptation, "cumulative")
  expect_true(is.matrix(res_lm_dense$metric[[1]]))
  expect_true(is.list(res_lm_dense$warmup_diagnostics))
  expect_true(is.data.frame(res_lm_dense$warmup_diagnostics[[1]]))
  expect_true(all(c("phase", "window", "metric_updated") %in% names(res_lm_dense$warmup_diagnostics[[1]])))

  res_lm_stan_window <- fit_lm$sample(
    chains = 1, sampling = 3, warmup = 20,
    metric_adaptation = "stan_window"
  )
  expect_identical(res_lm_stan_window$metric_adaptation, "stan_window")
  expect_true(is.data.frame(res_lm_stan_window$warmup_diagnostics[[1]]))
  expect_true("metric_updated" %in% names(res_lm_stan_window$warmup_diagnostics[[1]]))
})
