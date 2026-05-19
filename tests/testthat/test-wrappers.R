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
})
