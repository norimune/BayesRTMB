test_that("Wrappers optimize correctly", {
  
  # rtmb_lm
  fit_lm <- rtmb_lm(mpg ~ wt, data = mtcars)
  res_lm <- fit_lm$optimize()
  expect_true(inherits(res_lm, "map_fit"))
  expect_true(all(c("status", "selected") %in% names(res_lm$opt_history)))
  expect_true(any(res_lm$opt_history$selected))
  expect_true(all(res_lm$opt_history$status %in% c("converged", "not converged") | grepl("convergence", res_lm$opt_history$status)))
  
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
  expect_identical(res_lm$nuts_variant, "multinomial")
  expect_identical(res_lm$metric_requested, "auto")
  expect_true(res_lm$metric_type %in% c("diag", "hybrid", "mixed"))
  expect_identical(res_lm$metric_adaptation, "stan_window")

  res_lm_dense <- suppressWarnings(fit_lm$sample(
    chains = 1, sampling = 5, warmup = 25,
    metric = "dense", metric_init = "hessian",
    metric_adaptation = "cumulative",
    nuts_variant = "multinomial", delta = 0.99
  ))
  expect_true(inherits(res_lm_dense, "mcmc_fit"))
  expect_identical(res_lm_dense$metric_type, "dense")
  expect_identical(res_lm_dense$metric_init, "hessian")
  expect_identical(res_lm_dense$metric_adaptation, "cumulative")
  expect_true(is.matrix(res_lm_dense$metric[[1]]))
  expect_true(is.list(res_lm_dense$warmup_diagnostics))
  expect_true(is.data.frame(res_lm_dense$warmup_diagnostics[[1]]))
  expect_true(all(c("phase", "window", "metric_updated") %in% names(res_lm_dense$warmup_diagnostics[[1]])))

  dat_re <- list(Y = rnorm(8), group = rep(1:4, each = 2), G = 4)
  code_re <- rtmb_code(
    parameters = {
      mu = Dim()
      sigma = Dim(lower = 0)
      tau = Dim(lower = 0)
      r = Dim(G, random = TRUE)
    },
    model = {
      Y ~ normal(mu + r[group] * tau, sigma)
      r ~ normal(0, 1)
      tau ~ exponential(1)
      sigma ~ exponential(1)
    }
  )
  res_hybrid <- rtmb_model(dat_re, code_re)$sample(
    chains = 1, sampling = 3, warmup = 20,
    metric = "hybrid", nuts_variant = "multinomial", delta = 0.99
  )
  expect_identical(res_hybrid$metric_type, "hybrid")
  expect_true(is.list(res_hybrid$metric[[1]]))
  expect_true(length(res_hybrid$metric[[1]]$dense_idx) > 0)
  expect_true(length(res_hybrid$metric[[1]]$diag_idx) > 0)

  res_lm_stan_window <- fit_lm$sample(
    chains = 1, sampling = 3, warmup = 20,
    metric = "diag", delta = 0.99,
    metric_adaptation = "stan_window"
  )
  expect_identical(res_lm_stan_window$metric_adaptation, "stan_window")
  expect_true(is.data.frame(res_lm_stan_window$warmup_diagnostics[[1]]))
  expect_true("metric_updated" %in% names(res_lm_stan_window$warmup_diagnostics[[1]]))
})

test_that("Wrappers reject unsupported priors and unused dots", {
  Ybin <- matrix(sample(0:1, 30, replace = TRUE), nrow = 10)
  Yfa <- scale(mtcars[1:20, c("mpg", "disp", "hp")])

  expect_error(rtmb_irt(Ybin, prior = prior_rhs()), "supports only")
  expect_error(rtmb_fa(Yfa, nfactors = 1, prior = prior_jzs()), "supports only")
  expect_error(rtmb_ttest(mpg ~ am, data = mtcars, prior = prior_rhs()), "supports only")
  expect_error(rtmb_table(matrix(c(3, 4, 5, 6), 2), prior = prior_jzs()), "supports only")

  expect_error(rtmb_lmer(mpg ~ wt + (1 | cyl), data = mtcars, typo_arg = 1), "unused argument")
  expect_error(rtmb_ttest(mpg ~ am, data = mtcars, typo_arg = 1), "unused argument")
  expect_error(rtmb_table(matrix(c(1, 2, 3, 4), 2), typo_arg = 1), "unused argument")
})
