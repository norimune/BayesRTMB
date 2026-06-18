test_that("MCMC diagnose adds actionable recommendations", {
  fit_array <- array(
    rnorm(20 * 2 * 2),
    dim = c(20, 2, 2),
    dimnames = list(NULL, NULL, c("lp", "theta"))
  )
  divergent <- matrix(FALSE, nrow = 20, ncol = 2)
  divergent[1, 1] <- TRUE

  fake_fit <- list(
    fit = fit_array,
    accept = c(0.55, 0.58),
    treedepth = matrix(10, nrow = 20, ncol = 2),
    max_treedepth = 10,
    n_leapfrog = matrix(10, nrow = 20, ncol = 2),
    divergent = divergent,
    metric = list(
      list(
        type = "hybrid",
        dense = diag(c(1e9, 1)),
        dense_idx = 1:2,
        diag = numeric(0),
        diag_idx = integer(0)
      )
    ),
    metric_type = "hybrid",
    metric_requested = "auto",
    metric_effective = c("diag", "hybrid"),
    metric_auto = list(
      list(effective = "diag", reason = "posterior correlations are low"),
      list(effective = "hybrid", reason = "hybrid retained")
    ),
    metric_init = "identity",
    metric_adaptation = "stan_window",
    warmup_diagnostics = list(
      data.frame(
        accept = 0.8,
        eps = 0.1,
        metric_condition = 10,
        metric_updated = TRUE
      )
    ),
    pd_error_count = c(1L, 0L),
    laplace = FALSE,
    nuts_variant = "multinomial",
    summary = function(...) {
      data.frame(
        variable = "theta",
        rhat = 1.02,
        ess_bulk = 200,
        ess_tail = 300
      )
    }
  )

  diag <- diagnose_mcmc_fit(fake_fit)
  expect_s3_class(diag, "diagnose_BayesRTMB")
  expect_true(all(c("check", "priority", "message") %in% names(diag$recommendations)))
  expect_true(all(c("divergent", "ess", "rhat") %in% diag$recommendations$check))

  printed <- capture.output(print(diag))
  expect_true(any(grepl("^Recommendations:", printed)))

  diag_no_recommend <- diagnose_mcmc_fit(fake_fit, recommend = FALSE)
  expect_equal(nrow(diag_no_recommend$recommendations), 0L)
})
