test_that("fixed bounded parameters do not contribute a Jacobian", {
  code <- rtmb_code(
    parameters = {
      rho <- Dim(lower = -1, upper = 1)
    },
    model = {
      rho ~ lkj_corr(1)
    }
  )

  for (rho0 in c(0, 0.3, -0.6)) {
    model <- rtmb_model(
      data = list(),
      code = code,
      init = list(rho = rho0),
      fixed = list(rho = rho0),
      silent = TRUE
    )
    ad_obj <- model$build_ad_obj(jacobian_target = "all")$ad_obj

    expect_length(ad_obj$par, 0L)
    expect_equal(-ad_obj$fn(numeric(0)), 0, tolerance = 1e-12)
  }
})

test_that("calc_log_jacobian excludes only NA-mapped bounded components", {
  par_list <- list(theta = Dim(2, lower = -1, upper = 1))
  para_unc <- list(theta = c(0, 1))
  map <- list(theta = factor(c(NA, 1)))

  expected <- log(2) - 1 - 2 * log1p(exp(-1))

  expect_equal(
    calc_log_jacobian(para_unc, par_list, map = map),
    expected,
    tolerance = 1e-12
  )
})

test_that("sampling and rebuilt bridge targets agree with a fixed bound", {
  code <- rtmb_code(
    parameters = {
      mu <- Dim()
      rho <- Dim(lower = -1, upper = 1)
    },
    model = {
      mu ~ normal(0, 1)
      rho ~ lkj_corr(1)
    }
  )

  model <- rtmb_model(
    data = list(),
    code = code,
    init = list(mu = 0, rho = 0),
    fixed = list(rho = 0),
    silent = TRUE
  )
  fit <- model$sample(
    sampling = 20,
    warmup = 20,
    chains = 1,
    seed = 42,
    metric = "diag",
    progress = "none"
  )

  draws_uc <- fit$unconstrain_draws()
  lp_rebuilt <- apply(draws_uc, 1, fit$log_prob())
  lp_saved <- as.numeric(fit$fit[, , 1])

  expect_equal(lp_rebuilt, lp_saved, tolerance = 1e-10)
})
