test_that("optimize() marginal MAP and Satterthwaite DF work correctly", {
  set.seed(1)
  Y1 <- rnorm(40, mean = 10, sd = 2)
  Y2 <- rnorm(40, mean = 11, sd = 2)

  mdl_tt <- rtmb_ttest(Y1, Y2, prior = prior_jzs(r = 0.707))

  fit_tt <- mdl_tt$optimize(
    marginal = "auto",
    df_method = "satterthwaite"
  )

  expect_s3_class(fit_tt, "map_fit")
  expect_true(all(c("total_mean", "delta") %in% fit_tt$marginal_vars))
  expect_true(all(c("total_mean", "delta") %in% fit_tt$laplace_random_vars))
  expect_true(fit_tt$show_df)

  tab_tt <- fit_tt$summary()
  expect_true("df" %in% colnames(tab_tt))

  # Verify DF for derived quantities
  for (nm in c("diff", "mean0", "mean1")) {
    df_val <- tab_tt[nm, "df"]
    expect_false(is.na(df_val))
    expect_true(df_val < 1e12)
  }
  
  # delta should be Inf or large finite due to low sensitivity in this balanced case
  expect_true(is.infinite(tab_tt["delta", "df"]))
})

test_that("JZS unequal-variance t-test uses an explicit delta parameter", {
  set.seed(2)
  Y1 <- rnorm(20, mean = 10, sd = 1)
  Y2 <- rnorm(25, mean = 11, sd = 2)

  mdl_tt <- rtmb_ttest(Y1, Y2, prior = prior_jzs(r = 0.707), var.equal = FALSE)

  expect_true(all(c("total_mean", "sd", "delta") %in% names(mdl_tt$par_list)))
  expect_false(any(c("mean0", "mean1") %in% names(mdl_tt$par_list)))
  expect_true(all(c("total_mean", "delta") %in% mdl_tt$extra$marginal))

  fit_null <- mdl_tt$optimize(fixed = list(delta = 0), se_method = "none")
  expect_s3_class(fit_null, "map_fit")
})

test_that("Mixed model metadata separation and profile() logic", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  mdl_mm <- rtmb_lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

  fit_mm <- mdl_mm$optimize(laplace = TRUE, marginal = "auto")

  expect_true(all(c("Intercept", "b") %in% fit_mm$marginal_vars))
  expect_true(all(c("Intercept", "b", "r_re") %in% fit_mm$laplace_random_vars))
  expect_gt(length(fit_mm$laplace_random_vars), length(fit_mm$marginal_vars))

  # Marginalized vars should be blocked from profiling
  expect_error(fit_mm$profile(pars = "Intercept"), "marginalized")

  # Non-marginalized fixed parameter should be profile-able in a non-marginalized fit
  fit_mm_none <- mdl_mm$optimize(marginal = "none")
  prof_mm <- fit_mm_none$profile(pars = "b[Days]", quiet = TRUE)
  expect_false(is.null(prof_mm))
})

test_that("se_method = 'none' and df_method = 'bw' rejection", {
  Y1 <- rnorm(20); Y2 <- rnorm(20)
  mdl_tt <- rtmb_ttest(Y1, Y2)

  expect_warning(
    fit_none <- mdl_tt$optimize(se_method = "none", df_method = "satterthwaite"),
    "df_method is ignored"
  )
  tab_none <- fit_none$summary()
  
  expect_false(fit_none$show_df)
  expect_true(all(is.na(tab_none$`Std. Error`)))
  expect_false("df" %in% colnames(tab_none))

  expect_error(mdl_tt$optimize(df_method = "bw"), "not supported")
})

test_that("classic() regression and prior restriction", {
  mdl_lm <- rtmb_lm(mpg ~ wt, data = mtcars)
  fit_c <- mdl_lm$classic()
  expect_s3_class(fit_c, "Classic_Fit")

  mdl_p <- rtmb_lm(mpg ~ wt, data = mtcars, prior = prior_normal())
  expect_error(mdl_p$classic(), "prior_flat")
})
