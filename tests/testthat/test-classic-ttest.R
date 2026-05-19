library(testthat)
library(BayesRTMB)

test_that("classic ttest matches stats::t.test with equal variances", {
  set.seed(1)
  Y1 <- rnorm(100, mean = 4.1, sd = 1.0)
  Y2 <- rnorm(100, mean = 4.4, sd = 1.0)

  mdl <- rtmb_ttest(Y1, Y2, prior = prior_flat(), var.equal = TRUE)
  fit <- mdl$classic()

  # Use raw results for precise comparison (combine all components)
  all_dfs <- list()
  if (is.data.frame(fit$df_fixed)) all_dfs$fixed <- fit$df_fixed
  if (is.data.frame(fit$random_effects)) all_dfs$random <- fit$random_effects
  if (is.data.frame(fit$df_transform)) all_dfs$transform <- fit$df_transform
  tab <- do.call(rbind, unname(all_dfs))
  
  tt <- t.test(Y1, Y2, var.equal = TRUE)

  est <- unname(tt$estimate[1] - tt$estimate[2])
  se <- abs(est / unname(tt$statistic))

  expect_equal(tab["diff", "Estimate"], est, tolerance = 1e-6)
  expect_equal(tab["diff", "Std. Error"], se, tolerance = 1e-5)
  expect_equal(tab["diff", "Lower 95%"], unname(tt$conf.int[1]), tolerance = 1e-5)
  expect_equal(tab["diff", "Upper 95%"], unname(tt$conf.int[2]), tolerance = 1e-5)
  expect_equal(tab["diff", "df"], unname(tt$parameter), tolerance = 1e-8)
  
  t_val <- tab["diff", "Estimate"] / tab["diff", "Std. Error"]
  p_val <- 2 * pt(-abs(t_val), df = tab["diff", "df"])
  expect_equal(p_val, tt$p.value, tolerance = 1e-5)
  
  # Verify group means CI (derived quantities)
  n1 <- length(Y1); n2 <- length(Y2)
  s1 <- var(Y1); s2 <- var(Y2)
  sp <- sqrt(((n1-1)*s1 + (n2-1)*s2) / (n1+n2-2))
  se0 <- sp * sqrt(1/n1)
  se1 <- sp * sqrt(1/n2)
  df_t <- n1 + n2 - 2
  crit_t <- qt(0.975, df = df_t)
  
  expect_equal(tab["mean0", "Estimate"], mean(Y1), tolerance = 1e-6)
  expect_equal(tab["mean0", "Lower 95%"], mean(Y1) - crit_t * se0, tolerance = 1e-5)
  expect_equal(tab["mean1", "Estimate"], mean(Y2), tolerance = 1e-6)
  expect_equal(tab["mean1", "Lower 95%"], mean(Y2) - crit_t * se1, tolerance = 1e-5)

  coefs <- fit$summary(max_rows = 100)$coefficients
  expect_false("mean" %in% rownames(coefs))
  expect_true("total_mean" %in% rownames(coefs))
  expect_false(is.na(coefs["diff", "t value"]))
  expect_true(is.na(coefs["total_mean", "t value"]))
  expect_equal(coefs["total_mean", "Pr"], "")
  expect_equal(coefs["sd", "Pr"], "")
})

test_that("classic paired ttest matches stats::t.test paired", {
  set.seed(1)
  Y1 <- rnorm(30, mean = 1)
  Y2 <- Y1 + rnorm(30, mean = 0.2, sd = 1)

  mdl <- rtmb_ttest(Y1, Y2, prior = prior_flat(), paired = TRUE)
  fit <- mdl$classic()

  # Use raw results for precise comparison
  all_dfs <- list()
  if (is.data.frame(fit$df_fixed)) all_dfs$fixed <- fit$df_fixed
  if (is.data.frame(fit$random_effects)) all_dfs$random <- fit$random_effects
  if (is.data.frame(fit$df_transform)) all_dfs$transform <- fit$df_transform
  tab <- do.call(rbind, unname(all_dfs))
  
  tt <- t.test(Y1, Y2, paired = TRUE)

  est <- mean(Y1 - Y2)
  se <- abs(est / unname(tt$statistic))

  expect_equal(tab["diff", "Estimate"], est, tolerance = 1e-6)
  expect_equal(tab["diff", "Std. Error"], se, tolerance = 1e-5)
  expect_equal(tab["diff", "df"], unname(tt$parameter), tolerance = 1e-8)
})

test_that("classic auto df selects Satterthwaite for heteroscedastic t-test", {
  set.seed(1)
  Y1 <- rnorm(40, mean = 4.1, sd = 1.0)
  Y2 <- rnorm(70, mean = 4.4, sd = 2.0)

  mdl <- rtmb_ttest(Y1, Y2, prior = prior_flat(), var.equal = FALSE)
  fit <- mdl$classic()

  expect_equal(fit$df_method, "satterthwaite")
})

test_that("classic heteroscedastic t-test matches optimize with marginal Satterthwaite", {
  set.seed(1)
  Y1 <- rnorm(40, mean = 4.1, sd = 1.0)
  Y2 <- rnorm(70, mean = 4.4, sd = 2.0)

  mdl <- rtmb_ttest(Y1, Y2, prior = prior_flat(), var.equal = FALSE)

  fit_c <- mdl$classic()
  fit_o <- mdl$optimize(
    marginal = mdl$extra$marginal,
    df_method = "satterthwaite"
  )

  # Combine results for comparison
  df_c <- rbind(fit_c$df_fixed, fit_c$df_transform)
  df_o <- rbind(fit_o$df_fixed, fit_o$df_transform)

  expect_equal(df_c["diff", "Estimate"], df_o["diff", "Estimate"], tolerance = 1e-6)
  expect_equal(df_c["diff", "Std. Error"], df_o["diff", "Std. Error"], tolerance = 1e-5)
  expect_equal(df_c["diff", "Lower 95%"], df_o["diff", "Lower 95%"], tolerance = 1e-5)
  expect_equal(df_c["diff", "Upper 95%"], df_o["diff", "Upper 95%"], tolerance = 1e-5)
  expect_equal(df_c["diff", "df"], df_o["diff", "df"], tolerance = 1e-4)
})

test_that("classic heteroscedastic t-test is close to stats Welch t-test", {
  set.seed(1)
  Y1 <- rnorm(40, mean = 4.1, sd = 1.0)
  Y2 <- rnorm(70, mean = 4.4, sd = 2.0)

  mdl <- rtmb_ttest(Y1, Y2, prior = prior_flat(), var.equal = FALSE)
  fit <- mdl$classic()

  # Combine results
  df_c <- rbind(fit$df_fixed, fit$df_transform)
  tt <- t.test(Y1, Y2, var.equal = FALSE)

  expect_equal(df_c["diff", "Estimate"], unname(tt$estimate[1] - tt$estimate[2]), tolerance = 1e-6)
  # Welch and Satterthwaite are close but not identical
  expect_equal(df_c["diff", "Std. Error"], abs(unname(tt$estimate[1] - tt$estimate[2]) / unname(tt$statistic)), tolerance = 1e-3)
  expect_equal(df_c["diff", "df"], unname(tt$parameter), tolerance = 1.0)
})
