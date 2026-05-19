test_that("classic FA summary omits test columns", {
  fit <- rtmb_fa(scale(mtcars[, c("mpg", "disp", "hp")]), nfactors = 1)$classic()
  coefs <- fit$summary(max_rows = 100)$coefficients

  expect_false("df" %in% names(coefs))
  expect_false("t value" %in% names(coefs))
  expect_false("z value" %in% names(coefs))
  expect_false("Pr" %in% names(coefs))
})

test_that("classic wrapper scale rows omit test output", {
  fit_lm <- rtmb_lm(mpg ~ wt, data = mtcars)$classic()
  coefs_lm <- fit_lm$summary(max_rows = 100)$coefficients

  expect_true("sigma" %in% rownames(coefs_lm))
  expect_true(is.na(coefs_lm["sigma", "t value"]))
  expect_equal(coefs_lm["sigma", "Pr"], "")

  scale_like <- data.frame(
    Estimate = c(1, 1),
    `Std. Error` = c(0.1, 0.1),
    `Lower 95%` = c(0.8, 0.8),
    `Upper 95%` = c(1.2, 1.2),
    df = c(10, 10),
    check.names = FALSE,
    row.names = c("sigma1", "sd_group")
  )
  fit_mock <- Classic_Fit$new(model = list(type = "mediation"), df_fixed = scale_like, sd_rep = list())
  coefs_mock <- fit_mock$summary(max_rows = 100)$coefficients
  expect_true(all(is.na(coefs_mock[, "t value"])))
  expect_true(all(coefs_mock[, "Pr"] == ""))

  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")

  fit <- rtmb_lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)$classic()
  coefs <- fit$summary(max_rows = 100)$coefficients
  scale_rows <- grepl("^(sigma|sd)", rownames(coefs))

  expect_true(any(scale_rows))
  expect_true(all(is.na(coefs[scale_rows, "t value"])))
  expect_true(all(coefs[scale_rows, "Pr"] == ""))
})
