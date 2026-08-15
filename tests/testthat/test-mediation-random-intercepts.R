make_mediation_group_data <- function() {
  set.seed(814)
  id <- factor(rep(seq_len(8), each = 4))
  x <- rnorm(length(id))
  u_m <- rnorm(nlevels(id), sd = 0.4)
  u_y <- rnorm(nlevels(id), sd = 0.5)
  m <- 0.6 * x + u_m[id] + rnorm(length(id), sd = 0.5)
  y <- 0.2 * x + 0.7 * m + u_y[id] + rnorm(length(id), sd = 0.5)
  data.frame(id = id, x = x, m = m, y = y)
}

test_that("mediation permits a random intercept in only one equation", {
  dat <- make_mediation_group_data()
  model <- rtmb_mediation(
    list(m ~ x, y ~ x + m + (1 | id)),
    data = dat,
    prior = prior_normal()
  )

  parameter_text <- paste(deparse(model$code$parameters), collapse = "\n")
  model_text <- paste(deparse(model$code$model), collapse = "\n")

  expect_match(parameter_text, "sd_re <- Dim\\(1, lower = 0\\)")
  expect_match(parameter_text, "r_re <- Dim\\(mediation_num_groups, random = TRUE\\)")
  expect_false(grepl("CF_corr_re", parameter_text, fixed = TRUE))
  expect_match(
    model_text,
    "sd_re\\[1\\] \\* r_re\\[mediation_group_idx\\]"
  )
  expect_identical(model$extra$mediation$random_equations, 2L)
  expect_identical(model$extra$mediation$group, "id")
})

test_that("mediation correlates random intercepts across equations", {
  dat <- make_mediation_group_data()
  model <- rtmb_mediation(
    list(m ~ x + (1 | id), y ~ x + m + (1 | id)),
    data = dat,
    prior = prior_normal()
  )

  parameter_text <- paste(deparse(model$code$parameters), collapse = "\n")
  transform_text <- paste(deparse(model$code$transform), collapse = "\n")
  model_text <- paste(deparse(model$code$model), collapse = "\n")

  expect_match(parameter_text, "sd_re <- Dim\\(2, lower = 0\\)")
  expect_match(parameter_text, "CF_corr_re <- Dim")
  expect_match(
    parameter_text,
    "r_re <- Dim\\(c\\(mediation_num_groups, 2\\), random = TRUE\\)"
  )
  expect_match(transform_text, "corr_re <- CF_corr_re %\\*% t\\(CF_corr_re\\)")
  expect_match(model_text, "multi_normal_CF")
  expect_match(model_text, "CF_corr_re ~ lkj_CF_corr\\(1\\)")
  expect_identical(model$extra$mediation$random_equations, c(1L, 2L))
})

test_that("mediation rejects unsupported random-effect structures", {
  dat <- make_mediation_group_data()
  dat$id2 <- dat$id

  expect_error(
    rtmb_mediation(list(m ~ x + (1 + x | id), y ~ x + m), data = dat),
    "Only random intercepts"
  )
  expect_error(
    rtmb_mediation(
      list(m ~ x + (1 | id), y ~ x + m + (1 | id2)),
      data = dat
    ),
    "must use the same grouping variable"
  )
  expect_error(
    rtmb_mediation(
      list(m ~ x + (1 | id) + (1 | id2), y ~ x + m),
      data = dat
    ),
    "at most one random-effect term"
  )
})

test_that("mediation centers predictor uses without changing responses", {
  dat <- make_mediation_group_data()
  model <- rtmb_mediation(
    list(m ~ x, y ~ x + m),
    data = dat,
    centering = "x",
    cwc = list(id, "m")
  )

  expect_equal(as.numeric(model$data$Y_1), dat$m)
  expect_equal(
    as.numeric(model$data$X_1[, "x"]),
    center_grand_mean(dat$x)
  )
  expect_equal(
    as.numeric(model$data$X_2[, "m"]),
    center_within_cluster(dat$m, dat$id)
  )
  transform_text <- paste(deparse(model$code$transform), collapse = "\n")
  expect_match(transform_text, "IE_x_m_y")
  expect_identical(model$extra$mediation$centering, "x")
  expect_identical(
    model$extra$mediation$cwc,
    list(cluster = "id", pars = "m")
  )

  printed <- capture.output(model$print_code())
  expect_true(any(grepl("center_grand_mean", printed, fixed = TRUE)))
  expect_true(any(grepl("center_within_cluster", printed, fixed = TRUE)))

  rebuilt_code <- model$code
  rebuilt_code$setup_env <- NULL
  rebuilt <- rtmb_model(
    data = dat,
    code = rebuilt_code,
    init = model$init,
    silent = TRUE
  )
  expect_equal(rebuilt$data$X_1, model$data$X_1)
  expect_equal(rebuilt$data$X_2, model$data$X_2)
})

test_that("mediation validates centering specifications", {
  dat <- make_mediation_group_data()

  expect_error(
    rtmb_mediation(
      list(m ~ x, y ~ x + m),
      data = dat,
      gmc = "x",
      centering = "m"
    ),
    "Specify only one"
  )
  expect_error(
    rtmb_mediation(
      list(m ~ x, y ~ x + m),
      data = dat,
      cwc = list(id, "missing_variable")
    ),
    "were not found in data"
  )
  expect_error(
    rtmb_mediation(
      list(m ~ x, y ~ x + m),
      data = dat,
      cwc = list(id, "id")
    ),
    "not predictors"
  )
})
