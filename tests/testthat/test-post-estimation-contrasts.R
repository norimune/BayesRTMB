test_that("post-estimation methods preserve sum contrasts", {
  set.seed(1)
  dat <- expand.grid(
    subj = factor(seq_len(8)),
    cond = factor(c("A", "B")),
    load = factor(c("hi", "lo"))
  )
  subj_effect <- rnorm(8, sd = 0.3)
  dat$y <- 0.4 * (dat$cond == "B") +
    0.2 * (dat$load == "hi") +
    0.15 * (dat$cond == "B") * (dat$load == "lo") +
    subj_effect[dat$subj] + rnorm(nrow(dat), sd = 0.2)

  model <- rtmb_glmer(
    y ~ cond * load + (1 | subj),
    data = dat,
    family = "gaussian",
    contrasts = "sum"
  )
  fit <- model$classic()

  cell_means <- with(dat, tapply(y, list(load, cond), mean))
  expected_diff <- cell_means[, "B"] - cell_means[, "A"]
  simple <- simple_effects(fit, "cond:load")
  expect_equal(
    setNames(simple$estimate, as.character(simple$load)),
    expected_diff,
    tolerance = 1e-5
  )

  conditional <- conditional_effects(fit, "cond:load")$data
  observed_cells <- as.numeric(cell_means[cbind(
    match(as.character(conditional$load), rownames(cell_means)),
    match(as.character(conditional$cond), colnames(cell_means))
  )])
  expect_equal(conditional$estimate, observed_cells, tolerance = 1e-5)

  # Saved objects from older versions do not contain X_contrasts metadata.
  fit$model$extra$X_contrasts <- NULL
  legacy_simple <- simple_effects(fit, "cond:load")
  expect_equal(legacy_simple$estimate, simple$estimate, tolerance = 1e-8)
})

test_that("rtmb_lmer defaults to sum contrasts and honors treatment contrasts", {
  dat <- data.frame(
    y = c(1, 2, 2, 3, 3, 4, 4, 5),
    x = factor(rep(c("A", "B"), 4)),
    id = factor(rep(seq_len(4), each = 2))
  )

  default_model <- rtmb_lmer(y ~ x + (1 | id), data = dat)
  treatment_model <- rtmb_lmer(
    y ~ x + (1 | id),
    data = dat,
    contrasts = "treatment"
  )

  expect_identical(default_model$contrasts, "sum")
  expect_identical(treatment_model$contrasts, "treatment")
  expect_identical(default_model$extra$X_contrasts$x, "contr.sum")
  expect_identical(treatment_model$extra$X_contrasts$x, "contr.treatment")
})
