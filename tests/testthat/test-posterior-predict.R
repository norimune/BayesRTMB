test_that("regression wrappers provide automatic posterior predictions", {
  dat <- data.frame(
    y = c(-1.2, -0.7, -0.1, 0.2, 0.8, 1.3, 1.7, 2.1),
    x = seq(-1, 1, length.out = 8)
  )
  model <- rtmb_lm(y ~ x, data = dat)
  fit <- model$optimize(se_method = "none")

  yrep <- fit$posterior_predict(draws = 4, seed = 12)
  expect_s3_class(yrep, "rtmb_posterior_predict")
  expect_equal(dim(yrep), c(4, nrow(dat)))
  expect_equal(attr(yrep, "density"), "lpdf")

  yrep_s3 <- posterior_predict(fit, draws = 2, seed = 12)
  expect_equal(dim(yrep_s3), c(2, nrow(dat)))

  check <- fit$pp_check(type = "auto", draws = 4, seed = 12, plot = FALSE)
  expect_s3_class(check, "rtmb_pp_check")
  expect_equal(check$type, "dens")

  stat_check <- fit$pp_check(stat = mean, draws = 4, seed = 12, plot = FALSE)
  expect_length(stat_check$replicated_stat, 4)
  expect_length(stat_check$observed_stat, 1)

  named_stat_check <- fit$pp_check(stat = "mean", draws = 4, seed = 12, plot = FALSE)
  expect_equal(named_stat_check$stat$label, "mean")
  expect_length(named_stat_check$replicated_stat, 4)

  predict_code <- rtmb_code(
    generate = {
      y_rep <- stats::rnorm(N, mean = Intercept + X %*% b, sd = sigma)
      report(y_rep)
    }
  )
  yrep <- fit$posterior_predict(code = predict_code, draws = 3, seed = 44)

  expect_equal(dim(yrep), c(3, nrow(dat)))
  expect_equal(attr(yrep, "variable"), "y_rep")
})

test_that("hierarchical prediction modes use the requested random-effect level", {
  data <- list(
    Y = rep(0, 4),
    X = matrix(numeric(0), nrow = 4, ncol = 0),
    Z_mat = matrix(1, nrow = 4, ncol = 1),
    group_idx = c(1L, 1L, 2L, 2L)
  )
  state <- list(Intercept_c = 2, r_re = c(-1, 1), sd = 0.5)
  spec <- list(
    response = "Y",
    family = "gaussian",
    K = 0L,
    num_categories = NULL,
    has_intercept = TRUE,
    use_centering = TRUE,
    has_offset = FALSE,
    random_terms = list(list(
      z_name = "Z_mat",
      group_name = "group_idx",
      effect_name = "r_re",
      sd_name = "sd",
      corr_name = NULL,
      num_groups = 2L,
      num_ranef = 1L
    ))
  )

  expect_equal(.rtmb_glmer_eta(data, state, spec, "population"), rep(2, 4))
  expect_equal(.rtmb_glmer_eta(data, state, spec, "conditional"), c(1.5, 1.5, 2.5, 2.5))
  expect_length(.rtmb_glmer_eta(data, state, spec, "simulate"), 4)
})

test_that("auto checks distinguish probability densities and masses", {
  gaussian_model <- rtmb_lm(mpg ~ wt, data = mtcars)
  bernoulli_model <- rtmb_glm(am ~ wt, data = mtcars, family = "bernoulli")

  expect_equal(gaussian_model$extra$posterior_predict$density, "lpdf")
  expect_equal(bernoulli_model$extra$posterior_predict$density, "lpmf")

  custom_model <- list(
    extra = list(),
    data = list(Y = 1:3),
    code = list(
      setup = quote({ custom_lpdf <- function(x, location) x * 0 }),
      model = quote({ Y ~ custom(location) })
    )
  )
  expect_equal(.rtmb_infer_density(custom_model), "lpdf")
})
