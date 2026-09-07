make_mcmc_continuation_model <- function() {
  code <- rtmb_code(
    setup = {
      Y <- y
    },
    parameters = {
      mu <- Dim()
    },
    transform = {
      mu_twice <- 2 * mu
      report(mu_twice)
    },
    model = {
      Y ~ normal(mu, 1)
      mu ~ normal(0, 1)
    },
    generate = {
      mu_copy <- mu
      report(mu_copy)
    }
  )

  rtmb_model(
    data = list(y = c(-0.25, 0.1, 0.3)),
    code = code,
    init = list(mu = 0),
    silent = TRUE
  )
}

test_that("continue_sampling appends draws and reuses adapted NUTS state", {
  model <- make_mcmc_continuation_model()
  fit <- suppressWarnings(model$sample(
    sampling = 6,
    warmup = 15,
    chains = 2,
    seed = 101,
    metric = "diag",
    init_jitter = 0,
    progress = "none"
  ))

  old_fit <- fit$fit
  old_transform <- fit$transform_fit
  old_generate <- fit$generate_fit
  old_leapfrog <- fit$n_leapfrog
  old_metric <- fit$metric
  old_eps <- fit$eps
  fit$log_ml <- -10

  returned <- suppressWarnings(fit$continue_sampling(
    sampling = 4,
    seed = 202,
    progress = "none"
  ))

  expect_identical(returned, fit)
  expect_equal(dim(fit$fit), c(10L, 2L, 2L))
  expect_equal(fit$fit[seq_len(6), , , drop = FALSE], old_fit)
  expect_equal(fit$transform_fit[seq_len(6), , , drop = FALSE], old_transform)
  expect_equal(fit$generate_fit[seq_len(6), , , drop = FALSE], old_generate)
  expect_equal(
    unname(fit$n_leapfrog[seq_len(6), , drop = FALSE]),
    unname(old_leapfrog)
  )
  expect_equal(fit$metric, old_metric)
  expect_equal(fit$eps, old_eps)
  expect_length(fit$chain_state, 2L)
  expect_equal(fit$sampler_config$continuations, 1L)
  expect_equal(fit$sampler_config$retained_draws, 10L)
  expect_null(fit$log_ml)
})

test_that("continue_sampling can return an extended copy", {
  model <- make_mcmc_continuation_model()
  fit <- suppressWarnings(model$sample(
    sampling = 5,
    warmup = 15,
    chains = 1,
    seed = 303,
    metric = "diag",
    init_jitter = 0,
    progress = "none"
  ))
  original_fit <- fit$fit
  fit$chain_state <- NULL
  fit$sampler_config <- NULL
  original_state <- fit$chain_state

  extended <- suppressWarnings(fit$continue_sampling(
    sampling = 3,
    seed = 404,
    inplace = FALSE,
    progress = "none"
  ))

  expect_false(identical(extended, fit))
  expect_identical(fit$fit, original_fit)
  expect_identical(fit$chain_state, original_state)
  expect_null(fit$sampler_config)
  expect_equal(dim(extended$fit)[1L], 8L)
  expect_equal(extended$fit[seq_len(5), , , drop = FALSE], original_fit)
  expect_length(extended$chain_state, 1L)
})
