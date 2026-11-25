test_that("bridge draws and density use the map passed to sample", {
  code <- rtmb_code(
    parameters = { b <- Dim(4) },
    model = { b ~ normal(0, 1) }
  )
  model <- rtmb_model(list(), code, init = list(b = c(0, 0, 0, 0.7)), silent = TRUE)
  fit <- model$sample(
    map = list(b = factor(c(5, 2, 5, NA))), chains = 1,
    warmup = 40, sampling = 100, seed = 42, metric = "diag", progress = "none"
  )

  draws <- fit$unconstrain_draws()
  expect_equal(ncol(draws), 2L)
  expect_equal(unname(draws[100, ]), fit$chain_state[[1]], tolerance = 1e-12)
  expect_equal(apply(draws, 1, fit$log_prob()), as.numeric(fit$fit[, , 1]),
               tolerance = 1e-10)

  # Two normals share the first coordinate; the fourth has been fixed.
  expected <- -log(2 * sqrt(pi)) + dnorm(0.7, log = TRUE)
  expect_equal(as.numeric(fit$bridgesampling(seed = 42)), expected, tolerance = 0.05)
})

test_that("bridge density reconstructs the fixed coordinate of a nonlinear map", {
  code <- rtmb_code(
    parameters = { L <- Dim(c(3, 3), type = "CF_corr") },
    model = { L ~ lkj_CF_corr(2) }
  )
  init <- to_constrained(list(L = c(0.1, 0.2, 0.55)),
                         list(L = Dim(c(3, 3), type = "CF_corr")))
  model <- rtmb_model(list(), code, init = init, silent = TRUE)
  fit <- model$sample(
    map = list(L = factor(c(1, 2, NA))), chains = 1, warmup = 50,
    sampling = 50, seed = 42, metric = "diag", progress = "none"
  )
  expect_equal(apply(fit$unconstrain_draws(), 1, fit$log_prob()),
               as.numeric(fit$fit[, , 1]), tolerance = 1e-10)
})

test_that("bridge sampling rejects failure without updating the cache", {
  code <- rtmb_code(
    parameters = { x <- Dim(lower = 0) },
    model = { x ~ exponential(1) }
  )
  model <- rtmb_model(list(), code, init = list(x = 1), silent = TRUE)
  fit <- model$sample(
    chains = 1, warmup = 250, sampling = 1000,
    seed = 42, delta = 0.95, metric = "diag", progress = "none"
  )
  for (method in c("normal", "warp3")) {
    expect_equal(as.numeric(fit$bridgesampling(method, seed = 42)), 0, tolerance = 0.06)
    previous <- fit$log_ml
    expect_error(fit$bridgesampling(method, seed = 42, max_iter = 1), "did not converge")
    expect_identical(fit$log_ml, previous)
  }

  expect_error(fit$bridgesampling(max_iter = 0), "positive integer")
  short <- fit$clone()
  short$fit <- short$fit[1:2, , , drop = FALSE]
  expect_error(short$bridgesampling(), "at least four retained posterior draws")
  singular <- fit$clone()
  singular$fit[, , 2] <- 1
  expect_error(singular$bridgesampling(), "proposal covariance is singular")
})
