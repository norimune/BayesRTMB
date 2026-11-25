test_that("simplex transformations are inverses", {
  for (theta in list(c(0.3, 0.7), c(0.2, 0.3, 0.5), c(0.1, 0.2, 0.3, 0.4))) {
    par_list <- list(theta = Dim(length(theta), type = "simplex"))
    unconstrained <- to_unconstrained(list(theta = theta), par_list)

    expect_equal(to_constrained(unconstrained, par_list)$theta, theta,
                 tolerance = 1e-12)
  }
})

test_that("bridge target reproduces saved log densities for a simplex", {
  code <- rtmb_code(
    parameters = {
      theta <- Dim(3, type = "simplex")
    },
    model = {
      theta ~ dirichlet(c(2, 3, 4))
    }
  )
  model <- rtmb_model(
    data = list(),
    code = code,
    init = list(theta = c(0.2, 0.3, 0.5)),
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
  expect_equal(apply(draws_uc, 1, fit$log_prob()),
               as.numeric(fit$fit[, , 1]), tolerance = 1e-10)
})
