test_that("VB point estimates use the best estimate by default", {
  fit_array <- array(
    NA_real_,
    dim = c(3, 2, 2),
    dimnames = list(
      iteration = NULL,
      chain = paste0("est", 1:2),
      variable = c("lp", "theta")
    )
  )
  fit_array[, 1, "lp"] <- c(98, 99, 100)
  fit_array[, 2, "lp"] <- c(-3, -2, -1)
  fit_array[, 1, "theta"] <- c(100, 200, 300)
  fit_array[, 2, "theta"] <- c(1, 2, 3)

  model <- list(
    par_list = list(theta = list(dim = 1)),
    par_names = list(),
    data = list(),
    view = NULL
  )

  vb_fit <- VB_Fit$new(
    model = model,
    fit = fit_array,
    random_fit = NULL,
    elbo_history = list(c(0), c(1)),
    laplace = FALSE,
    posterior_mean = c(theta = 2),
    ELBO = c(0, 10),
    rel_obj_vals = c(1e-4, 1e-5),
    best_chain = 2,
    mu_history = NULL
  )

  expect_equal(unname(vb_fit$EAP("theta")), 2)
  expect_equal(unname(vb_fit$estimate("theta", type = "EAP", drop = TRUE)), 2)
  expect_equal(unname(vb_fit$MAP("theta", type = "joint")), 3)
  expect_equal(unname(vb_fit$get_point_estimate("theta")), 3)

  expect_equal(unname(vb_fit$EAP("theta", chains = 1:2)), 101)
  expect_equal(unname(vb_fit$MAP("theta", chains = 1:2, type = "joint")), 300)
  expect_true(is.list(vb_fit$EAP("theta", drop = FALSE)))
})
