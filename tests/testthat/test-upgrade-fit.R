test_that("upgrade_fit refreshes MCMC fit objects", {
  mdl <- rtmb_lm(mpg ~ wt, data = mtcars)
  draws <- array(
    c(0, 1, 2, 3),
    dim = c(2, 1, 2),
    dimnames = list(NULL, "chain1", c("lp", "b"))
  )

  fit <- MCMC_Fit$new(
    model = mdl,
    fit = draws,
    random_fit = NULL,
    eps = c(chain1 = 0.1),
    accept = c(chain1 = 0.9),
    treedepth = c(chain1 = 2),
    laplace = FALSE,
    posterior_mean = c(b = 2)
  )
  fit$transform_dims <- list(extra = 1)

  upgraded <- upgrade_fit(fit, upgrade_model = FALSE)

  expect_true(inherits(upgraded, "mcmc_fit"))
  expect_true(identical(upgraded$model, fit$model))
  expect_identical(upgraded$fit, fit$fit)
  expect_identical(upgraded$transform_dims, fit$transform_dims)
  expect_true("EAP" %in% names(upgraded))
})

test_that("upgrade_fit refreshes VB fit objects", {
  draws <- array(
    NA_real_,
    dim = c(3, 2, 2),
    dimnames = list(NULL, paste0("est", 1:2), c("lp", "theta"))
  )
  draws[, 1, "lp"] <- c(0, 1, 2)
  draws[, 2, "lp"] <- c(10, 11, 12)
  draws[, 1, "theta"] <- c(100, 200, 300)
  draws[, 2, "theta"] <- c(1, 2, 3)

  model <- list(par_list = list(theta = list(dim = 1)))
  fit <- VB_Fit$new(
    model = model,
    fit = draws,
    random_fit = NULL,
    elbo_history = list(c(1), c(2)),
    laplace = FALSE,
    posterior_mean = c(theta = 2),
    ELBO = c(1, 2),
    rel_obj_vals = c(1e-4, 1e-5),
    best_chain = 2,
    mu_history = NULL
  )
  fit$generate_dims <- list(g = 1)

  upgraded <- upgrade_fit(fit)

  expect_true(inherits(upgraded, "advi_fit"))
  expect_identical(upgraded$best_chain, 2)
  expect_identical(upgraded$generate_dims, fit$generate_dims)
  expect_equal(unname(upgraded$EAP("theta")), 2)
})

test_that("upgrade_fit refreshes MAP and model objects", {
  mdl <- rtmb_lm(mpg ~ wt, data = mtcars)
  fit <- mdl$optimize(se_method = "none")

  upgraded <- upgrade_fit(fit)

  expect_true(inherits(upgraded, "map_fit"))
  expect_true(inherits(upgraded$model, "RTMB_Model"))
  expect_false(identical(upgraded$model, fit$model))
  expect_identical(names(upgraded$par), names(fit$par))
})

test_that("upgrade_fit refreshes classic fit objects and preserves extras", {
  mdl <- rtmb_lm(mpg ~ wt, data = mtcars)
  fit <- mdl$classic()
  fit$se_method <- "custom"
  fit$cluster <- "id"

  upgraded <- upgrade_fit(fit, upgrade_model = FALSE)

  expect_true(inherits(upgraded, "Classic_Fit"))
  expect_true(identical(upgraded$model, fit$model))
  expect_identical(upgraded$se_method, "custom")
  expect_identical(upgraded$cluster, "id")
  expect_identical(names(upgraded$par), names(fit$par))
})
