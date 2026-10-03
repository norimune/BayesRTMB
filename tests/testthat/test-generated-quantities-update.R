make_generated_quantity_test_fit <- function(type = c("mcmc", "vb")) {
  type <- match.arg(type)
  model <- list(
    data = list(),
    par_list = list(mu = list(length = 1, dim = 1, bounds = "none")),
    par_names = list(),
    transform = NULL
  )
  draws <- array(
    c(-1, -2, 1, 2),
    dim = c(2, 1, 2),
    dimnames = list(NULL, "chain1", c("lp", "mu"))
  )

  if (type == "mcmc") {
    return(MCMC_Fit$new(
      model = model,
      fit = draws,
      random_fit = NULL,
      eps = NULL,
      accept = NULL,
      treedepth = NULL,
      laplace = FALSE,
      posterior_mean = NULL
    ))
  }

  VB_Fit$new(
    model = model,
    fit = draws,
    random_fit = NULL,
    elbo_history = list(),
    laplace = FALSE,
    posterior_mean = NULL,
    ELBO = 0,
    rel_obj_vals = 0,
    best_chain = 1,
    mu_history = NULL
  )
}

for (fit_type in c("mcmc", "vb")) {
  test_that(paste("generated_quantities replaces existing variables for", fit_type), {
    fit <- make_generated_quantity_test_fit(fit_type)

    fit$generated_quantities({
      existing <- mu + 1
      keep <- mu + 10
    }, progress = "none", tape = "none")

    fit$generated_quantities({
      existing <- c(mu + 2, mu + 3)
      added <- keep + 1
    }, progress = "none", tape = "none")

    expect_identical(
      dimnames(fit$generate_fit)[[3]],
      c("keep", "existing[1]", "existing[2]", "added")
    )
    expect_false(anyDuplicated(dimnames(fit$generate_fit)[[3]]) > 0)
    expect_equal(as.numeric(fit$generate_fit[, 1, "keep"]), c(11, 12))
    expect_equal(as.numeric(fit$generate_fit[, 1, "existing[1]"]), c(3, 4))
    expect_equal(as.numeric(fit$generate_fit[, 1, "existing[2]"]), c(4, 5))
    expect_equal(as.numeric(fit$generate_fit[, 1, "added"]), c(12, 13))
    expect_identical(fit$generate_dims$existing, 2L)
  })
}
