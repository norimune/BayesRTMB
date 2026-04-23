\donttest{
  # Fit a linear regression model using the built-in mtcars dataset
  fit_model <- rtmb_lm(mpg ~ wt + cyl, data = mtcars)

  # 1. Maximum A Posteriori (MAP) estimation
  map_res <- fit_model$optimize()

  # Output the summary of the MAP estimates
  map_res$summary()

  # 2. MCMC sampling
  # The number of chains and iterations are reduced to minimize execution time
  mcmc_res <- fit_model$sample(sampling = 500, warmup = 500, chains = 2)

  # Output the summary of the posterior distribution
  mcmc_res$summary()

  # 3. Variational Inference (ADVI)
  vb_res <- fit_model$variational()

  # Output the summary of the variational approximation
  vb_res$summary()
}
