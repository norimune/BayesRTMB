 \donttest{
  # Prepare a subset of the mtcars dataset for factor analysis
  # Scaling is recommended for variables with different units
  fa_data <- scale(mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")])

  # --- 1. Standard Exploratory Factor Analysis (1 Factor) ---
  fit_fa1 <- rtmb_fa(data = fa_data, nfactors = 1)

  # Maximum A Posteriori (MAP) estimation
  map_fa1 <- fit_fa1$optimize()
  map_fa1$summary()

  # --- 2. Factor Analysis with Rotation and Factor Scores (2 Factors) ---
  # Extract 2 factors, apply Promax rotation, and calculate factor scores
  fit_fa2 <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "promax", score = TRUE)

  map_fa2 <- fit_fa2$optimize()
  # The summary prioritizes rotated loadings (L_promax), standard deviations, and factor correlations
  map_fa2$summary()

  # --- 3. Regularized Factor Analysis using Spike-and-Slab Prior (SSP) ---
  # Specifying 'rotate = "ssp"' enables sparse loading matrix estimation
  fit_ssp <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "ssp")

  map_ssp <- fit_ssp$optimize()
  map_ssp$summary()

  # MCMC sampling for the SSP model (chains and iterations reduced for faster execution)
  mcmc_ssp <- fit_ssp$sample(sampling = 500, warmup = 500, chains = 2)

  # --- 4. Resolving Label Switching in MCMC ---
  # In Bayesian factor analysis, MCMC chains often suffer from label switching or sign flipping.
  # This can be resolved by applying post-hoc Procrustes rotation to the posterior samples.

  # Summary of unrotated loadings (may show poor convergence / large SE due to switching)
  mcmc_ssp$summary("L")

  # Apply Procrustes rotation targeting the loading matrix "L"
  mcmc_ssp$rotate(target = "L")

  # Summary of the rotated loadings (L_rot) with stabilized estimates
  mcmc_ssp$summary("L_rot")
}
