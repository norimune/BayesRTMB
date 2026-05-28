
  # Prepare a subset of the BigFive dataset for factor analysis
  data(BigFive, package = "BayesRTMB")
  fa_data <- BigFive[, 1:10]

  # --- 1. Standard Exploratory Factor Analysis (1 Factor) ---
  fit_fa1 <- rtmb_fa(data = fa_data, nfactors = 1)

  # Maximum A Posteriori (MAP) estimation
  map_fa1 <- fit_fa1$optimize()
  map_fa1$summary()



  # --- 2. Factor Analysis with Rotation and Factor Scores (2 Factors) ---
  # Extract 2 factors, apply Promax rotation during model fitting, and calculate factor scores
  fit_fa2 <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "promax", score = TRUE)

  # Setting se_method = "sampling" enables standard errors and 95% CIs
  # for transformed and generated quantities, such as factor scores and post-hoc rotations.
  # (It uses multivariate normal sampling from the unconstrained parameter space)
  map_fa2 <- fit_fa2$optimize(se_method = "sampling")
  map_fa2$summary()

  # Post-hoc rotation using the fa_rotate() method
  # You can also apply other rotation methods (e.g., "varimax" from the stats package)
  # to the unrotated loading matrix ("L") after estimation.
  map_fa2$fa_rotate(target = "L", rotate = "varimax")

  # The post-hoc rotated loadings are automatically stored with the method's
  # suffix (e.g., "L_varimax")
  map_fa2$summary("L_varimax")



  # --- 3. Regularized Factor Analysis using Spike-and-Slab Prior (SSP) ---
  # Specifying 'rotate = "ssp"' enables sparse loading matrix estimation
  fit_ssp <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "ssp")

  map_ssp <- fit_ssp$optimize()
  map_ssp$summary()
