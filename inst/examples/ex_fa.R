
  # Prepare a subset of the BigFive dataset for factor analysis
  data(BigFive, package = "BayesRTMB")
  fa_data <- BigFive[, 1:10]

  # --- 1. Standard Exploratory Factor Analysis (1 Factor) ---
  fit_fa1 <- rtmb_fa(data = fa_data, nfactors = 1)

  # Maximum A Posteriori (MAP) estimation
  map_fa1 <- fit_fa1$optimize()
  map_fa1$summary()
