
  # Simulate mediation data
  set.seed(123)
  N <- 100
  X <- rnorm(N)
  # M is influenced by X
  M <- 0.5 * X + rnorm(N, 0, 0.5)
  # Y is influenced by both X and M
  Y <- 0.3 * X + 0.8 * M + rnorm(N, 0, 0.5)
  dat <- data.frame(X, M, Y)

  # Fit a mediation model
  # The formula list specifies the two regression equations
  fit_med <- rtmb_mediation(list(M ~ X, Y ~ X + M), data = dat)
  
  # Maximum A Posteriori (MAP) estimation
  map_med <- fit_med$optimize()
  # The summary automatically calculates indirect, direct, and total effects
  map_med$summary()
