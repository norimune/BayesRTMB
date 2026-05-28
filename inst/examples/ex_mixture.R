
  # Simulate 1D mixture data (2 components)
  set.seed(123)
  N <- 100
  group <- rbinom(N, 1, 0.4)
  y <- ifelse(group == 1, rnorm(N, 5, 1), rnorm(N, 0, 1))
  dat <- data.frame(y)

  # Fit a 1D Gaussian mixture model with 2 components
  fit_mix <- rtmb_mixture(y ~ 1, k = 2, data = dat)
  
  # MAP estimation
  map_mix <- fit_mix$optimize()
  map_mix$summary()
