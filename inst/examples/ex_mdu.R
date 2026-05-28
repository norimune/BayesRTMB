
  # Simulate rating data for Multidimensional Unfolding (MDU)
  set.seed(123)
  N <- 50  # Number of persons
  M <- 10  # Number of items
  D <- 2   # Number of dimensions
  
  # True person and item coordinates in a 2D space
  theta <- matrix(rnorm(N * D), N, D)
  delta <- matrix(rnorm(M * D), M, D)
  
  # Generate distance-like ratings (smaller means more preferred)
  Y <- matrix(NA, N, M)
  for(i in 1:N) {
    for(j in 1:M) {
      Y[i, j] <- sum((theta[i,] - delta[j,])^2) + rnorm(1, 0, 0.5)
    }
  }

  # Fit a 2-dimensional MDU model
  fit_mdu <- rtmb_mdu(Y, ndim = 2, distance = "squared", method = "rating")
  
  # MAP estimation
  map_mdu <- fit_mdu$optimize()
  map_mdu$summary()
  
  # Note: MDU models have many parameters, so MCMC sampling might take time.
  \donttest{
  # mcmc_mdu <- fit_mdu$sample(sampling = 500, warmup = 500, chains = 2)
  }
