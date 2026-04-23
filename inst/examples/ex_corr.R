\donttest{
  # Simulate bivariate normal data with a true correlation of 0.5
  set.seed(123)
  N <- 50
  rho <- 0.5
  cov_mat <- matrix(c(1, rho, rho, 1), nrow = 2)

  if (requireNamespace("MASS", quietly = TRUE)) {
    data_corr <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = cov_mat)
    colnames(data_corr) <- c("X1", "X2")

    fit_corr <- rtmb_corr(data = data_corr)

    mcmc_corr <- fit_corr$sample(sampling = 500, warmup = 500, chains = 2)
    mcmc_corr$summary()

    bf_corr <- mcmc_corr$bayes_factor(null_model = "corr")
    print(bf_corr)
  }
}
