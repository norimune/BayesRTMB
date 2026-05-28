
  # Simulate two-sample data with a true effect size
  set.seed(123)
  y1 <- rnorm(30, mean = 0.5, sd = 1)
  y2 <- rnorm(30, mean = 0.0, sd = 1)

  # Fit the Bayesian two-sample t-test model
  # r = 0.707 is the standard scale for the Cauchy prior on the effect size
  fit_ttest <- rtmb_ttest(y1, y2, prior = prior_jzs(r = 0.707))

  \donttest{
  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_ttest <- fit_ttest$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_ttest$summary()
  }
