\donttest{
  # Simulate two-sample data with a true effect size
  set.seed(123)
  y1 <- rnorm(30, mean = 0.5, sd = 1)
  y2 <- rnorm(30, mean = 0.0, sd = 1)

  # Fit the Bayesian two-sample t-test model
  # r = 0.707 is the standard scale for the Cauchy prior on the effect size
  fit_ttest <- rtmb_ttest(y1, y2, r = 0.707)

  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_ttest <- fit_ttest$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_ttest$summary()

  # Calculate Bayes factor against the null hypothesis (effect size delta = 0)
  # Specifying "delta" automatically fixes the parameter to 0 and drops its prior
  bf_ttest <- mcmc_ttest$bayes_factor(null_model = "delta")
  print(bf_ttest)
}
