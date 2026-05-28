
  # Estimate the correlation between two variables in the debate dataset
  data(debate, package = "BayesRTMB")

  fit_corr <- rtmb_corr(cbind(sat, perf), data = debate)

  \donttest{
  mcmc_corr <- fit_corr$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_corr$summary()

  bf_corr <- mcmc_corr$bayes_factor(fixed = list(corr = 0))
  print(bf_corr)
  }
