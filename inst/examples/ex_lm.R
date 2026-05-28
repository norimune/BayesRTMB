  # --- 1. Linear Regression (rtmb_lm) ---
  # Fit a linear regression model using the debate dataset
  data(debate, package = "BayesRTMB")
  fit_lm <- rtmb_lm(sat ~ talk + perf, data = debate)
  map_lm <- fit_lm$optimize()
  map_lm$summary()

  # --- 2. Generalized Linear Model (rtmb_glm) ---
  # Fit a logistic regression model using the debate dataset
  data(debate, package = "BayesRTMB")
  fit_glm <- rtmb_glm(cond ~ talk + sat, data = debate, family = "bernoulli")
  map_glm <- fit_glm$optimize()
  map_glm$summary()

  # --- 3. Generalized Linear Mixed Model (rtmb_glmer) ---
  # Fit a linear mixed-effects model using the debate dataset
  data(debate, package = "BayesRTMB")
  fit_glmer <- rtmb_glmer(talk ~ cond + (1 | group), data = debate, family = "gaussian")

  # MAP estimation using Laplace approximation for random effects
  map_glmer <- fit_glmer$optimize(laplace = TRUE)
  map_glmer$summary()

  # MCMC sampling (chains and iterations reduced for faster execution)
  \donttest{
  mcmc_glmer <- fit_glmer$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_glmer$summary()
  }

  # --- 4. Linear Mixed Model (rtmb_lmer) ---
  # A convenient wrapper for Gaussian mixed models (identical to rtmb_glmer with family="gaussian")
  fit_lmer <- rtmb_lmer(sat ~ talk + (1 | group), data = debate)
  map_lmer <- fit_lmer$optimize()
  map_lmer$summary()

  # --- 5. Regularized Regression (Variable Selection) ---
  # You can apply regularization to the fixed effects to shrink noise variables towards zero.
  # Use prior = prior_rhs() for the Regularized Horseshoe prior,
  # or prior_ssp() for the Spike-and-Slab prior.
  # Note: When using regularization, you must specify 'y_range' (the theoretical minimum and maximum
  # values of the response variable) to automatically set up the required weakly informative priors.

  # Fit a linear regression using debate predictors with the Horseshoe prior
  fit_rhs <- rtmb_lm(
    sat ~ talk + perf + skill,
    data = debate,
    prior = prior_rhs(),
    y_range = c(1, 5)
  )
  map_rhs <- fit_rhs$optimize()
  # Summarize only the fixed effects (slopes)
  map_rhs$summary("b")

  # Fit a linear regression with the Spike-and-Slab prior
  fit_ssp <- rtmb_lm(
    sat ~ talk + perf + skill,
    data = debate,
    prior = prior_ssp(),
    y_range = c(1, 5)
  )
  map_ssp <- fit_ssp$optimize()
  map_ssp$summary("b")
