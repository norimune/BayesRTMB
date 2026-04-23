\donttest{
  # --- 1. Linear Regression (rtmb_lm) ---
  # Fit a linear regression model using the mtcars dataset
  fit_lm <- rtmb_lm(mpg ~ wt + cyl, data = mtcars)
  map_lm <- fit_lm$optimize()
  map_lm$summary()

  # --- 2. Generalized Linear Model (rtmb_glm) ---
  # Fit a logistic regression model
  fit_glm <- rtmb_glm(am ~ mpg + hp, data = mtcars, family = "binomial")
  map_glm <- fit_glm$optimize()
  map_glm$summary()

  # --- 3. Generalized Linear Mixed Model (rtmb_glmer) ---
  # Fit a linear mixed-effects model with random intercepts
  fit_glmer <- rtmb_glmer(weight ~ Time + (1 | Chick), data = ChickWeight, family = "gaussian")

  # MAP estimation using Laplace approximation for random effects
  map_glmer <- fit_glmer$optimize(laplace = TRUE)
  map_glmer$summary()

  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_glmer <- fit_glmer$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_glmer$summary()
}
