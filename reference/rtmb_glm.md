# RTMB-based GLM wrapper function (no random effects)

RTMB-based GLM wrapper function (no random effects)

## Usage

``` r
rtmb_glm(
  formula,
  data,
  family = "gaussian",
  penalty = c("none", "rhs", "ssp"),
  y_range = NULL,
  use_weak_info = FALSE,
  prior = list(Intercept_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, nu_rate = 0.1,
    cutpoint_sd = 2.5, shape_rate = 1, phi_rate = 1, lkj_eta = 1),
  weak_info_prior = list(max_beta = 1, sd_ratio = 0.5, expected_vars = 3, slab_scale = 2,
    slab_df = 4, ssp_ratio = 0.25),
  init = NULL,
  null = NULL
)
```

## Arguments

- formula:

  Formula

- data:

  Data frame

- family:

  Character string of the distribution family (e.g., "gaussian",
  "binomial", "poisson")

- penalty:

  Type of regularization for fixed effects: "none", "rhs" (Regularized
  Horseshoe), or "ssp" (Spike and Slab Prior). Default is "none".

- y_range:

  Theoretical minimum and maximum values of the response variable as a
  vector c(min, max). Specifying this automatically enables weakly
  informative priors.

- use_weak_info:

  Logical; whether to explicitly use weakly informative priors (requires
  y_range for continuous models).

- prior:

  List of hyperparameters for the default fixed priors.

- weak_info_prior:

  List of hyperparameters for the weakly informative priors and
  regularization.

- init:

  List of initial values

- null:

  Character string specifying the target parameter for the null model.

## Examples

``` r
if (FALSE) { # \dontrun{
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

  # --- 4. Linear Mixed Model (rtmb_lmer) ---
  # A convenient wrapper for Gaussian mixed models (identical to rtmb_glmer with family="gaussian")
  fit_lmer <- rtmb_lmer(mpg ~ wt + (1 | cyl), data = mtcars)
  map_lmer <- fit_lmer$optimize()
  map_lmer$summary()

  # --- 5. Regularized Regression (Variable Selection) ---
  # You can apply regularization to the fixed effects to shrink noise variables towards zero.
  # Use penalty = "rhs" for the Regularized Horseshoe prior, or "ssp" for the Spike-and-Slab prior.
  # Note: When using regularization, you must specify 'y_range' (the theoretical minimum and maximum
  # values of the response variable) to automatically set up the required weakly informative priors.

  # Fit a linear regression using all predictors in mtcars with the Horseshoe prior
  # 'mpg' theoretically ranges roughly between 0 and 40
  fit_rhs <- rtmb_lm(mpg ~ ., data = mtcars, penalty = "rhs", y_range = c(0, 40))
  map_rhs <- fit_rhs$optimize()
  # Summarize only the fixed effects (slopes)
  map_rhs$summary("b")

  # Fit a linear regression with the Spike-and-Slab prior
  fit_ssp <- rtmb_lm(mpg ~ ., data = mtcars, penalty = "ssp", y_range = c(0, 40))
  map_ssp <- fit_ssp$optimize()
  map_ssp$summary("b")

  # For models with complex penalties, MCMC often provides more reliable credible intervals
  mcmc_ssp <- fit_ssp$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_ssp$summary("b")
} # }
```
