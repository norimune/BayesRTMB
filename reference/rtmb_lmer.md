# RTMB-based Linear Mixed Model (LMM) wrapper function

RTMB-based Linear Mixed Model (LMM) wrapper function

## Usage

``` r
rtmb_lmer(
  formula,
  data,
  laplace = TRUE,
  prior = prior_uniform(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  null = NULL,
  gmc = NULL,
  cwc = NULL,
  view = NULL,
  sigma_by = NULL,
  factors = NULL,
  contrasts = "treatment",
  resid_corr = NULL,
  resid_time = NULL,
  resid_group = NULL,
  within = NULL,
  ...
)
```

## Arguments

- formula:

  Formula

- data:

  Data frame

- laplace:

  Logical; whether to marginalize random effects using Laplace
  approximation

- prior:

  An object of class "rtmb_prior" specifying the prior distribution.
  Default is NULL (flat prior).

- y_range:

  Theoretical minimum and maximum values of the response variable

- init:

  Initial values

- fixed:

  Optional named list of fixed values for specific parameters.

- null:

  Null model parameters

- gmc:

  Character vector of variable names for GMC

- cwc:

  List for CWC

- view:

  Character vector of parameter names to prioritize in summary.

- sigma_by:

  Character vector specifying variables to group residual variance by
  (heteroscedasticity).

- factors:

  Character vector of variable names to be treated as factors.

- contrasts:

  Character string specifying the contrast type ("treatment" or "sum").

- resid_corr:

  Residual correlation structure (e.g., "ar1", "cs", "un", "toep").

- resid_time:

  Variable name for time points in residual correlation.

- resid_group:

  Variable name for grouping in residual correlation.

- within:

  Optional list for wide-to-long conversion. For repeated measures data
  in wide format, specify the factor names and their levels, e.g.,
  `list(Time = 4)` or `list(A = 2, B = 3)`. The total number of levels
  must match the number of columns in
  [`cbind()`](https://rdrr.io/r/base/cbind.html) on the LHS. If omitted
  and the LHS is [`cbind()`](https://rdrr.io/r/base/cbind.html), the
  within-factor name is inferred from RHS variables not present in the
  data.

- ...:

  Additional arguments passed to
  [`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md).

## Value

RTMB_Model object

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
