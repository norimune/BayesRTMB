# RTMB-based GLM wrapper function (no random effects)

RTMB-based GLM wrapper function (no random effects)

## Usage

``` r
rtmb_glm(
  formula,
  data,
  family = "gaussian",
  prior = prior_flat(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  gmc = NULL,
  view = NULL,
  factors = NULL,
  contrasts = "treatment",
  missing = c("listwise", "fiml"),
  ...
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

- prior:

  An object of class \`"rtmb_prior"\`. Use \`prior_flat()\` for no
  prior, \`prior_normal()\` for default normal/exponential priors, or
  \`prior_weak()\`, \`prior_rhs()\`, \`prior_ssp()\` for weakly
  informative or regularized Bayesian inference. Default is
  \`prior_flat()\`.

- y_range:

  Theoretical minimum and maximum values of the response variable as a
  vector c(min, max). Specifying this automatically enables weakly
  informative priors.

- init:

  List of initial values

- fixed:

  Optional named list of fixed values for specific parameters.

- gmc:

  Character vector of variable names for GMC

- view:

  Optional character vector of parameter names to show first in summary.

- factors:

  Character vector of variable names to be treated as factors.

- contrasts:

  Character string specifying the contrast type ("treatment" or "sum").

- missing:

  Missing value handling strategy: "listwise".

- ...:

  Additional arguments passed to
  [`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md).

## Examples

``` r
  # --- 1. Linear Regression (rtmb_lm) ---
  # Fit a linear regression model using the mtcars dataset
  fit_lm <- rtmb_lm(mpg ~ wt + cyl, data = mtcars)
#> Pre-checking model code...
#> Checking RTMB setup...
  map_lm <- fit_lm$optimize()
#> Starting RTMB optimization...
#> 
  map_lm$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 74.01
#> Approx. Log Marginal Likelihood (Laplace): -74.09
#> 
#> Point Estimates and 95% Wald CI:
#>  variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept  39.68626     1.63262   36.48639   42.88613 
#> b[wt]      -3.19097     0.72055   -4.60323   -1.77871 
#> b[cyl]     -1.50779     0.39477   -2.28153   -0.73406 
#> sigma       2.44420     0.30553    1.91310    3.12275 
#> 

  # --- 2. Generalized Linear Model (rtmb_glm) ---
  # Fit a logistic regression model using the debate dataset
  data(debate, package = "BayesRTMB")
  fit_glm <- rtmb_glm(cond ~ talk + sat, data = debate, family = "bernoulli")
#> Pre-checking model code...
#> Checking RTMB setup...
  map_glm <- fit_glm$optimize()
#> Starting RTMB optimization...
#> 
  map_glm$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 183.31
#> Approx. Log Marginal Likelihood (Laplace): -186.56
#> 
#> Point Estimates and 95% Wald CI:
#>  variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept  -3.16409     0.57473   -4.29054   -2.03765 
#> b[talk]     0.84795     0.15040    0.55317    1.14273 
#> b[sat]      0.17405     0.13431   -0.08920    0.43730 
#> 

  # --- 3. Generalized Linear Mixed Model (rtmb_glmer) ---
  # Fit a linear mixed-effects model using the debate dataset
  data(debate, package = "BayesRTMB")
  fit_glmer <- rtmb_glmer(talk ~ cond + (1 | group), data = debate, family = "gaussian")
#> Pre-checking model code...
#> Checking RTMB setup...

  # MAP estimation using Laplace approximation for random effects
  map_glmer <- fit_glmer$optimize(laplace = TRUE)
#> Starting RTMB optimization...
#> 
  map_glmer$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 393.51
#> Approx. Log Marginal Likelihood (Laplace): -400.58
#> Note: Random effects are stored in $random_effects (use ranef = TRUE to show them)
#> 
#> Point Estimates and 95% Wald CI:
#>      variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept       2.64000     0.08790    2.46772    2.81228 
#> b[cond]         0.76000     0.12431    0.51636    1.00364 
#> sigma           0.82057     0.04103    0.74397    0.90506 
#> sd[group:Int]   0.40233     0.07340    0.28137    0.57528 
#> 

  # MCMC sampling (chains and iterations reduced for faster execution)
  if (FALSE) { # \dontrun{
  mcmc_glmer <- fit_glmer$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_glmer$summary()
  } # }

  # --- 4. Linear Mixed Model (rtmb_lmer) ---
  # A convenient wrapper for Gaussian mixed models (identical to rtmb_glmer with family="gaussian")
  fit_lmer <- rtmb_lmer(sat ~ talk + (1 | group), data = debate)
#> Pre-checking model code...
#> Checking RTMB setup...
  map_lmer <- fit_lmer$optimize()
#> Starting RTMB optimization...
#> 
  map_lmer$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 394.85
#> Approx. Log Marginal Likelihood (Laplace): -402.71
#> Note: Random effects are stored in $random_effects (use ranef = TRUE to show them)
#> 
#> Point Estimates and 95% Wald CI:
#>      variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept       2.60860     0.18122    2.25341    2.96379 
#> b[talk]         0.27309     0.05535    0.16461    0.38157 
#> sigma           0.77837     0.03900    0.70557    0.85869 
#> sd[group:Int]   0.53681     0.06762    0.41936    0.68715 
#> 

  # --- 5. Regularized Regression (Variable Selection) ---
  # You can apply regularization to the fixed effects to shrink noise variables towards zero.
  # Use prior = prior_rhs() for the Regularized Horseshoe prior,
  # or prior_ssp() for the Spike-and-Slab prior.
  # Note: When using regularization, you must specify 'y_range' (the theoretical minimum and maximum
  # values of the response variable) to automatically set up the required weakly informative priors.

  # Fit a linear regression using all predictors in mtcars with the Horseshoe prior
  # 'mpg' theoretically ranges roughly between 0 and 40
  fit_rhs <- rtmb_lm(mpg ~ ., data = mtcars, prior = prior_rhs(), y_range = c(0, 40))
#> Pre-checking model code...
#> Checking RTMB setup...
  map_rhs <- fit_rhs$optimize()
#> Starting RTMB optimization...
#> 
#> Warning: Optimization did not converge ( convergence code = 1; message = function evaluation limit reached without convergence (9)). Estimates may be unreliable; consider increasing num_estimate, changing initial values, or adjusting optimizer control settings.
#> SE warning from sdreport(): NaNs produced
#> SE warning: sdreport() returned pdHess = FALSE; Hessian-based fallback will be attempted.
  # Summarize only the fixed effects (slopes)
  map_rhs$summary("b")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 95.77
#> Approx. Log Marginal Likelihood (Laplace): NA
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[cyl]    -0.43391     0.74997   -1.90383    1.03601 
#> b[disp]   -0.00328     0.01159   -0.02600    0.01945 
#> b[hp]     -0.00000     0.00014   -0.00027    0.00027 
#> b[drat]    0.87950     1.13178   -1.33874    3.09775 
#> b[wt]     -2.29835     1.28136   -4.80978    0.21307 
#> b[qsec]    0.49484     0.55340   -0.58980    1.57948 
#> b[vs]     -0.09523     1.30357   -2.65017    2.45972 
#> b[am]      1.66383     1.34271   -0.96784    4.29550 
#> b[gear]    0.71143     1.03019   -1.30770    2.73055 
#> b[carb]   -0.84622     0.49989   -1.82599    0.13355 
#> 

  # Fit a linear regression with the Spike-and-Slab prior
  fit_ssp <- rtmb_lm(mpg ~ ., data = mtcars, prior = prior_ssp(), y_range = c(0, 40))
#> Pre-checking model code...
#> Checking RTMB setup...
  map_ssp <- fit_ssp$optimize()
#> Starting RTMB optimization...
#> 
#> Warning: Optimization did not converge ( convergence code = 1; message = singular convergence (7)). Estimates may be unreliable; consider increasing num_estimate, changing initial values, or adjusting optimizer control settings.
#> SE warning: sdreport() produced non-finite standard errors; Hessian-based fallback will be attempted.
#> SE warning: Hessian matrix was singular; using MASS::ginv() to approximate the covariance matrix.
  map_ssp$summary("b")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 79.98
#> Approx. Log Marginal Likelihood (Laplace): NA
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[cyl]    -2.84120     0.31331   -3.45528   -2.22713 
#> b[disp]    0.00000     0.00000    0.00000    0.00000 
#> b[hp]      0.00000     0.00000    0.00000    0.00000 
#> b[drat]    0.00000     0.00000    0.00000    0.00000 
#> b[wt]      0.00000     0.00000    0.00000    0.00000 
#> b[qsec]    0.00000     0.00000    0.00000    0.00000 
#> b[vs]      0.00000     0.00000   -0.00000    0.00000 
#> b[am]      0.00000     0.00000    0.00000    0.00000 
#> b[gear]    0.00000     0.00000    0.00000    0.00000 
#> b[carb]    0.00000     0.00000    0.00000    0.00000 
#> 

  # For models with complex penalties, MCMC often provides more reliable credible intervals
  if (FALSE) { # \dontrun{
  mcmc_ssp <- fit_ssp$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_ssp$summary("b")
  mcmc_ssp$summary("b")
  } # }
```
