# RTMB-based GLMM wrapper function

RTMB-based GLMM wrapper function

## Usage

``` r
rtmb_glmer(
  formula,
  data = NULL,
  family = "gaussian",
  laplace = FALSE,
  prior = prior_flat(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  gmc = NULL,
  centering = NULL,
  cwc = NULL,
  std = FALSE,
  view = NULL,
  within = NULL,
  factors = NULL,
  contrasts = "treatment",
  sigma_by = NULL,
  resid_corr = NULL,
  resid_time = NULL,
  resid_group = NULL,
  generate = NULL,
  missing = c("listwise", "fiml"),
  WAIC = FALSE,
  .force_sum = FALSE
)
```

## Arguments

- formula:

  lme4-style formula (e.g., Y ~ X + (1 \| GID))

- data:

  Optional data frame. If omitted, variables are resolved from the
  formula environment. In this mode, use bare variable names; formulas
  using `$`, `[[`, or `.` require an explicit data argument.

- family:

  Character string of the distribution family (e.g., "gaussian",
  "binomial", "poisson", "ordered", "sequential")

- laplace:

  Logical; whether to marginalize random effects using Laplace
  approximation

- prior:

  An object of class "rtmb_prior" specifying the prior distribution. Use
  \`prior_flat()\`, \`prior_normal()\`, \`prior_weak()\`,
  \`prior_rhs()\`, or \`prior_ssp()\`. \`prior_jzs()\` is available for
  continuous families. Default is \`prior_flat()\`.

- y_range:

  Theoretical minimum and maximum values of the response variable as a
  vector c(min, max). Required when using weakly informative or
  regularized priors with continuous models.

- init:

  List of initial values (generated automatically based on glm if
  omitted)

- fixed:

  Optional named list of fixed values for specific parameters.

- gmc:

  Character vector of variable names for Grand Mean Centering (GMC). If
  "all", all numeric variables are centered.

- centering:

  Alias for \`gmc\`.

- cwc:

  List for Centering Within Cluster (CWC). Should contain `cluster`
  (group variable) and `pars` (variable names to center). You can also
  use `cwc = list(ID, "x")` or `cwc = list(ID, "all")`; `"all"` centers
  all numeric fixed-effect variables within the cluster.

- std:

  Logical; if \`TRUE\`, add post-hoc standardized fixed-effect
  coefficients as \`b_std\`. Every column of the fixed-effect design
  matrix is standardized, including columns generated from factors and
  interactions. For Gaussian models, coefficients are scaled by both the
  predictor and response standard deviations. For other families, only
  the predictor standard deviations are used, so coefficients remain on
  the link scale.

- view:

  Character vector of parameter names to prioritize in summary.

- within:

  Optional list for wide-to-long conversion.

- factors:

  Character vector of variable names to be treated as factors.

- contrasts:

  Character string specifying the contrast type ("treatment" or "sum").

- sigma_by:

  Character vector specifying variables to group residual variance by
  (heteroscedasticity).

- resid_corr:

  Residual correlation structure: "ar1" (Autoregressive), "cs" (Compound
  Symmetry), "toep" (Toeplitz), or "un" (Unstructured).

- resid_time:

  Variable name for time points in residual correlation.

- resid_group:

  Variable name for grouping in residual correlation.

- generate:

  Optional expression for generated quantities.

- missing:

  Missing value handling strategy: "listwise".

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

- .force_sum:

  Logical; internal use only.

## Examples

``` r
  # --- 1. Linear Regression (rtmb_lm) ---
  # Fit a linear regression model using the debate dataset
  data(debate, package = "BayesRTMB")
  fit_lm <- rtmb_lm(sat ~ talk + perf, data = debate)
#> Pre-checking model code...
#> Checking RTMB setup...
  map_lm <- fit_lm$optimize()
#> Starting RTMB optimization...
  map_lm$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 395.58
#> Approx. Log Marginal Likelihood (Laplace): -404.61
#> 
#> Point Estimates and 95% Wald CI:
#>  variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept   1.83366     0.21105    1.42000    2.24732 
#> b[talk]     0.28694     0.05291    0.18323    0.39064 
#> b[perf]     0.15632     0.02987    0.09777    0.21486 
#> sigma       0.90453     0.03693    0.83497    0.97988 
#> 

  # --- 2. Generalized Linear Model (rtmb_glm) ---
  # Fit a logistic regression model using the debate dataset
  data(debate, package = "BayesRTMB")
  fit_glm <- rtmb_glm(cond ~ talk + sat, data = debate, family = "bernoulli")
#> Pre-checking model code...
#> Checking RTMB setup...
  map_glm <- fit_glm$optimize()
#> Starting RTMB optimization...
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
#> Intercept  -3.16409     0.57473   -4.29053   -2.03765 
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
  # \donttest{
  mcmc_glmer <- fit_glmer$sample(sampling = 500, warmup = 500, chains = 2)
#> Starting sequential sampling (chains = 2)...
#> chain 1 started...
#> chain 1: iter 200/1000 (20%) warmup
#> chain 1: iter 400/1000 (40%) warmup
#> chain 1: iter 600/1000 (60%) sampling
#> chain 1: iter 800/1000 (80%) sampling
#> chain 1: iter 1000/1000 (100%) sampling
#> chain 1 done (100%)
#> chain 2 started...
#> chain 2: iter 200/1000 (20%) warmup
#> chain 2: iter 400/1000 (40%) warmup
#> chain 2: iter 600/1000 (60%) sampling
#> chain 2: iter 800/1000 (80%) sampling
#> chain 2: iter 1000/1000 (100%) sampling
#> chain 2 done (100%)
#> sampling: 100%
  mcmc_glmer$summary()
#>      variable     mean     sd      map     q2.5    q97.5  ess_bulk  ess_tail  rhat 
#> lp             -513.26  12.98  -512.75  -543.41  -490.50       137        94  1.01 
#> Intercept         2.64   0.09     2.65     2.46     2.81      1147       730  1.00 
#> b[cond]           0.76   0.13     0.77     0.51     0.99      1095       790  1.01 
#> sigma             0.83   0.05     0.81     0.75     0.93       328       368  1.01 
#> sd[group:Int]     0.39   0.10     0.41     0.06     0.55       173        66  1.01 
#> r_re[1]          -0.67   0.78    -0.61    -2.29     0.84      1421       617  1.00 
#> r_re[2]          -0.94   0.79    -1.02    -2.48     0.60       977       714  1.00 
#> r_re[3]           0.03   0.83     0.09    -1.48     1.64      1491       601  1.00 
#> r_re[4]           0.59   0.85     0.67    -1.07     2.32      1181       634  1.01 
#> r_re[5]          -0.36   0.78    -0.42    -1.86     1.21      1726       779  1.00 
  # }

  # --- 4. Linear Mixed Model (rtmb_lmer) ---
  # A convenient wrapper for Gaussian mixed models (identical to rtmb_glmer with family="gaussian")
  fit_lmer <- rtmb_lmer(sat ~ talk + (1 | group), data = debate)
#> Pre-checking model code...
#> Checking RTMB setup...
  map_lmer <- fit_lmer$optimize()
#> Starting RTMB optimization...
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

  # Fit a linear regression using debate predictors with the Horseshoe prior
  fit_rhs <- rtmb_lm(
    sat ~ talk + perf + skill,
    data = debate,
    prior = prior_rhs(),
    y_range = c(1, 5)
  )
#> Pre-checking model code...
#> Checking RTMB setup...
  map_rhs <- fit_rhs$optimize()
#> Starting RTMB optimization...
#> Warning: Best optimization run ended with singular convergence. Estimates may be usable, but check opt_history or try more starts.
#> SE warning: sdreport() returned pdHess = FALSE; Hessian-based fallback will be attempted.
#> SE warning: sdreport() produced non-finite standard errors; Hessian-based fallback will be attempted.
#> SE warning: Hessian matrix was singular; using MASS::ginv() to approximate the covariance matrix.
  # Summarize only the fixed effects (slopes)
  map_rhs$summary("b")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 398.12
#> Approx. Log Marginal Likelihood (Laplace): NA
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[talk]    0.26657     0.05253    0.16361    0.36954 
#> b[perf]    0.15374     0.02944    0.09605    0.21144 
#> b[skill]   0.19129     0.06369    0.06645    0.31613 
#> 

  # Fit a linear regression with the Spike-and-Slab prior
  fit_ssp <- rtmb_lm(
    sat ~ talk + perf + skill,
    data = debate,
    prior = prior_ssp(),
    y_range = c(1, 5)
  )
#> Pre-checking model code...
#> Checking RTMB setup...
  map_ssp <- fit_ssp$optimize()
#> Starting RTMB optimization...
#> SE warning: sdreport() produced non-finite standard errors; Hessian-based fallback will be attempted.
#> SE warning: Hessian matrix was singular; using MASS::ginv() to approximate the covariance matrix.
  map_ssp$summary("b")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 400.19
#> Approx. Log Marginal Likelihood (Laplace): NA
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[talk]    0.00000     0.00000    0.00000    0.00000 
#> b[perf]    0.16011     0.03082    0.09970    0.22052 
#> b[skill]   0.22134     0.06679    0.09043    0.35225 
#> 
```
