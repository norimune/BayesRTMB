# RTMB-based Linear Mixed Model (LMM) wrapper function

RTMB-based Linear Mixed Model (LMM) wrapper function

## Usage

``` r
rtmb_lmer(
  formula,
  data,
  laplace = TRUE,
  prior = prior_flat(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  gmc = NULL,
  centering = NULL,
  cwc = NULL,
  view = NULL,
  sigma_by = NULL,
  factors = NULL,
  contrasts = "treatment",
  resid_corr = NULL,
  resid_time = NULL,
  resid_group = NULL,
  within = NULL,
  missing = c("listwise", "fiml"),
  WAIC = FALSE,
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

  An object of class \`"rtmb_prior"\`. Use \`prior_flat()\` for no
  prior, \`prior_normal()\` for default normal/exponential priors,
  \`prior_jzs()\` for JZS regression priors, or \`prior_weak()\`,
  \`prior_rhs()\`, \`prior_ssp()\` for weakly informative or regularized
  Bayesian inference. Default is \`prior_flat()\`.

- y_range:

  Theoretical minimum and maximum values of the response variable

- init:

  Initial values

- fixed:

  Optional named list of fixed values for specific parameters.

- gmc:

  Character vector of variable names for GMC

- centering:

  Alias for \`gmc\`.

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

- missing:

  Missing value handling strategy: "listwise".

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

- ...:

  Additional arguments passed to
  [`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md).

## Value

RTMB_Model object

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
#> 
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
#> lp             -510.23  12.54  -507.43  -535.93  -488.01       148       351  1.01 
#> Intercept         2.64   0.09     2.62     2.47     2.81      1162       665  1.00 
#> b[cond]           0.76   0.13     0.75     0.51     1.01      1232       525  1.00 
#> sigma             0.83   0.04     0.82     0.75     0.91       538       648  1.00 
#> sd[group:Int]     0.41   0.08     0.40     0.24     0.57       209       256  1.01 
#> r_re[1]          -0.75   0.77    -0.79    -2.18     0.76      1578       704  1.00 
#> r_re[2]          -0.98   0.66    -0.97    -2.36     0.30      1593       894  1.01 
#> r_re[3]           0.02   0.77    -0.12    -1.50     1.54      2322       746  1.00 
#> r_re[4]           0.60   0.77     0.63    -0.91     2.10      1945       756  1.01 
#> r_re[5]          -0.39   0.75    -0.67    -1.78     1.05      2142       969  1.00 
  # }

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
#> 
#> Warning: Best optimization run ended with singular convergence. Estimates may be usable, but check opt_history or try more starts.
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
#> b[skill]   0.19129     0.06369    0.06646    0.31613 
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
#> 
#> SE warning: sdreport() produced non-finite standard errors; Hessian-based fallback will be attempted.
#> SE warning: Hessian matrix was singular; using MASS::ginv() to approximate the covariance matrix.
  map_ssp$summary("b")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 390.80
#> Approx. Log Marginal Likelihood (Laplace): NA
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[talk]    0.26100     0.05297    0.15717    0.36483 
#> b[perf]    0.15004     0.02966    0.09192    0.20817 
#> b[skill]   0.18076     0.06487    0.05362    0.30790 
#> 
```
