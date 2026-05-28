# RTMB-based GLMM wrapper function

RTMB-based GLMM wrapper function

## Usage

``` r
rtmb_glmer(
  formula,
  data,
  family = "gaussian",
  laplace = FALSE,
  prior = prior_flat(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  gmc = NULL,
  centering = NULL,
  cwc = NULL,
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

  Data frame

- family:

  Character string of the distribution family (e.g., "gaussian",
  "binomial", "poisson", "ordered", "sequential")

- laplace:

  Logical; whether to marginalize random effects using Laplace
  approximation

- prior:

  An object of class "rtmb_prior" specifying the prior distribution. Use
  prior_weak(), prior_rhs(), or prior_ssp(). Default is
  \`prior_flat()\`.

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
  (group variable) and `pars` (variable names to center).

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
#> chain 1: iter 200 warmup 
#> chain 1: iter 400 warmup 
#> chain 1: iter 600 sampling 
#> chain 1: iter 800 sampling 
#> chain 1: iter 1000 sampling 
#> chain 2 started... 
#> chain 2: iter 200 warmup 
#> chain 2: iter 400 warmup 
#> chain 2: iter 600 sampling 
#> chain 2: iter 800 sampling 
#> chain 2: iter 1000 sampling 
  mcmc_glmer$summary()
#>      variable     mean     sd      map     q2.5    q97.5  ess_bulk  ess_tail  rhat 
#> lp             -512.04  13.30  -512.86  -541.01  -485.69        73       121  1.05 
#> Intercept         2.64   0.09     2.65     2.47     2.82       537       524  1.00 
#> b[cond]           0.76   0.12     0.77     0.51     0.99       572       356  1.00 
#> sigma             0.83   0.05     0.82     0.75     0.92       244       647  1.02 
#> sd[group:Int]     0.40   0.10     0.42     0.16     0.56       146        67  1.03 
#> r_re[1]          -0.71   0.78    -0.74    -2.22     0.97      1105       574  1.00 
#> r_re[2]          -0.96   0.84    -0.81    -2.58     0.64       959       553  1.00 
#> r_re[3]           0.02   0.79    -0.17    -1.62     1.72      2007       521  1.00 
#> r_re[4]           0.58   0.76     0.39    -0.94     2.01      1566       743  1.01 
#> r_re[5]          -0.40   0.78    -0.25    -1.95     1.15      1741       562  1.00 
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
#> Warning: Optimization did not converge ( convergence code = 1; message = singular convergence (7)). Estimates may be unreliable; consider increasing num_estimate, changing initial values, or adjusting optimizer control settings.
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
#> 
#> Warning: Optimization did not converge ( convergence code = 1; message = singular convergence (7)). Estimates may be unreliable; consider increasing num_estimate, changing initial values, or adjusting optimizer control settings.
#> SE warning: sdreport() produced non-finite standard errors; Hessian-based fallback will be attempted.
#> SE warning: Hessian matrix was singular; using MASS::ginv() to approximate the covariance matrix.
  map_ssp$summary("b")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 400.78
#> Approx. Log Marginal Likelihood (Laplace): NA
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[talk]    0.27952     0.05518    0.17136    0.38768 
#> b[perf]    0.00000     0.00000    0.00000    0.00000 
#> b[skill]   0.18934     0.06775    0.05656    0.32211 
#> 

  # For models with complex penalties, MCMC often provides more reliable credible intervals
  # \donttest{
  mcmc_ssp <- fit_ssp$sample(sampling = 500, warmup = 500, chains = 2)
#> Starting sequential sampling (chains = 2)...
#> chain 1 started... 
#> chain 1: iter 200 warmup 
#> chain 1: iter 400 warmup 
#> chain 1: iter 600 sampling 
#> chain 1: iter 800 sampling 
#> chain 1: iter 1000 sampling 
#> chain 2 started... 
#> chain 2: iter 200 warmup 
#> chain 2: iter 400 warmup 
#> chain 2: iter 600 sampling 
#> chain 2: iter 800 sampling 
#> chain 2: iter 1000 sampling 
#> Calculating transformed parameters...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
  mcmc_ssp$summary("b")
#> variable  mean    sd   map   q2.5  q97.5  ess_bulk  ess_tail  rhat 
#> b[talk]   0.26  0.05  0.27   0.15   0.36      1022       707  1.00 
#> b[perf]   0.15  0.03  0.15   0.09   0.21      1080       628  1.00 
#> b[skill]  0.15  0.08  0.17  -0.00   0.29       205       293  1.01 
  mcmc_ssp$summary("b")
#> variable  mean    sd   map   q2.5  q97.5  ess_bulk  ess_tail  rhat 
#> b[talk]   0.26  0.05  0.27   0.15   0.36      1022       707  1.00 
#> b[perf]   0.15  0.03  0.15   0.09   0.21      1080       628  1.00 
#> b[skill]  0.15  0.08  0.17  -0.00   0.29       205       293  1.01 
  # }
```
