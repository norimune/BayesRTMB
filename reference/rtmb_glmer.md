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

- .force_sum:

  Logical; internal use only.

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
#> Intercept  39.68627     1.63262   36.48640   42.88614 
#> b[wt]      -3.19097     0.72055   -4.60323   -1.77871 
#> b[cyl]     -1.50780     0.39477   -2.28153   -0.73406 
#> sigma       2.44420     0.30553    1.91310    3.12275 
#> 

  # --- 2. Generalized Linear Model (rtmb_glm) ---
  # Fit a logistic regression model
  fit_glm <- rtmb_glm(am ~ mpg + hp, data = mtcars, family = "binomial")
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
#> Negative Log-Posterior: 9.62
#> Approx. Log Marginal Likelihood (Laplace): -12.65
#> 
#> Point Estimates and 95% Wald CI:
#>  variable   Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept  -33.60517    15.12417  -63.24799   -3.96235 
#> b[mpg]       1.25961     0.56938    0.14365    2.37558 
#> b[hp]        0.05504     0.02698    0.00217    0.10792 
#> 

  # --- 3. Generalized Linear Mixed Model (rtmb_glmer) ---
  # Fit a linear mixed-effects model with random intercepts
  fit_glmer <- rtmb_glmer(weight ~ Time + (1 | Chick), data = ChickWeight, family = "gaussian")
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
#> Negative Log-Posterior: 2811.17
#> Approx. Log Marginal Likelihood (Laplace): -2806.94
#> Note: Random effects are stored in $random_effects (use ranef = TRUE to show them)
#> 
#> Point Estimates and 95% Wald CI:
#>      variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept      27.84415     4.35087   19.31661   36.37170 
#> b[Time]         8.72626     0.17536    8.38256    9.06996 
#> sigma          28.24713     0.86881   26.59460   30.00234 
#> sd[Chick:Int]  26.49976     2.91129   21.36626   32.86665 
#> 

  # MCMC sampling (chains and iterations reduced for faster execution)
  # \donttest{
  mcmc_glmer <- fit_glmer$sample(sampling = 500, warmup = 500, chains = 2)
#> Starting sequential sampling (chains = 2)...
#> chain 1 started... 
#> chain 1: iter 100 warmup 
#> chain 1: iter 200 warmup 
#> chain 1: iter 300 warmup 
#> chain 1: iter 400 warmup 
#> chain 1: iter 500 warmup 
#> chain 1: iter 600 sampling 
#> chain 1: iter 700 sampling 
#> chain 1: iter 800 sampling 
#> chain 1: iter 900 sampling 
#> chain 1: iter 1000 sampling 
#> chain 2 started... 
#> chain 2: iter 100 warmup 
#> chain 2: iter 200 warmup 
#> chain 2: iter 300 warmup 
#> chain 2: iter 400 warmup 
#> chain 2: iter 500 warmup 
#> chain 2: iter 600 sampling 
#> chain 2: iter 700 sampling 
#> chain 2: iter 800 sampling 
#> chain 2: iter 900 sampling 
#> chain 2: iter 1000 sampling 
  mcmc_glmer$summary()
#>      variable      mean    sd       map      q2.5     q97.5  ess_bulk  ess_tail  rhat 
#> lp             -2817.40  7.05  -2817.51  -2832.93  -2804.62       170       307  1.01 
#> Intercept         27.26  4.50     27.89     18.39     35.50        70       245  1.06 
#> b[Time]            8.71  0.17      8.67      8.39      9.05       361       504  1.00 
#> sigma             28.34  0.88     28.29     26.57     30.00       415       491  1.00 
#> sd[Chick:Int]     27.17  3.00     26.44     21.83     33.60       116       252  1.03 
#> r_re[1]            0.03  0.63      0.10     -1.21      1.18       353       471  1.01 
#> r_re[2]           -0.95  0.38     -0.95     -1.68     -0.18       320       544  1.01 
#> r_re[3]           -0.92  0.39     -0.75     -1.71     -0.21       337       625  1.01 
#> r_re[4]           -1.85  0.37     -1.82     -2.56     -1.06       291       260  1.02 
#> r_re[5]           -1.41  0.35     -1.30     -2.12     -0.74       334       370  1.00 
  # }

  # --- 4. Linear Mixed Model (rtmb_lmer) ---
  # A convenient wrapper for Gaussian mixed models (identical to rtmb_glmer with family="gaussian")
  fit_lmer <- rtmb_lmer(mpg ~ wt + (1 | cyl), data = mtcars)
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
#> Negative Log-Posterior: 78.24
#> Approx. Log Marginal Likelihood (Laplace): -75.65
#> Note: Random effects are stored in $random_effects (use ranef = TRUE to show them)
#> 
#> Point Estimates and 95% Wald CI:
#>    variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept    32.12009     2.87856   26.47821   37.76196 
#> b[wt]        -3.70755     0.81856   -5.31191   -2.10320 
#> sigma         2.52988     0.33760    1.94765    3.28616 
#> sd[cyl:Int]   2.09732     1.13099    0.72888    6.03498 
#> 

  # --- 5. Regularized Regression (Variable Selection) ---
  # You can apply regularization to the fixed effects to shrink noise variables towards zero.
  # Use prior = prior_rhs() for the Regularized Horseshoe prior, or prior_ssp() for the Spike-and-Slab prior.
  # Note: When using regularization, you must specify 'y_range' (the theoretical minimum and maximum
  # values of the response variable) to automatically set up the required weakly informative priors.

  # Fit a linear regression using all predictors in mtcars with the Horseshoe prior
  # 'mpg' theoretically ranges roughly between 0 and 40
  fit_rhs <- rtmb_lm(mpg ~ ., data = mtcars, prior = prior_rhs(), y_range = c(0, 40))
#> Error: Variables not found in data: .
  map_rhs <- fit_rhs$optimize()
#> Error: object 'fit_rhs' not found
  # Summarize only the fixed effects (slopes)
  map_rhs$summary("b")
#> Error: object 'map_rhs' not found

  # Fit a linear regression with the Spike-and-Slab prior
  fit_ssp <- rtmb_lm(mpg ~ ., data = mtcars, prior = prior_ssp(), y_range = c(0, 40))
#> Error: Variables not found in data: .
  map_ssp <- fit_ssp$optimize()
#> Error: object 'fit_ssp' not found
  map_ssp$summary("b")
#> Error: object 'map_ssp' not found

  # For models with complex penalties, MCMC often provides more reliable credible intervals
  # \donttest{
  mcmc_ssp <- fit_ssp$sample(sampling = 500, warmup = 500, chains = 2)
#> Error: object 'fit_ssp' not found
  mcmc_ssp$summary("b")
#> Error: object 'mcmc_ssp' not found
  mcmc_ssp$summary("b")
#> Error: object 'mcmc_ssp' not found
  # }
```
