# RTMB-based Mediation Analysis Wrapper

\`rtmb_mediation\` performs mediation analysis by simultaneously
estimating multiple GLM regression equations. It automatically
identifies mediation paths and calculates indirect, direct, and total
effects.

## Usage

``` r
rtmb_mediation(
  formula,
  data,
  family = "gaussian",
  prior = prior_flat(),
  y_range = NULL,
  fixed = NULL,
  view = NULL,
  WAIC = FALSE,
  ...
)
```

## Arguments

- formula:

  A list of formulas defining the regression paths (e.g., \`list(M ~ X,
  Y ~ X + M)\`).

- data:

  A data frame containing the variables.

- family:

  A single character string or a list of character strings specifying
  the error distribution for each equation (e.g., \`family =
  list("gaussian", "binomial")\`). Default is "gaussian".

- prior:

  An object of class "rtmb_prior" specifying the prior distribution.

- y_range:

  Theoretical minimum and maximum values of the response variable.

- fixed:

  A named list of parameter values to fix (optional).

- view:

  Character vector of parameter names to prioritize in summary.

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

- ...:

  Reserved; unused arguments are rejected.

## Value

An \`RTMB_Model\` object.

## Details

The function identifies mediation paths by looking for variables that
are responses in one equation and predictors in another. Indirect
effects are calculated as the product of coefficients along these paths
(\\a \* b\\).

**Uncertainty Estimation**: When using \`\$optimize(ci_method =
"sampling")\`, the function provides asymmetric confidence intervals for
indirect effects based on the distribution of products, which is more
accurate than the standard Sobel test (Delta Method).

## Examples

``` r

  # Simulate mediation data
  set.seed(123)
  N <- 100
  X <- rnorm(N)
  # M is influenced by X
  M <- 0.5 * X + rnorm(N, 0, 0.5)
  # Y is influenced by both X and M
  Y <- 0.3 * X + 0.8 * M + rnorm(N, 0, 0.5)
  dat <- data.frame(X, M, Y)

  # Fit a mediation model
  # The formula list specifies the two regression equations
  fit_med <- rtmb_mediation(list(M ~ X, Y ~ X + M), data = dat)
#> Pre-checking model code...
#> Checking RTMB setup...
  
  # Maximum A Posteriori (MAP) estimation
  map_med <- fit_med$optimize()
#> Starting RTMB optimization...
#> 
  # The summary automatically calculates indirect, direct, and total effects
  map_med$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 134.66
#> Approx. Log Marginal Likelihood (Laplace): -149.35
#> 
#> Point Estimates and 95% Wald CI:
#>      variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b1[Intercept]  -0.05140     0.04829   -0.14604    0.04324 
#> b1[X]           0.47376     0.05290    0.37008    0.57745 
#> b2[Intercept]   0.06753     0.04734   -0.02526    0.16032 
#> b2[X]           0.22151     0.06924    0.08580    0.35721 
#> b2[M]           0.82381     0.09750    0.63272    1.01490 
#> IE_X_M_Y        0.39029     0.06351    0.26582    0.51476 
#> DE_X_Y          0.22151     0.06924    0.08580    0.35721 
#> TE_X_M_Y        0.61180     0.06753    0.47945    0.74415 
#> sigma1          0.48048     0.03398    0.41830    0.55191 
#> sigma2          0.46846     0.03313    0.40783    0.53810 
#> 
```
