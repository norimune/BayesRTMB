# Mixture Model Wrapper for RTMB

Provides a user-friendly interface for fitting Gaussian mixture models
with optional covariates on class membership probabilities and various
covariance structures.

## Usage

``` r
rtmb_mixture(
  formula,
  k = 2,
  data = NULL,
  covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
  prior = prior_flat(),
  y_range = NULL,
  fixed = NULL,
  WAIC = FALSE,
  ...
)
```

## Arguments

- formula:

  A formula specifying the response variable(s). For multivariate, use
  `cbind(y1, y2) ~ 1`.

- k:

  Number of mixture components.

- data:

  A data frame containing the variables in the model.

- covariance:

  Covariance structure: "diagonal" (default), "diagonal_equal", "full",
  "full_equal", or "full_equal_corr".

- prior:

  Prior configuration: \`prior_flat()\`, \`prior_normal()\`,
  \`prior_weak()\`, \`prior_rhs()\`, or \`prior_ssp()\`. Default is
  \`prior_flat()\`. If \`y_range\` is supplied with the default flat
  prior, the wrapper automatically switches to \`prior_weak()\`.

- y_range:

  Optional numeric vector or matrix defining the theoretical range (min,
  max) of response variables. Specifying this automatically enables
  weakly informative priors if \`prior\` is \`prior_flat()\`.

- fixed:

  Optional named list of fixed values for specific parameters.

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

- ...:

  Additional arguments passed to \`rtmb_model\`.

## Value

A `RTMB_Model` object.

## Examples

``` r

  # Simulate 1D mixture data (2 components)
  set.seed(123)
  N <- 100
  group <- rbinom(N, 1, 0.4)
  y <- ifelse(group == 1, rnorm(N, 5, 1), rnorm(N, 0, 1))
  dat <- data.frame(y)

  # Fit a 1D Gaussian mixture model with 2 components
  fit_mix <- rtmb_mixture(y ~ 1, k = 2, data = dat)
#> Pre-checking model code...
#> Checking RTMB setup...
  
  # MAP estimation
  map_mix <- fit_mix$optimize()
#> Starting RTMB optimization...
#> 
  map_mix$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 200.57
#> Approx. Log Marginal Likelihood (Laplace): -207.18
#> 
#> Point Estimates and 95% Wald CI:
#>     variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> prob_mean[1]   0.41451     0.05103    0.31450    0.51452 
#> prob_mean[2]   0.58549     0.05103    0.48548    0.68550 
#> mu[C1]         4.81787     0.17776    4.46946    5.16628 
#> mu[C2]         0.04764     0.12571   -0.19874    0.29402 
#> sigma[C1]      1.00445     0.14877    0.75137    1.34278 
#> sigma[C2]      0.88706     0.09688    0.71612    1.09880 
#> theta[1]       0.41451     0.05103    0.31920    0.51668 
#> theta[2]       0.58549     0.05103    0.48332    0.68080 
#> 
```
