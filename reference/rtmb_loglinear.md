# RTMB-based Log-linear analysis (Poisson regression)

Performs Bayesian or Frequentist log-linear analysis (Poisson
regression) on a contingency table or raw data.

## Usage

``` r
rtmb_loglinear(
  formula,
  data,
  prior = prior_flat(),
  y_range = NULL,
  fixed = NULL,
  WAIC = FALSE,
  ...
)
```

## Arguments

- formula:

  A formula (e.g., \`~ A + B + A:B\`) or a contingency table.

- data:

  A data frame (required if \`formula\` is used).

- prior:

  An object of class "rtmb_prior" specifying the prior distribution.

- y_range:

  Optional theoretical minimum and maximum values of the count response.
  If supplied with the default flat prior, \`prior_weak()\` is used.

- fixed:

  Optional named list of fixed values for specific parameters.

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

- ...:

  Additional arguments passed to \`rtmb_glm()\`.

## Value

An \`RTMB_Model\` object.

## Examples

``` r

  # Create a contingency table
  tab <- matrix(c(10, 20, 30, 40), nrow = 2)
  dimnames(tab) <- list(A = c("A1", "A2"), B = c("B1", "B2"))
  
  # Fit a log-linear model (independence model: ~ A + B)
  fit_log <- rtmb_loglinear(~ A + B, data = tab)
#> Pre-checking model code...
#> Checking RTMB setup...
  
  # MAP estimation
  map_log <- fit_log$optimize()
#> Starting RTMB optimization...
  map_log$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 10.29
#> Approx. Log Marginal Likelihood (Laplace): -12.95
#> 
#> Point Estimates and 95% Wald CI:
#>  variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept   2.48491     0.21985    2.05401    2.91580 
#> b[AA2]      0.40547     0.20412    0.00539    0.80554 
#> b[BB2]      0.84730     0.21822    0.41960    1.27500 
#> mu[1,1]    12.00000     2.63819    6.82923   17.17077 
#> mu[2,1]    18.00000     3.60002   10.94409   25.05591 
#> mu[3,1]    28.00000     4.79168   18.60847   37.39153 
#> mu[4,1]    42.00000     6.07950   30.08439   53.91561 
#> 
```
