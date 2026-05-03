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
  prior_type = c("uniform", "weakly_informative"),
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

- prior_type:

  Prior type: "uniform" (default) or "weakly_informative".

- ...:

  Additional arguments passed to `rtmb_model`.

- prob_formula:

  Optional formula for latent class regression (covariates for
  probabilities).

## Value

A `RTMB_Model` object.
