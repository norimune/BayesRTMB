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
  prior = prior_uniform(),
  fixed = NULL,
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

  Prior configuration object:
  [`prior_uniform()`](https://norimune.github.io/BayesRTMB/reference/prior_uniform.md)
  (default),
  [`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md),
  [`prior_rhs()`](https://norimune.github.io/BayesRTMB/reference/prior_rhs.md),
  or
  [`prior_ssp()`](https://norimune.github.io/BayesRTMB/reference/prior_ssp.md).

- ...:

  Additional arguments passed to `rtmb_model`.

## Value

A `RTMB_Model` object.
