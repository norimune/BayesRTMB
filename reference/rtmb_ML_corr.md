# Multilevel Correlation Analysis (Kenny's Model)

This function fits a multilevel correlation model to estimate
between-group and within-group covariance/correlation structures, along
with Intraclass Correlation Coefficients (ICC).

## Usage

``` r
rtmb_ML_corr(
  formula,
  ID,
  data = NULL,
  prior = prior_uniform(),
  y_range = NULL,
  ...
)
```

## Arguments

- formula:

  A formula (e.g., `cbind(V1, V2) ~ 1`) or a matrix of response
  variables.

- ID:

  A character string or expression specifying the group ID variable.

- data:

  A data frame.

- prior:

  Prior configuration object:
  [`prior_uniform()`](https://norimune.github.io/BayesRTMB/reference/prior_uniform.md)
  or
  [`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md).
  Default is
  [`prior_uniform()`](https://norimune.github.io/BayesRTMB/reference/prior_uniform.md).

- y_range:

  Optional numeric vector or matrix defining the theoretical range (min,
  max) of response variables. Required when using
  [`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md).
  Can be a vector of length 2 (applies to all variables) or a
  matrix/list of length P.

- ...:

  Additional arguments passed to `rtmb_model`.

## Value

A `RTMB_Model` object.
