# Fit a Correlation Model using RTMB

\`rtmb_corr\` fits a correlation model to estimate means, standard
deviations, and correlation structures. It supports simple correlation,
multilevel correlation, and classical frequentist estimation.

## Usage

``` r
rtmb_corr(
  x = NULL,
  data = NULL,
  ID = NULL,
  covariates = NULL,
  prior = prior_uniform(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  null = NULL,
  ...
)
```

## Arguments

- x:

  A matrix, data frame, formula, or expression (e.g., `cbind(V1, V2)`)
  of response variables.

- data:

  An optional data frame containing the variables.

- ID:

  A character string or expression specifying the group ID variable for
  multilevel models.

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

- init:

  Optional list of initial values.

- null:

  Optional list specifying parameters to fix to null values.

- ...:

  Additional arguments passed to `rtmb_model`.

## Value

A `RTMB_Model` object.
