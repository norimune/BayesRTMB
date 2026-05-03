# Fit a Correlation Model using RTMB'

Fits a correlation model to estimate means, standard deviations, and
correlation structures. Supports both bivariate and multivariate data.

## Usage

``` r
rtmb_corr(
  data,
  prior = prior_uniform(),
  y_range = NULL,
  init = NULL,
  null = NULL,
  ...
)
```

## Arguments

- data:

  A data frame or matrix of response variables.

- prior:

  Prior configuration object:
  [`prior_uniform()`](https://norimune.github.io/BayesRTMB/reference/prior_uniform.md)
  (default) or
  [`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md).

- y_range:

  Optional numeric vector or matrix defining the theoretical range (min,
  max) of response variables. Required when using
  [`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md).

- init:

  Optional list of initial values.

- null:

  Optional list specifying parameters to fix to null values.

- ...:

  Additional arguments passed to `rtmb_model`.

## Value

A `RTMB_Model` object.
