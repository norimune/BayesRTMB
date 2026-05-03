# Fit a Mixture Model

Fits a mixture model using the RTMB engine. Supports a formula interface
(response ~ covariates) for latent class regression.

## Usage

``` r
rtmb_mixture(
  formula,
  k,
  data = NULL,
  covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
  prior = prior_uniform(),
  multivariate = NULL,
  ...
)
```

## Arguments

- formula:

  A formula (e.g., `y ~ x1 + x2`) for latent class regression, or a data
  object (matrix/vector) for a basic mixture.

- k:

  Integer; the number of mixture components (clusters).

- data:

  A data frame (required if formula is used).

- prior:

  An `rtmb_prior` object. Default is
  [`prior_uniform()`](https://norimune.github.io/BayesRTMB/reference/prior_uniform.md).

- multivariate:

  Logical; whether to treat the data as multivariate.

- ...:

  Additional arguments passed to `rtmb_model`.

## Value

An `RTMB_Model` object.
