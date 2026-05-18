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

- ...:

  Additional arguments passed to \`rtmb_model\`.

## Value

A `RTMB_Model` object.
