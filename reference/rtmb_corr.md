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
  method = c("pearson", "spearman", "reml"),
  prior = prior_flat(),
  y_range = NULL,
  init = NULL,
  fixed = NULL
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

- covariates:

  Optional numeric matrix or data frame of covariates to be included in
  the joint MVN model.

- method:

  Correlation method for `classic()`: `"pearson"`, `"spearman"`, or
  `"reml"`.

- prior:

  Prior configuration object: \`prior_flat()\`, \`prior_normal()\`, or
  \`prior_weak()\`. Default is \`prior_flat()\`.

- y_range:

  Optional numeric vector or matrix defining the theoretical range (min,
  max) of response variables. Required when using
  [`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md).
  Can be a vector of length 2 (applies to all variables) or a
  matrix/list of length P.

- init:

  Optional list of initial values.

- fixed:

  Optional named list of fixed values for specific parameters.

## Value

A `RTMB_Model` object.
