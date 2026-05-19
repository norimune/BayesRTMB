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
  fixed = NULL,
  missing = c("listwise", "fiml", "pairwise")
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

- missing:

  Missing value handling strategy: "listwise" (default), "pairwise", or
  "fiml" (Full Information Maximum Likelihood).

## Examples

``` r

  # Estimate the correlation between two variables in the debate dataset
  data(debate, package = "BayesRTMB")

  fit_corr <- rtmb_corr(cbind(sat, perf), data = debate)
#> Pre-checking model code...
#> Checking RTMB setup...

  if (FALSE) { # \dontrun{
  mcmc_corr <- fit_corr$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_corr$summary()

  bf_corr <- mcmc_corr$bayes_factor(fixed = list(corr = 0))
  print(bf_corr)
  } # }
```
