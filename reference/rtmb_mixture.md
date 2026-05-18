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

## Examples

``` r
if (FALSE) { # \dontrun{
  # Simulate 1D mixture data (2 components)
  set.seed(123)
  N <- 100
  group <- rbinom(N, 1, 0.4)
  y <- ifelse(group == 1, rnorm(N, 5, 1), rnorm(N, 0, 1))
  dat <- data.frame(y)

  # Fit a 1D Gaussian mixture model with 2 components
  fit_mix <- rtmb_mixture(y ~ 1, k = 2, data = dat)
  
  # MAP estimation
  map_mix <- fit_mix$optimize()
  map_mix$summary()

  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_mix <- fit_mix$sample(sampling = 500, warmup = 500, chains = 2)
  # MCMC summary provides estimates for component means, standard deviations, and mixture probabilities
  mcmc_mix$summary()
} # }
```
