# Fit a Latent Rank Theory (LRT) Model

Fits a Latent Rank Theory model, which is a mixture model with ordered
ranks and Gaussian Process smoothing on the mean profiles.

## Usage

``` r
rtmb_lrt(
  formula,
  k = 3,
  data = NULL,
  rank_coords = NULL,
  covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
  magnitude = NULL,
  smoothing = NULL,
  noise = 0.01,
  prob_smoothing = FALSE,
  link = c("ordered", "sequential"),
  prior = prior_uniform(),
  ...
)
```

## Arguments

- formula:

  A formula specifying the response variable(s).

- k:

  Number of ranks (mixture components).

- data:

  A data frame containing the variables.

- magnitude:

  Signal standard deviation for the GP prior. If NULL, it is estimated.

- smoothing:

  Length-scale for the GP prior. If NULL, it is estimated.

- noise:

  Measurement noise for the GP prior (default is 0.01).

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
