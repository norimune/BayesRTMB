# Fit a Latent Rank Theory (LRT) Model

Fits a Latent Rank Theory model, which is a mixture model with ordered
ranks and Gaussian Process smoothing on the mean profiles.

## Usage

``` r
rtmb_lrt(
  formula,
  k = 3,
  data = NULL,
  covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
  magnitude = NULL,
  smoothing = NULL,
  noise = 0.01,
  prob_smoothing = FALSE,
  link = c("ordered", "sequential"),
  prior_type = c("weakly_informative", "uniform"),
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

- prior_type:

  Prior type: "weakly_informative" or "uniform".

- ...:

  Additional arguments passed to `rtmb_model`.

- prob_formula:

  Optional formula for latent class regression.

## Value

A `RTMB_Model` object.
