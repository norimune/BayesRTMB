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
  prior = prior_flat(),
  y_range = NULL,
  fixed = NULL,
  two_stage = FALSE,
  WAIC = FALSE,
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

- rank_coords:

  Optional numeric vector of coordinates for each rank. Default is 1:k.

- covariance:

  Covariance structure: "diagonal", "diagonal_equal", "full",
  "full_equal", or "full_equal_corr".

- magnitude:

  Signal standard deviation for the GP prior. If NULL, it is estimated.

- smoothing:

  Length-scale for the GP prior. If NULL, it is estimated.

- noise:

  Measurement noise for the GP prior (default is 0.01).

- prob_smoothing:

  Logical; whether to apply smoothing to the class membership
  probabilities.

- link:

  Link function for class probabilities: "ordered" or "sequential".

- prior:

  Prior configuration: \`prior_flat()\`, \`prior_normal()\`,
  \`prior_weak()\`, \`prior_rhs()\`, or \`prior_ssp()\`. Default is
  \`prior_flat()\`. If \`y_range\` is supplied with the default flat
  prior, the wrapper automatically switches to \`prior_weak()\`.

- y_range:

  Optional theoretical range of the response variables. Use a numeric
  vector \`c(min, max)\` for a common range, or a P x 2 matrix/list for
  response-specific ranges. Specifying this automatically enables
  \`prior_weak()\` when \`prior\` is \`prior_flat()\` and calibrates the
  residual standard-deviation prior.

- fixed:

  Optional named list of fixed values for specific parameters.

- two_stage:

  Logical; if TRUE, estimate the latent-rank measurement model first and
  then estimate the rank regression with delta-method uncertainty
  propagation. Currently supported for \`\$optimize()\` only.

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

- ...:

  Additional arguments passed to \`rtmb_model\`.

## Value

A `RTMB_Model` object.
