# Diffusion model log-probability density function

Log-likelihood for a two-boundary drift diffusion model. The first
argument is the response time and the second argument is the binary
response. Positive response values are treated as upper-boundary
responses; zero or negative values are treated as lower-boundary
responses.

If `x` is a matrix, rows are treated as independent units and columns as
repeated trials. In that case, `sum = FALSE` returns one log-likelihood
contribution per row. If `x` is a vector, `sum = FALSE` returns one
contribution per observation.

## Usage

``` r
diffusion_lpdf(x, response, alpha, tau, beta, delta, K_diff = 20, sum = TRUE)
```

## Arguments

- x:

  Response times. A vector, or a matrix with units in rows and trials in
  columns.

- response:

  Binary responses with the same shape as `x`. Positive values indicate
  the upper boundary.

- alpha:

  Boundary separation. A scalar, row-wise vector, trial-wise vector, or
  matrix compatible with `x`.

- tau:

  Non-decision time. A scalar, row-wise vector, trial-wise vector, or
  matrix compatible with `x`.

- beta:

  Initial bias as a proportion of the boundary separation. A scalar,
  row-wise vector, trial-wise vector, or matrix compatible with `x`.

- delta:

  Drift rate. A scalar, row-wise vector, trial-wise vector, or matrix
  compatible with `x`.

- K_diff:

  Number of terms in the truncated infinite-series approximation.

- sum:

  Logical; if `TRUE` (default), returns the summed log-likelihood.

## Value

The summed log-likelihood, row-wise log-likelihoods for matrix input, or
observation-wise log-likelihoods for vector input.
