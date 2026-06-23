# Gaussian Process Log-Density (Squared Exponential Kernel)

Calculates the log-density of a Gaussian Process with a Squared
Exponential (RBF) kernel.

## Usage

``` r
gaussian_process_lpdf(
  y,
  x,
  mean = 0,
  magnitude = 1,
  smoothing = 1,
  noise = 0.01,
  sum = TRUE
)
```

## Arguments

- y:

  Observation vector (N), or matrix (M x N) whose rows are independent
  GP realizations on the same coordinates.

- x:

  Coordinate vector or matrix (N x D), where N matches `length(y)` for
  vector `y` or `ncol(y)` for matrix `y`.

- mean:

  Mean vector (scalar, length N, or M x N matrix for matrix `y`).

- magnitude:

  Signal standard deviation (alpha).

- smoothing:

  Length-scale (rho).

- noise:

  Measurement noise standard deviation (sigma).

- sum:

  Logical; whether to return the sum of log-densities. If `FALSE` and
  `y` is a matrix, returns one log-density per row.

## Value

Log-density value.
