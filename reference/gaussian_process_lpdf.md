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

  Observation vector (N).

- x:

  Coordinate vector or matrix (N x D).

- mean:

  Mean vector (scalar or length N).

- magnitude:

  Signal standard deviation (alpha).

- smoothing:

  Length-scale (rho).

- noise:

  Measurement noise standard deviation (sigma).

- sum:

  Logical; whether to return the sum of log-densities.

## Value

Log-density value.
