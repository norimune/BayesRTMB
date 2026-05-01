# Normal log-probability density function

Normal log-probability density function

## Usage

``` r
normal_lpdf(x, mean, sd, sum = TRUE)
```

## Arguments

- x:

  Vector of quantiles.

- mean:

  Vector of means.

- sd:

  Vector of standard deviations.

- sum:

  Logical; if TRUE (default), returns the sum of log-densities. If
  FALSE, returns element-wise log-densities.

## Value

The sum of the log-density (scalar) or a vector of log-densities.
