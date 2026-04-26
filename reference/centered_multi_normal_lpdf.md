# Centered / Centered matrix multivariate normal log-probability density function

Centered / Centered matrix multivariate normal log-probability density
function

## Usage

``` r
centered_multi_normal_lpdf(x, sigma, K = NULL)
```

## Arguments

- x:

  Vector or matrix of quantiles. If a matrix, it assumes each column
  sums to zero.

- sigma:

  Standard deviation parameter. Can be a scalar, or a vector of length
  ncol(x) if x is a matrix.

- K:

  Dimension. If NULL, inferred from the length or nrow of x.

## Value

The log-density.
