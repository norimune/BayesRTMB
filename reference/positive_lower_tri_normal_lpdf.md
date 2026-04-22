# \#' Positive lower-triangular normal log-probability density function

\#' Positive lower-triangular normal log-probability density function

## Usage

``` r
positive_lower_tri_normal_lpdf(x, mean = 0, sd = 1)
```

## Arguments

- x:

  A matrix of lower-triangular parameters (Cholesky factor with positive
  diagonals).

- mean:

  Mean of the normal distribution (assumed to be 0 for the half-normal
  correction).

- sd:

  Standard deviation of the normal distribution.

## Value

The log-density calculated for the non-zero lower triangular elements.
