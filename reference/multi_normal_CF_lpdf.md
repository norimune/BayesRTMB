# Multivariate normal log-probability density function parameterized by Cholesky factor of correlation matrix

Multivariate normal log-probability density function parameterized by
Cholesky factor of correlation matrix

## Usage

``` r
multi_normal_CF_lpdf(x, mean, sd, CF_Omega)
```

## Arguments

- x:

  Vector or matrix of quantiles.

- mean:

  Vector or matrix of means.

- sd:

  Vector of standard deviations.

- CF_Omega:

  Cholesky factor of the correlation matrix.

## Value

The sum of the log-density.
