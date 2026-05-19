# Multivariate normal log-probability density function parameterized by Cholesky factor of correlation matrix (FIML)

Multivariate normal log-probability density function parameterized by
Cholesky factor of correlation matrix (FIML)

## Usage

``` r
fiml_multi_normal_CF_lpdf(x, mean, sd, CF_Omega)
```

## Arguments

- x:

  Vector or matrix of quantiles containing NAs.

- mean:

  Vector of means.

- sd:

  Vector of standard deviations.

- CF_Omega:

  Cholesky factor of the correlation matrix.

## Value

The sum of the log-density using Full Information Maximum Likelihood.
