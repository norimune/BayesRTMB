# Wishart log-probability density function parameterized by Cholesky factor of correlation matrix

Wishart log-probability density function parameterized by Cholesky
factor of correlation matrix

## Usage

``` r
wishart_CF_lpdf(X, n, sd, CF_Omega)
```

## Arguments

- X:

  Symmetric positive-definite matrix (scatter matrix).

- n:

  Degrees of freedom (sample size N).

- sd:

  Vector of standard deviations.

- CF_Omega:

  Cholesky factor of the correlation matrix.

## Value

The log-density with all constant terms.
