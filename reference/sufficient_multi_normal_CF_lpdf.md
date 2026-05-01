# Sufficient statistics multivariate normal log-probability density function

Sufficient statistics multivariate normal log-probability density
function

## Usage

``` r
sufficient_multi_normal_CF_lpdf(S_mat, N, y_bar, mean, sd, CF_Omega)
```

## Arguments

- S_mat:

  Deviation sum of squares matrix.

- N:

  Sample size.

- y_bar:

  Sample mean vector.

- mean:

  Mean parameter vector.

- sd:

  Standard deviation parameter vector.

- CF_Omega:

  Cholesky factor of correlation matrix.

## Value

The exact log-likelihood of the N raw observations.
