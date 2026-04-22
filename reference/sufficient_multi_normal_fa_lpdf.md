# Sufficient statistics factor analysis multivariate normal log-probability density function

Sufficient statistics factor analysis multivariate normal
log-probability density function

## Usage

``` r
sufficient_multi_normal_fa_lpdf(S_mat, N, y_bar, mu, psi, Lambda)
```

## Arguments

- S_mat:

  Deviation sum of squares matrix.

- N:

  Sample size.

- y_bar:

  Sample mean vector.

- mu:

  Mean parameter vector.

- psi:

  Vector of unique variances (P).

- Lambda:

  Factor loading matrix (P x K).

## Value

The exact log-likelihood of the N raw observations.
