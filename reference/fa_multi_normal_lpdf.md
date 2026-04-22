# Factor analysis multivariate normal log-probability density function

Woodbury matrix identity is used for efficient computation.

## Usage

``` r
fa_multi_normal_lpdf(x, mu, Lambda, psi)
```

## Arguments

- x:

  Vector or matrix of quantiles.

- mu:

  Vector of means.

- Lambda:

  Factor loading matrix (P x K).

- psi:

  Vector of unique variances (P).

## Value

The sum of the log-density.
