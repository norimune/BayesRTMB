# Factor analysis multivariate normal log-probability density function (FIML)

Factor analysis multivariate normal log-probability density function
(FIML)

## Usage

``` r
fiml_multi_normal_fa_lpdf(x, mu, Lambda, psi)
```

## Arguments

- x:

  Matrix of quantiles containing NAs.

- mu:

  Vector of means.

- Lambda:

  Factor loading matrix (P x K).

- psi:

  Vector of unique variances (P).

## Value

The sum of the log-density.
