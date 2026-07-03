# Exponentially modified normal log-probability density function

Log-density for \\X = Z + E\\, where \\Z\\ follows a normal distribution
with mean `mu` and standard deviation `sigma`, and \\E\\ follows an
exponential distribution with rate `lambda`.

## Usage

``` r
exp_mod_normal_lpdf(x, mu, sigma, lambda, sum = TRUE)
```

## Arguments

- x:

  Vector of quantiles.

- mu:

  Mean of the normal component.

- sigma:

  Standard deviation of the normal component.

- lambda:

  Rate of the exponential component.

- sum:

  Logical; if `TRUE` (default), returns the sum of log-densities. If
  `FALSE`, returns element-wise log-densities.

## Value

The sum of the log-density or a vector of log-densities.
