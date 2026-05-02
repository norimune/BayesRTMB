# Specify a weakly informative prior

Specify a weakly informative prior

## Usage

``` r
prior_weak(sd_ratio = 0.5, max_beta = 1, ...)
```

## Arguments

- sd_ratio:

  Ratio of the prior standard deviation to the half-range of the
  response variable. Default is 0.5.

- max_beta:

  Maximum expected effect size. Default is 1.0.

- ...:

  Optional hyperparameters for other distributions.

## Value

A list with class "rtmb_prior"
