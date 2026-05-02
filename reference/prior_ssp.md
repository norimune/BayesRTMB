# Specify a Spike-and-Slab prior for variable selection

Specify a Spike-and-Slab prior for variable selection

## Usage

``` r
prior_ssp(ssp_ratio = 0.25, max_beta = 1, ...)
```

## Arguments

- ssp_ratio:

  Prior probability of inclusion for each variable. Default is 0.25.

- max_beta:

  Maximum expected effect size. Default is 1.0.

- ...:

  Optional hyperparameters for other distributions.

## Value

A list with class "rtmb_prior"
