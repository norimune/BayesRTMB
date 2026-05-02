# Specify a Regularized Horseshoe prior for continuous shrinkage

Specify a Regularized Horseshoe prior for continuous shrinkage

## Usage

``` r
prior_rhs(expected_vars = 3, slab_scale = 2, slab_df = 4, ...)
```

## Arguments

- expected_vars:

  Expected number of non-zero variables. Default is 3.

- slab_scale:

  Scale parameter for the slab distribution. Default is 2.0.

- slab_df:

  Degrees of freedom for the slab distribution. Default is 4.0.

- ...:

  Optional hyperparameters for other distributions.

## Value

A list with class "rtmb_prior"
