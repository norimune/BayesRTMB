# Transformed Code Wrapper for RTMB

Transformed Code Wrapper for RTMB

## Usage

``` r
transform_code(expr, env = parent.frame(), ad_seed_name = NULL)
```

## Arguments

- expr:

  A block of code containing calculations for transformed parameters.

- env:

  Environment to assign to the generated function.

- ad_seed_name:

  Optional parameter name used internally to seed RTMB's AD type.

## Value

A function taking (dat, par) that returns a named list.
