# Model Code Wrapper for RTMB

Model Code Wrapper for RTMB

## Usage

``` r
model_code(expr, env = parent.frame())
```

## Arguments

- expr:

  A block of code containing model description.

- env:

  Environment to assign to the generated function.

## Value

A standard R function object taking (dat, par).
