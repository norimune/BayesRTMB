# Calculate Conditional Effects

Calculate Conditional Effects

## Usage

``` r
conditional_effects(fit, effect, ...)
```

## Arguments

- fit:

  Model fit object.

- effect:

  Name of the variable to visualize the effect.

- ...:

  Additional arguments.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- rtmb_lm(mpg ~ wt + hp, data = mtcars)
  mcmc_fit <- fit$sample(sampling = 500, warmup = 500)
  ce <- conditional_effects(mcmc_fit, effect = "wt")
  plot(ce)
} # }
```
