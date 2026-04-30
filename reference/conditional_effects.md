# Calculate Conditional Effects

Calculate and visualize the predicted values (marginal effects) of a
variable, potentially conditional on the levels of another variable
(interaction).

## Usage

``` r
conditional_effects(fit, effect, sd_multiplier = 1, ...)
```

## Arguments

- fit:

  Model fit object (e.g., \`mcmc_fit\`).

- effect:

  Name of the variable to visualize (e.g., "X1" or "X1:X2").

- sd_multiplier:

  Numeric. Multiplier for standard deviation when splitting continuous
  moderators (default is 1).

- ...:

  Additional arguments.

## Value

An object of class \`ce_rtmb\` containing the predicted values and their
credible intervals.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- rtmb_lm(mpg ~ wt * hp, data = mtcars)
  mcmc_fit <- fit$sample()
  ce <- conditional_effects(mcmc_fit, effect = "wt:hp")
  plot(ce)
  summary(ce)
} # }
```
