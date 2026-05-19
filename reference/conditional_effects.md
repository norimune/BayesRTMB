# Calculate Conditional Effects

Calculate and visualize the predicted values (marginal effects) of a
variable, potentially conditional on the levels of another variable
(interaction).

## Usage

``` r
conditional_effects(fit, effect, prob = 0.95, sd_multiplier = 1, ...)
```

## Arguments

- fit:

  Model fit object (e.g., \`mcmc_fit\`, \`map_fit\`).

- effect:

  Name of the variable to visualize (e.g., "X1" or "X1:X2").

- prob:

  Probability for the credible/confidence interval (default is 0.95).

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
  data(debate, package = "BayesRTMB")
  fit <- rtmb_lm(sat ~ talk * perf, data = debate)
  mcmc_fit <- fit$sample()
  ce <- conditional_effects(mcmc_fit, effect = "talk:perf")
  plot(ce)
  summary(ce)
} # }
```
