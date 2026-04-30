# Calculate Simple Effects

Calculate the effect of a focal variable at different levels of a
moderator. For categorical focal variables, it calculates pairwise
differences (contrasts). For continuous focal variables, it calculates
the slope (simple slopes).

## Usage

``` r
simple_effects(fit, effect, sd_multiplier = 1, ...)
```

## Arguments

- fit:

  Model fit object (mcmc_fit).

- effect:

  Character string of the interaction (e.g., "A:B"). The first variable
  is the focal variable.

- sd_multiplier:

  Multiplier for SD for continuous moderators (default is 1).

- ...:

  Additional arguments.

## Value

A \`ce_simple\` object (data frame) containing the estimated effects and
their credible intervals.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- rtmb_lm(len ~ supp * dose, data = ToothGrowth)
  mcmc_fit <- fit$sample()
  # Effect of supplement at each dose level
  se <- simple_effects(mcmc_fit, effect = "supp:dose")
  print(se)
} # }
```
