# Calculate conditional effects for MCMC fit objects

Calculate conditional effects for MCMC fit objects

## Usage

``` r
# S3 method for class 'mcmc_fit'
conditional_effects(
  fit,
  effect,
  prob = 0.95,
  sd_multiplier = 1,
  sd_slice = NULL,
  resolution = 100,
  ...
)
```

## Arguments

- fit:

  An object of class \`MCMC_Fit\`.

- effect:

  Name of the explanatory variable to visualize (e.g., "X1" or "X1:X2").

- prob:

  Probability for the credible/confidence interval (default is 0.95).

- sd_multiplier:

  Numeric. Multiplier for standard deviation when splitting continuous
  moderators (default is 1).

- sd_slice:

  Logical or NULL. If TRUE, continuous moderators are evaluated at
  mean - SD, mean, and mean + SD. If FALSE, all observed moderator
  values are used. If NULL (default), sd slicing is used automatically
  when the moderator has 6 or more unique values.

- resolution:

  Grid resolution to calculate for continuous variables (default is
  100).

- ...:

  Additional arguments.

## Value

A \`ce_rtmb\` object.
