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

- resolution:

  Grid resolution to calculate for continuous variables (default is
  100).

- ...:

  Additional arguments.

## Value

A \`ce_rtmb\` object.
