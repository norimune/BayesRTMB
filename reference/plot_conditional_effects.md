# Plot conditional effects

Calculate and plot conditional effects for a fitted model.

## Usage

``` r
plot_conditional_effects(
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

  A fitted model object (e.g., \`MCMC_Fit\`, \`MAP_Fit\`, \`VB_Fit\`).

- effect:

  Name of the variable to visualize (e.g., "X1" or "X1:X2").

- prob:

  Probability for the credible/confidence interval.

- sd_multiplier:

  Numeric multiplier for standard deviation when splitting continuous
  moderators.

- sd_slice:

  Logical or NULL; controls whether continuous moderators are evaluated
  at mean - SD, mean, and mean + SD.

- resolution:

  Grid resolution to calculate for continuous variables.

- ...:

  Additional arguments passed to the plot method.

## Value

Invisibly returns the plotted object of class `ce_rtmb`.
