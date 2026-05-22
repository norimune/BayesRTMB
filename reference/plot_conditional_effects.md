# Plot conditional effects

Calculate and plot conditional effects for a fitted model.

## Usage

``` r
plot_conditional_effects(fit, effect, ...)
```

## Arguments

- fit:

  A fitted model object (e.g., \`MCMC_Fit\`, \`MAP_Fit\`, \`VB_Fit\`).

- effect:

  Name of the variable to visualize (e.g., "X1" or "X1:X2").

- ...:

  Additional arguments passed to
  [`conditional_effects`](https://norimune.github.io/BayesRTMB/reference/conditional_effects.md)
  or the plot method.

## Value

Invisibly returns the plotted object of class `ce_rtmb`.
