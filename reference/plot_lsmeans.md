# Plot least-squares marginal means

Calculate and plot marginal means for a fitted model.

## Usage

``` r
plot_lsmeans(fit, specs, ...)
```

## Arguments

- fit:

  A fitted model object (e.g. \`Classic_Fit\`).

- specs:

  Character vector specifying the variables for which marginal means are
  calculated.

- ...:

  Additional arguments passed to
  [`lsmeans`](https://norimune.github.io/BayesRTMB/reference/lsmeans.md)
  or the plot method.

## Value

Invisibly returns the plotted object.
