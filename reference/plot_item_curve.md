# Plot item/category response curves

Calculate and plot item/category response curves for a fitted IRT model.

## Usage

``` r
plot_item_curve(fit, ...)
```

## Arguments

- fit:

  A fitted IRT model object.

- ...:

  Additional arguments passed to
  [`item_curve`](https://norimune.github.io/BayesRTMB/reference/item_curve.md)
  or the plot method.

## Value

Invisibly returns the plotted object of class `rtmb_item_curve`.
