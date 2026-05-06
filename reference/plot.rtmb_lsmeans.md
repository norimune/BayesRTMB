# Plot marginal means with error bars

## Usage

``` r
# S3 method for class 'rtmb_lsmeans'
plot(
  x,
  y,
  error_bar = "se",
  col = NULL,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  ...
)
```

## Arguments

- x:

  An rtmb_lsmeans object.

- y:

  Ignored.

- error_bar:

  Type of error bar. Either se (Standard Error) or ci (95

  colColor of the bars.

  mainPlot title.

  xlabX axis label.

  ylabY axis label.

  ...Additional arguments passed to barplot().

The bar centers (invisibly). plot method for rtmb_lsmeans objects.
Creates a bar plot of marginal means with error bars.
