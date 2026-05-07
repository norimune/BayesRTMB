# Plot marginal means with error bars

Plot method for rtmb_lsmeans objects.

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

  Type of error bar ("se" or "ci").

- col:

  Color of the bars.

- main:

  Plot title.

- xlab:

  X axis label.

- ylab:

  Y axis label.

- ...:

  Additional arguments passed to barplot.
