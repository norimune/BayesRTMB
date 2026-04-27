# Item Response Curve for RTMB_Fit_Base

Item Response Curve for RTMB_Fit_Base

## Usage

``` r
# S3 method for class 'RTMB_Fit_Base'
item_curve(x, theta_seq = seq(-4, 4, length.out = 100), items = NULL, ...)
```

## Arguments

- x:

  An object of class RTMB_Fit_Base.

- theta_seq:

  Sequence of trait values (ability) to evaluate.

- items:

  Index or item names to restrict the calculation to specific items
  (optional).

- ...:

  Additional arguments.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- rtmb_irt(data = BigFive[, 1:5], model = "2pl")
  map_fit <- fit$optimize()
  ic <- item_curve(map_fit)
  plot(ic)
} # }
```
