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
# \donttest{
  fit <- rtmb_irt(data = BigFive[, 1:5], model = "2PL")
#> Error: Binary IRT models (type = 'binary') require responses to be exactly 0, 1, or NA. Found other values.
  map_fit <- fit$optimize()
#> Error: object 'fit' not found
  ic <- item_curve(map_fit)
#> Error: object 'map_fit' not found
  plot(ic)
#> Error: object 'ic' not found
# }
```
