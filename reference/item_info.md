# Calculate Item Information Function

Calculate Item Information Function

## Usage

``` r
item_info(x, ...)
```

## Arguments

- x:

  An object of class RTMB_Fit_Base

- ...:

  Additional arguments.

## Examples

``` r
# \donttest{
  fit <- rtmb_irt(data = BigFive[, 1:5], model = "2PL")
#> Error: Binary IRT models (type = 'binary') require responses to be exactly 0, 1, or NA. Found other values.
  map_fit <- fit$optimize()
#> Error: object 'fit' not found
  ii <- item_info(map_fit)
#> Error: object 'map_fit' not found
  plot(ii)
#> Error: object 'ii' not found
# }
```
