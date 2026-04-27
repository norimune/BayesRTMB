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
if (FALSE) { # \dontrun{
  fit <- rtmb_irt(data = BigFive[, 1:5], model = "2pl")
  map_fit <- fit$optimize()
  ii <- item_info(map_fit)
  plot(ii)
} # }
```
