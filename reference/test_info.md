# Calculate Test Information Function

Calculate Test Information Function

## Usage

``` r
test_info(x, ...)
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
  ti <- test_info(map_fit)
  plot(ti)
} # }
```
