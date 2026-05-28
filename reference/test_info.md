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
# \donttest{
  fit <- rtmb_irt(data = BigFive[, 1:5], model = "2PL", type = "ordered")
#> Pre-checking model code...
#> Checking RTMB setup...
  map_fit <- fit$optimize()
#> Starting RTMB optimization...
#> 
  ti <- test_info(map_fit)
  plot(ti)

# }
```
