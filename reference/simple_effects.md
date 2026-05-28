# Calculate Simple Effects

Calculate the effect of a focal variable at different levels of a
moderator. For categorical focal variables, it calculates pairwise
differences (contrasts). For continuous focal variables, it calculates
the slope (simple slopes).

## Usage

``` r
simple_effects(fit, effect, prob = 0.95, sd_multiplier = 1, ...)
```

## Arguments

- fit:

  Model fit object (e.g., \`map_fit\`, \`mcmc_fit\`).

- effect:

  Character string of the interaction (e.g., "A:B"). The first variable
  is the focal variable.

- prob:

  Probability for the credible/confidence interval (default is 0.95).

- sd_multiplier:

  Multiplier for SD for continuous moderators (default is 1).

- ...:

  Additional arguments.

## Value

A \`ce_simple\` object (data frame) containing the estimated effects and
their credible intervals.

## Examples

``` r
# \donttest{
  data(debate, package = "BayesRTMB")
  fit <- rtmb_lm(sat ~ talk * perf, data = debate)
#> Pre-checking model code...
#> Checking RTMB setup...
  map_fit <- fit$optimize()
#> Starting RTMB optimization...
#> 
  # Effect of talk at different levels of performance
  se <- simple_effects(map_fit, effect = "talk:perf")
  print(se)
#> --- Simple Effects Analysis ---
#>  moderator  perf          term estimate lower upper
#>       perf 2.930 Slope of talk    0.037 0.037 0.037
#>       perf 4.690 Slope of talk    0.266 0.266 0.266
#>       perf 6.450 Slope of talk    0.495 0.495 0.495
#> 
# }
```
