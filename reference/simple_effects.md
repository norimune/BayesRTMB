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
  mcmc_fit <- fit$sample()
#> Starting sequential sampling (chains = 4)...
#> chain 1 started... 
#> chain 1: iter 200 warmup 
#> chain 1: iter 400 warmup 
#> chain 1: iter 600 warmup 
#> chain 1: iter 800 warmup 
#> chain 1: iter 1000 warmup 
#> chain 1: iter 1200 sampling 
#> chain 1: iter 1400 sampling 
#> chain 1: iter 1600 sampling 
#> chain 1: iter 1800 sampling 
#> chain 1: iter 2000 sampling 
#> chain 2 started... 
#> chain 2: iter 200 warmup 
#> chain 2: iter 400 warmup 
#> chain 2: iter 600 warmup 
#> chain 2: iter 800 warmup 
#> chain 2: iter 1000 warmup 
#> chain 2: iter 1200 sampling 
#> chain 2: iter 1400 sampling 
#> chain 2: iter 1600 sampling 
#> chain 2: iter 1800 sampling 
#> chain 2: iter 2000 sampling 
#> chain 3 started... 
#> chain 3: iter 200 warmup 
#> chain 3: iter 400 warmup 
#> chain 3: iter 600 warmup 
#> chain 3: iter 800 warmup 
#> chain 3: iter 1000 warmup 
#> chain 3: iter 1200 sampling 
#> chain 3: iter 1400 sampling 
#> chain 3: iter 1600 sampling 
#> chain 3: iter 1800 sampling 
#> chain 3: iter 2000 sampling 
#> chain 4 started... 
#> chain 4: iter 200 warmup 
#> chain 4: iter 400 warmup 
#> chain 4: iter 600 warmup 
#> chain 4: iter 800 warmup 
#> chain 4: iter 1000 warmup 
#> chain 4: iter 1200 sampling 
#> chain 4: iter 1400 sampling 
#> chain 4: iter 1600 sampling 
#> chain 4: iter 1800 sampling 
#> chain 4: iter 2000 sampling 
  # Effect of talk at different levels of performance
  se <- simple_effects(mcmc_fit, effect = "talk:perf")
  print(se)
#> --- Simple Effects Analysis ---
#>  moderator  perf          term estimate  lower upper
#>       perf 2.930 Slope of talk    0.034 -0.103 0.172
#>       perf 4.690 Slope of talk    0.264  0.171 0.354
#>       perf 6.450 Slope of talk    0.494  0.358 0.628
#> 
# }
```
