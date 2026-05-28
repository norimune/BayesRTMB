# Calculate Bayes Factor

Compare two models by calculating the Bayes Factor based on their
marginal likelihoods.

## Usage

``` r
bayes_factor(logml1, logml2, error_threshold = 0.2)
```

## Arguments

- logml1:

  The first model fit (e.g., \`mcmc_fit\`) or its log-marginal
  likelihood.

- logml2:

  The second model fit (e.g., \`mcmc_fit\`) or its log-marginal
  likelihood.

- error_threshold:

  Threshold for warning about high bridge-sampling error (default 0.2).

## Value

An object of class \`bayes_factor\` containing the Bayes factor, log
Bayes factor, approximate estimation error, and interpretation.

## Examples

``` r
# \donttest{
  # Compare two models using Bayes Factor
  data(debate, package = "BayesRTMB")
  fit1 <- rtmb_lm(sat ~ talk, data = debate)
#> Pre-checking model code...
#> Checking RTMB setup...
  fit2 <- rtmb_lm(sat ~ talk + perf, data = debate)
#> Pre-checking model code...
#> Checking RTMB setup...
  map1 <- fit1$optimize()
#> Starting RTMB optimization...
#> 
  map2 <- fit2$optimize()
#> Starting RTMB optimization...
#> 
  bf <- bayes_factor(map1, map2)
  print(bf)
#> Bayes Factor (BF12) : 0 
#> Log Bayes Factor    : -10.3811
#> Evidence            : Decisive evidence for Model 2 
# }
```
