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
  mcmc1 <- fit1$sample(sampling = 500, warmup = 500)
#> Starting sequential sampling (chains = 4)...
#> chain 1 started... 
#> chain 1: iter 200 warmup 
#> chain 1: iter 400 warmup 
#> chain 1: iter 600 sampling 
#> chain 1: iter 800 sampling 
#> chain 1: iter 1000 sampling 
#> chain 2 started... 
#> chain 2: iter 200 warmup 
#> chain 2: iter 400 warmup 
#> chain 2: iter 600 sampling 
#> chain 2: iter 800 sampling 
#> chain 2: iter 1000 sampling 
#> chain 3 started... 
#> chain 3: iter 200 warmup 
#> chain 3: iter 400 warmup 
#> chain 3: iter 600 sampling 
#> chain 3: iter 800 sampling 
#> chain 3: iter 1000 sampling 
#> chain 4 started... 
#> chain 4: iter 200 warmup 
#> chain 4: iter 400 warmup 
#> chain 4: iter 600 sampling 
#> chain 4: iter 800 sampling 
#> chain 4: iter 1000 sampling 
  mcmc2 <- fit2$sample(sampling = 500, warmup = 500)
#> Starting sequential sampling (chains = 4)...
#> chain 1 started... 
#> chain 1: iter 200 warmup 
#> chain 1: iter 400 warmup 
#> chain 1: iter 600 sampling 
#> chain 1: iter 800 sampling 
#> chain 1: iter 1000 sampling 
#> chain 2 started... 
#> chain 2: iter 200 warmup 
#> chain 2: iter 400 warmup 
#> chain 2: iter 600 sampling 
#> chain 2: iter 800 sampling 
#> chain 2: iter 1000 sampling 
#> chain 3 started... 
#> chain 3: iter 200 warmup 
#> chain 3: iter 400 warmup 
#> chain 3: iter 600 sampling 
#> chain 3: iter 800 sampling 
#> chain 3: iter 1000 sampling 
#> chain 4 started... 
#> chain 4: iter 200 warmup 
#> chain 4: iter 400 warmup 
#> chain 4: iter 600 sampling 
#> chain 4: iter 800 sampling 
#> chain 4: iter 1000 sampling 
  bf <- bayes_factor(mcmc1, mcmc2)
#> Calculating marginal likelihood for Model 1...
#> Bridge Sampling Converged: LogML = -414.976 (Error = 0.0044, ESS = 383.1)
#> Calculating marginal likelihood for Model 2...
#> Bridge Sampling Converged: LogML = -404.592 (Error = 0.0053, ESS = 326.0)
  print(bf)
#> Bayes Factor (BF12) : 0 
#> Log Bayes Factor    : -10.3844 (Approx. Error = 0.0069)
#> Evidence            : Decisive evidence for Model 2 
# }
```
