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
if (FALSE) { # \dontrun{
  # Compare two models using Bayes Factor
  fit1 <- rtmb_lm(mpg ~ wt, data = mtcars)
  fit2 <- rtmb_lm(mpg ~ wt + hp, data = mtcars)
  mcmc1 <- fit1$sample(sampling = 500, warmup = 500)
  mcmc2 <- fit2$sample(sampling = 500, warmup = 500)
  bf <- bayes_factor(mcmc1, mcmc2)
  print(bf)
} # }
```
