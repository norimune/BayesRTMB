# Calculate Bayes factor from log marginal likelihoods

Calculate Bayes factor from log marginal likelihoods

## Usage

``` r
bayes_factor(logml1, logml2)
```

## Arguments

- logml1:

  Log marginal likelihood of Model 1 (e.g., target model)

- logml2:

  Log marginal likelihood of Model 2 (e.g., reference/null model)

## Value

An object containing Bayes factor, log Bayes factor, estimation error,
and interpretation

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
