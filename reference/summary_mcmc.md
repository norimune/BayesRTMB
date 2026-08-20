# Summarize MCMC Draws Stored in an Array

Summarizes posterior draws stored in the same array layout used by
[`MCMC_Fit`](https://norimune.github.io/BayesRTMB/reference/MCMC_Fit.md):
iterations by chains by variables. The returned object uses the same
columns and print method as `MCMC_Fit$summary()`.

## Usage

``` r
summary_mcmc(draws, pars = NULL, chains = NULL, max_rows = 10, digits = 2)
```

## Arguments

- draws:

  A numeric three-dimensional array with dimensions iterations by chains
  by variables.

- pars:

  Optional numeric or character vector selecting variables. Character
  values may be full variable names (for example, \`"beta\[1\]"\`) or
  base names (for example, \`"beta"\`). Prefix character names with
  \`"-"\` to exclude them.

- chains:

  Optional numeric vector selecting chains. Positive and negative
  integer indexing are supported.

- max_rows:

  Maximum number of variables to include. Use \`NULL\` to include all
  variables.

- digits:

  Number of decimal places used when printing the result.

## Value

A data frame with class \`"summary_BayesRTMB"\` containing posterior
means, standard deviations, marginal MAP estimates, 95 percent
intervals, bulk and tail effective sample sizes, and split R-hat values.

## Examples

``` r
set.seed(123)
draws <- array(
  rnorm(200 * 2 * 2),
  dim = c(200, 2, 2),
  dimnames = list(NULL, c("chain1", "chain2"), c("alpha", "beta"))
)
summary_mcmc(draws)
#> variable  mean    sd    map   q2.5  q97.5  ess_bulk  ess_tail  rhat 
#> alpha     0.02  0.97  -0.23  -1.75   2.04       434       426  1.00 
#> beta      0.00  0.99   0.06  -2.04   1.88       443       359  1.00 
```
