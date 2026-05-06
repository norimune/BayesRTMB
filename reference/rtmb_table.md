# RTMB-based Contingency Table Analysis (Chi-squared Test)

\`rtmb_table\` performs a chi-squared test of independence between two
categorical variables. It provides both classic (frequentist) Pearson
chi-squared tests and Bayesian multinomial-style models.

## Usage

``` r
rtmb_table(
  x,
  y = NULL,
  data = NULL,
  classic = FALSE,
  correct = TRUE,
  prior = prior_uniform(),
  ...
)
```

## Arguments

- x:

  Variable name or formula.

- y:

  Variable name (optional if x is a formula).

- data:

  A data frame.

- classic:

  Logical; if TRUE, perform frequentist chi-squared and Fisher's exact
  tests.

- correct:

  Logical; if TRUE, apply Yates' continuity correction (for 2x2 classic
  only).

- prior:

  Prior specification (Bayesian mode). Default is \`prior_uniform()\`.

- ...:

  Additional arguments.

## Value

A \`Classic_Fit\` or \`MCMC_Fit\` object.

## Examples

``` r
# \donttest{
# Classic chi-squared test
rtmb_table(skill, cond, data = debate, classic = TRUE)
#> <Classic_Model>
#>   Public:
#>     .perform_fit: function (data) 
#>     .resample_data: function () 
#>     clone: function (deep = FALSE) 
#>     data: data.frame
#>     estimate: function (bootstrap = FALSE, n_boot = 1000) 
#>     extra: list
#>     family: gaussian
#>     formula: NULL
#>     initialize: function (type, formula, data, family = "gaussian", view = NULL, 
#>     obj: NULL
#>     refit_fn: NULL
#>     type: table
#>     view: NULL
# }
```
