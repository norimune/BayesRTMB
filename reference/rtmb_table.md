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
  fixed = NULL,
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

- fixed:

  Optional named list of fixed values for specific parameters.

- ...:

  Additional arguments.

## Value

A \`Classic_Fit\` or \`MCMC_Fit\` object.

## Examples

``` r
# \donttest{
# Classic chi-squared test
rtmb_table(skill, cond, data = debate, classic = TRUE)
#> 
#> Contingency Table Analysis
#> 
#> ---
#> 
#> Test Results:
#>           Estimate     Pr      Statistic  
#> X-squared  1.20479 .54750    Chi-squared  
#>                 NA .54310 Fisher's Exact  
# }
```
