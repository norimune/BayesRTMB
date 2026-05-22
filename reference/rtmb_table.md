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
  correct = TRUE,
  prior = prior_flat(),
  fixed = NULL,
  WAIC = FALSE,
  ...
)
```

## Arguments

- x:

  Variable name, formula, table, or matrix.

- y:

  Variable name (optional if x is a formula).

- data:

  A data frame.

- correct:

  Logical; if TRUE, apply Yates' continuity correction for 2x2 classic
  analyses.

- prior:

  Prior specification (Bayesian mode). Default is \`prior_flat()\`.

- fixed:

  Optional named list of fixed values for specific parameters.

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

- ...:

  Additional arguments.

## Value

An \`RTMB_Model\` object.

## Examples

``` r
# \donttest{
# Classic chi-squared test
rtmb_table(skill, cond, data = debate)$classic()
#> Pre-checking model code...
#> Checking RTMB setup...
#> Starting RTMB optimization...
#> 
#> 
#> Contingency Table Analysis
#> 
#> Log-Likelihood: -8.673, AIC: 29.346, BIC: 17.346
#> 
#> ---
#> Cell Probabilities (p) and Confidence Intervals:
#>                    Estimate Std. Error Lower 95% Upper 95% z value     Pr
#> p[skill:1, cond:0]  0.14333    0.02023   0.10368   0.18299 7.08479 <.0001
#> p[skill:2, cond:0]  0.17667    0.02202   0.13351   0.21982 8.02322 <.0001
#> p[skill:3, cond:0]  0.18000    0.02218   0.13653   0.22347 8.11501 <.0001
#> p[skill:1, cond:1]  0.14333    0.02023   0.10368   0.18299 7.08480 <.0001
#> p[skill:2, cond:1]  0.15000    0.02062   0.10959   0.19041 7.27604 <.0001
#> p[skill:3, cond:1]  0.20667    0.02338   0.16085   0.25249 8.84032 <.0001
#> 
#> Expected Counts (Independence) and Pearson Residuals:
#>                    Expected Residual
#> E[skill:1, cond:0] 43.00002  0.00000
#> E[skill:2, cond:0] 48.99992  0.57144
#> E[skill:3, cond:0] 58.00005 -0.52523
#> E[skill:1, cond:1] 43.00002  0.00000
#> E[skill:2, cond:1] 48.99992 -0.57142
#> E[skill:3, cond:1] 58.00005  0.52522
rtmb_table(table(debate$skill, debate$cond))$classic()
#> Pre-checking model code...
#> Checking RTMB setup...
#> Starting RTMB optimization...
#> 
#> 
#> Contingency Table Analysis
#> 
#> Log-Likelihood: -8.673, AIC: 29.346, BIC: 17.346
#> 
#> ---
#> Cell Probabilities (p) and Confidence Intervals:
#>                                                                                   Estimate
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:0]  0.14333
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:0]  0.17667
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:0]  0.18000
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:1]  0.14333
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:1]  0.15000
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:1]  0.20667
#>                                                                                   Std. Error
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:0]    0.02023
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:0]    0.02202
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:0]    0.02218
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:1]    0.02023
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:1]    0.02062
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:1]    0.02338
#>                                                                                   Lower 95%
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:0]   0.10368
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:0]   0.13351
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:0]   0.13653
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:1]   0.10368
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:1]   0.10959
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:1]   0.16085
#>                                                                                   Upper 95%
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:0]   0.18299
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:0]   0.21982
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:0]   0.22347
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:1]   0.18299
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:1]   0.19041
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:1]   0.25249
#>                                                                                   z value
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:0] 7.08478
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:0] 8.02323
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:0] 8.11502
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:1] 7.08479
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:1] 7.27606
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:1] 8.84030
#>                                                                                       Pr
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:0] <.0001
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:0] <.0001
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:0] <.0001
#> p[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:1] <.0001
#> p[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:1] <.0001
#> p[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:1] <.0001
#> 
#> Expected Counts (Independence) and Pearson Residuals:
#>                                                                                   Expected
#> E[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:0] 42.99998
#> E[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:0] 49.00007
#> E[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:0] 58.00008
#> E[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:1] 42.99990
#> E[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:1] 48.99999
#> E[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:1] 57.99998
#>                                                                                   Residual
#> E[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:0]  0.00000
#> E[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:0]  0.57142
#> E[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:0] -0.52524
#> E[table(debate$skill, debate$cond)_row:1, table(debate$skill, debate$cond)_col:1]  0.00001
#> E[table(debate$skill, debate$cond)_row:2, table(debate$skill, debate$cond)_col:1] -0.57143
#> E[table(debate$skill, debate$cond)_row:3, table(debate$skill, debate$cond)_col:1]  0.52523
# }
```
