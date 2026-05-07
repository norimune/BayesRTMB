# RTMB-based Log-linear Analysis (Contingency Table)

\`rtmb_loglinear\` performs log-linear analysis of contingency tables
using RTMB. It models cell frequencies using a Poisson distribution and
a log link function. For Bayesian inference, it supports weakly
informative priors (Stan-style). For frequentist inference (\`classic =
TRUE\`), it provides Wald-based ANOVA for independence tests.

## Usage

``` r
rtmb_loglinear(formula, data, prior = prior_uniform(), ...)
```

## Arguments

- formula:

  A formula describing the table structure (e.g., \`~ A + B\` for
  independence, \`~ A \* B\` for dependence).

- data:

  A data frame, table, or matrix.

- prior:

  Prior specification (e.g., \`prior_uniform()\` or \`prior_weak()\`).
  Default is \`prior_uniform()\`.

- ...:

  Additional arguments passed to the internal estimation engine.

- classic:

  Logical; if TRUE, perform frequentist estimation and ANOVA tests.

- fisher:

  Logical; if TRUE, perform Fisher's exact test (for 2D tables).

## Value

An \`RTMB_Model\`, \`MCMC_Fit\`, or \`Classic_Fit\` object depending on
the settings.

## Examples

``` r
# \donttest{
# 2x2 table analysis
df <- as.data.frame(Titanic)
fit <- rtmb_loglinear(Survived ~ Sex + Age, data = df, classic = TRUE)
#> Error in rtmb_glmer(formula = formula, data = data, family = "poisson",     prior = prior, generate = gen_block, ...): unused argument (classic = TRUE)
fit$anova()
#> Error: object 'fit' not found

# Bayesian 4-way table interaction analysis
fit_bayas <- rtmb_loglinear(~ Class * Sex * Age * Survived, data = Titanic, prior = prior_weak())
#> Pre-checking model code...
#> Checking RTMB setup...
# }
```
