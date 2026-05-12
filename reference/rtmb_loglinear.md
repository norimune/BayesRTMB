# RTMB-based Log-linear analysis (Poisson regression)

Performs Bayesian or Frequentist log-linear analysis (Poisson
regression) on a contingency table or raw data.

## Usage

``` r
rtmb_loglinear(formula, data, prior = prior_uniform(), fixed = NULL, ...)
```

## Arguments

- formula:

  A formula (e.g., \`~ A + B + A:B\`) or a contingency table.

- data:

  A data frame (required if \`formula\` is used).

- prior:

  An object of class "rtmb_prior" specifying the prior distribution.

- fixed:

  Optional named list of fixed values for specific parameters.

- ...:

  Additional arguments passed to \`rtmb_glm()\`.

## Value

An \`RTMB_Model\` object.
