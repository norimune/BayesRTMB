# RTMB-based Log-linear analysis (Poisson regression)

Performs Bayesian or Frequentist log-linear analysis (Poisson
regression) on a contingency table or raw data.

## Usage

``` r
rtmb_loglinear(formula, data, prior = prior_flat(), fixed = NULL, ...)
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

## Examples

``` r
if (FALSE) { # \dontrun{
  # Create a contingency table
  tab <- matrix(c(10, 20, 30, 40), nrow = 2)
  dimnames(tab) <- list(A = c("A1", "A2"), B = c("B1", "B2"))
  
  # Fit a log-linear model (independence model: ~ A + B)
  fit_log <- rtmb_loglinear(~ A + B, data = tab)
  
  # MAP estimation
  map_log <- fit_log$optimize()
  map_log$summary()

  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_log <- fit_log$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_log$summary()
} # }
```
