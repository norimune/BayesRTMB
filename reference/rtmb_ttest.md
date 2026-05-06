# RTMB-based Bayesian two-sample t-test wrapper function

Performs a Bayesian two-sample t-test using RTMB. It estimates the
effect size (delta) with a Cauchy prior, allowing for robust inference
and calculation of Bayes factors.

## Usage

``` r
rtmb_ttest(
  x,
  y = NULL,
  data = NULL,
  r = 0.707,
  paired = FALSE,
  ID = NULL,
  y_range = NULL,
  prior = prior_uniform(),
  init = NULL,
  null = NULL,
  classic = FALSE,
  var.equal = FALSE,
  ...
)
```

## Arguments

- x:

  Numeric vector of responses for group 1, a formula (e.g., \`y ~
  group\`), or a column name (unquoted) if \`data\` is provided.

- y:

  Numeric vector of responses for group 2, or a column name (unquoted)
  if \`data\` is provided. Required if \`x\` is not a formula.

- data:

  Data frame containing the variables.

- r:

  Numeric; Cauchy prior scale for the effect size (delta). Default is
  0.707.

- paired:

  Logical; whether to perform a paired t-test.

- ID:

  Character; name of the ID variable for paired t-tests (required for
  formula input with paired = TRUE).

- y_range:

  Theoretical minimum and maximum values of the response variable as a
  vector c(min, max). Required when using weakly informative priors.

- prior:

  An object of class "rtmb_prior" (e.g., \`prior_uniform()\` or
  \`prior_weak()\`).

- init:

  List of initial values.

- null:

  Character string specifying the target parameter for the null model
  (e.g., "delta" or "delta ~ cauchy(0, r)").

- classic:

  Logical; if TRUE, use classical (frequentist) estimation.

- var.equal:

  Logical; for independent t-tests in classic mode, whether to assume
  equal variances (Student's t-test) or not (Welch's t-test). Default is
  FALSE.

- ...:

  Additional arguments.

## Value

An `RTMB_Model` object.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Simulate two-sample data with a true effect size
  set.seed(123)
  y1 <- rnorm(30, mean = 0.5, sd = 1)
  y2 <- rnorm(30, mean = 0.0, sd = 1)

  # Fit the Bayesian two-sample t-test model
  # r = 0.707 is the standard scale for the Cauchy prior on the effect size
  fit_ttest <- rtmb_ttest(y1, y2, r = 0.707)

  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_ttest <- fit_ttest$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_ttest$summary()

  # Calculate Bayes factor against the null hypothesis (effect size delta = 0)
  # Specifying "delta" automatically fixes the parameter to 0 and drops its prior
  bf_ttest <- mcmc_ttest$bayes_factor(null_model = "delta")
  print(bf_ttest)
} # }
```
