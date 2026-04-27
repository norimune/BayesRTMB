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
  y_range = NULL,
  use_weak_info = FALSE,
  prior = list(mean_sd = 10, sd_rate = 0.1),
  weak_info_prior = list(sd_ratio = 0.5),
  init = NULL,
  null = NULL,
  y1 = NULL,
  y2 = NULL
)
```

## Arguments

- x:

  Numeric vector of responses for group 1, or a formula (e.g., \`y ~
  group\`).

- y:

  Numeric vector of responses for group 2. Required if \`x\` is not a
  formula.

- data:

  Data frame containing the variables in the formula.

- r:

  Numeric; Cauchy prior scale for the effect size (delta). Default is
  0.707.

- y_range:

  Theoretical minimum and maximum values of the response variable as a
  vector c(min, max). Specifying this automatically enables weakly
  informative priors.

- use_weak_info:

  Logical; whether to explicitly use weakly informative priors.

- prior:

  List of hyperparameters for the default fixed priors.

- weak_info_prior:

  List of hyperparameters for the weakly informative priors.

- init:

  List of initial values.

- null:

  Character string specifying the target parameter for the null model
  (e.g., "delta" or "delta ~ cauchy(0, r)").

- y1:

  Deprecated. Use \`x\` instead.

- y2:

  Deprecated. Use \`y\` instead.

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
