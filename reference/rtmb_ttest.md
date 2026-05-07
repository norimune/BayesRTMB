# RTMB-based Bayesian two-sample t-test wrapper function

Performs a Bayesian or Frequentist two-sample t-test using RTMB.

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
  fixed = NULL,
  null = NULL,
  var.equal = TRUE,
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

  An object of class "rtmb_prior" (e.g., \`prior_uniform()\`,
  \`prior_jzs()\`, or \`prior_weak()\`).

- init:

  List of initial values.

- null:

  Character string specifying the target parameter for the null model
  (e.g., "delta").

- var.equal:

  Logical; whether to assume equal variances. Default is TRUE.

- ...:

  Additional arguments.

## Value

An \`RTMB_Model\` object.
