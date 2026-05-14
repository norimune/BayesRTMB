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
  prior = prior_flat(),
  init = NULL,
  fixed = NULL,
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

  An object of class \`"rtmb_prior"\`. Use \`prior_flat()\` for no
  prior, \`prior_normal()\` for default normal/exponential priors, or
  \`prior_weak()\` for weakly informative Bayesian inference. Default is
  \`prior_flat()\`.

- init:

  List of initial values.

- fixed:

  Optional named list of fixed values for specific parameters.

- var.equal:

  Logical; whether to assume equal variances. Default is TRUE.

- ...:

  Additional arguments.

## Value

An \`RTMB_Model\` object.

## Details

For classic inference, heteroscedastic two-sample t-tests use the same
RTMB Satterthwaite machinery as \`optimize(marginal = ..., df_method =
"satterthwaite")\`. The result is model-based and reproducible from the
printed model code. This corresponds to the Welch-type unequal-variance
t-test, but the degrees of freedom are computed by the package's
internal Satterthwaite procedure rather than by a separate closed-form
formula.
