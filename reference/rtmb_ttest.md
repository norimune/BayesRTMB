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
  missing = c("listwise", "fiml"),
  WAIC = FALSE,
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

- missing:

  Missing value handling strategy: "listwise".

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

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

## Examples

``` r

  # Simulate two-sample data with a true effect size
  set.seed(123)
  y1 <- rnorm(30, mean = 0.5, sd = 1)
  y2 <- rnorm(30, mean = 0.0, sd = 1)

  # Fit the Bayesian two-sample t-test model
  # r = 0.707 is the standard scale for the Cauchy prior on the effect size
  fit_ttest <- rtmb_ttest(y1, y2, r = 0.707)
#> Pre-checking model code...
#> Checking RTMB setup...

  if (FALSE) { # \dontrun{
  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_ttest <- fit_ttest$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_ttest$summary()

  # Calculate Bayes factor against the null hypothesis (effect size delta = 0)
  # Specifying fixed = list(delta = 0) compares against a model with delta fixed to 0
  bf_ttest <- mcmc_ttest$bayes_factor(fixed = list(delta = 0))
  print(bf_ttest)
  } # }
```
