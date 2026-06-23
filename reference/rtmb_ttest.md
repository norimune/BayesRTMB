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
  prior, \`prior_normal()\` for default normal/exponential priors,
  \`prior_jzs()\` for JZS t-test priors, or \`prior_weak()\` for weakly
  informative Bayesian inference. Default is \`prior_flat()\`.

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

  Reserved; unused arguments are rejected.

## Value

An \`RTMB_Model\` object.

## Details

For classic inference, heteroscedastic two-sample t-tests use the same
RTMB Satterthwaite machinery as \`optimize(marginal = ..., df_method =
"satterthwaite")\`. The result is model-based and reproducible from the
printed model code. This corresponds to the Welch-type unequal-variance
t-test, but the degrees of freedom are computed by the package's
internal Satterthwaite procedure rather than by a separate closed-form
formula. When \`prior_jzs()\` is combined with \`var.equal = FALSE\`,
BayesRTMB uses a Welch-style JZS extension: the effect size \`delta\` is
an explicit parameter with a Cauchy prior, and the group mean difference
is scaled by the root-mean-square of the two group standard deviations.

## Examples

``` r

  # Simulate two-sample data with a true effect size
  set.seed(123)
  y1 <- rnorm(30, mean = 0.5, sd = 1)
  y2 <- rnorm(30, mean = 0.0, sd = 1)

  # Fit the Bayesian two-sample t-test model
  # r = 0.707 is the standard scale for the Cauchy prior on the effect size
  fit_ttest <- rtmb_ttest(y1, y2, prior = prior_jzs(r = 0.707))
#> Pre-checking model code...
#> Checking RTMB setup...

  # \donttest{
  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_ttest <- fit_ttest$sample(sampling = 500, warmup = 500, chains = 2)
#> Starting sequential sampling (chains = 2)...
#> chain 1 started...
#> chain 1: iter 200/1000 (20%) warmup
#> chain 1: iter 400/1000 (40%) warmup
#> chain 1: iter 600/1000 (60%) sampling
#> chain 1: iter 800/1000 (80%) sampling
#> chain 1: iter 1000/1000 (100%) sampling
#> chain 1 done (100%)
#> chain 2 started...
#> chain 2: iter 200/1000 (20%) warmup
#> chain 2: iter 400/1000 (40%) warmup
#> chain 2: iter 600/1000 (60%) sampling
#> chain 2: iter 800/1000 (80%) sampling
#> chain 2: iter 1000/1000 (100%) sampling
#> chain 2 done (100%)
#> sampling: 100%
#> Calculating transformed parameters...
#> transformed parameters: 20%
#> transformed parameters: 40%
#> transformed parameters: 60%
#> transformed parameters: 80%
#> transformed parameters: 100%
  mcmc_ttest$summary()
#>   variable    mean    sd     map    q2.5   q97.5  ess_bulk  ess_tail  rhat 
#> lp          -81.01  1.27  -80.20  -84.28  -79.67       521       586  1.00 
#> diff          0.24  0.23    0.16   -0.21    0.71       759       463  1.00 
#> delta         0.26  0.25    0.17   -0.22    0.78       723       445  1.00 
#> total_mean    0.32  0.11    0.32    0.10    0.54       854       618  1.00 
#> sd            0.93  0.08    0.89    0.78    1.10       952       779  1.00 
#> mean0         0.44  0.16    0.47    0.12    0.79       802       602  1.00 
#> mean1         0.20  0.16    0.22   -0.13    0.50       748       604  1.00 
  # }
```
