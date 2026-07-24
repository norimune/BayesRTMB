# RTMB-based Factor Analysis Wrapper

Performs exploratory factor analysis (EFA) using RTMB. It supports
standard rotation methods (e.g., varimax, promax) as well as regularized
factor analysis using a Spike-and-Slab Prior (SSP) for estimating sparse
loading matrices.

## Usage

``` r
rtmb_fa(
  data,
  nfactors = 1,
  rotate = NULL,
  score = FALSE,
  prior = prior_flat(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  missing = c("listwise", "fiml"),
  WAIC = FALSE
)
```

## Arguments

- data:

  Observation data frame or matrix (N x J).

- nfactors:

  Number of factors (K).

- rotate:

  String specifying the rotation method (e.g., "varimax", "promax",
  "ssp"). If NULL, no rotation is applied. Specifying "ssp" performs
  regularized factor analysis.

- score:

  Logical; if TRUE, factor scores are calculated in the generate block
  (default is FALSE).

- prior:

  Prior configuration: \`prior_flat()\`, \`prior_normal()\`, or
  \`prior_weak()\`. \`prior_ssp()\` is accepted when \`rotate = "ssp"\`.
  Hyperparameters can be specified within these functions (e.g.,
  \`prior_normal(mean_sd = 10, sd_rate = 1 / 5, loadings_sd = 10)\`).
  Available FA aliases are \`mean_sd\`, \`sd_rate\`, and
  \`loadings_sd\`; they override the corresponding common arguments
  \`mu_sd\`, \`sigma_rate\`, and \`b_sd\`. \`sd_rate\` is an exponential
  rate, so \`1 / 5\` gives a prior mean of 5. \`ssp_ratio\` is available
  when \`rotate = "ssp"\`.

- y_range:

  A numeric vector of length 2 specifying the theoretical min and max
  values of the items. Used to construct weakly informative priors when
  \`prior = prior_weak()\`.

- init:

  List of initial values.

- fixed:

  A named list of parameter values to fix (optional).

- missing:

  Missing value handling strategy: "listwise" (default) or "fiml" (Full
  Information Maximum Likelihood).

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

## Examples

``` r

  # Prepare a subset of the BigFive dataset for factor analysis
  data(BigFive, package = "BayesRTMB")
  fa_data <- BigFive[, 1:10]

  # --- 1. Standard Exploratory Factor Analysis (1 Factor) ---
  fit_fa1 <- rtmb_fa(data = fa_data, nfactors = 1)
#> Pre-checking model code...
#> Checking RTMB setup...

  # Maximum A Posteriori (MAP) estimation
  map_fa1 <- fit_fa1$optimize()
#> Starting RTMB optimization...
  map_fa1$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 2485.47
#> Approx. Log Marginal Likelihood (Laplace): -2535.33
#> 
#> Point Estimates and 95% Wald CI:
#>        variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L[BF1,Factor1]    0.05250     0.08868   -0.12130    0.22630 
#> L[BF2,Factor1]    0.79972     0.07145    0.65968    0.93976 
#> L[BF3,Factor1]   -0.13478     0.08707   -0.30542    0.03587 
#> L[BF4,Factor1]   -0.20465     0.08257   -0.36649   -0.04282 
#> L[BF5,Factor1]   -0.14443     0.08455   -0.31015    0.02128 
#> L[BF6,Factor1]   -0.00682     0.08660   -0.17655    0.16291 
#> L[BF7,Factor1]    0.84338     0.07318    0.69994    0.98681 
#> L[BF8,Factor1]   -0.24816     0.08298   -0.41080   -0.08551 
#> L[BF9,Factor1]   -0.11264     0.08648   -0.28214    0.05686 
#> L[BF10,Factor1]   0.21358     0.08794    0.04122    0.38593 
#> 
```
