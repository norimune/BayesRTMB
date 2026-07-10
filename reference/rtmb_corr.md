# Fit a Correlation Model using RTMB

\`rtmb_corr\` fits a correlation model to estimate means, standard
deviations, and correlation structures. It supports simple correlation,
multilevel correlation, and classical frequentist estimation.

## Usage

``` r
rtmb_corr(
  x = NULL,
  data = NULL,
  ID = NULL,
  covariates = NULL,
  method = c("pearson", "spearman", "reml"),
  prior = prior_flat(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  missing = c("listwise", "fiml", "pairwise"),
  WAIC = FALSE
)
```

## Arguments

- x:

  A matrix, data frame, formula, or expression (e.g., `cbind(V1, V2)`)
  of response variables.

- data:

  An optional data frame containing the variables.

- ID:

  A character string or expression specifying the group ID variable for
  multilevel models.

- covariates:

  Optional numeric matrix or data frame of covariates to be included in
  the joint MVN model.

- method:

  Correlation method for `classic()`: `"pearson"`, `"spearman"`, or
  `"reml"`.

- prior:

  Prior configuration object: \`prior_flat()\`, \`prior_normal()\`, or
  \`prior_weak()\`. Default is \`prior_flat()\`.

- y_range:

  Optional numeric vector or matrix defining the theoretical range (min,
  max) of response variables. Required when using
  [`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md).
  Can be a vector of length 2 (applies to all variables) or a
  matrix/list of length P.

- init:

  Optional list of initial values.

- fixed:

  Optional named list of fixed values for specific parameters.

- missing:

  Missing value handling strategy: "listwise" (default), "pairwise", or
  "fiml" (Full Information Maximum Likelihood).

- WAIC:

  Logical; if TRUE, add pointwise \`log_lik\` to the generate block for
  WAIC.

## Examples

``` r

  # Estimate the correlation between two variables in the debate dataset
  data(debate, package = "BayesRTMB")

  fit_corr <- rtmb_corr(cbind(sat, perf), data = debate)
#> Pre-checking model code...
#> Checking RTMB setup...

  # \donttest{
  mcmc_corr <- fit_corr$sample(sampling = 500, warmup = 500, chains = 2)
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
  mcmc_corr$summary()
#>   variable      mean    sd       map      q2.5     q97.5  ess_bulk  ess_tail  rhat 
#> lp          -1006.54  1.67  -1005.57  -1010.64  -1004.37       492       732  1.00 
#> corr[rho]       0.29  0.06      0.28      0.19      0.40       926       643  1.00 
#> mean[sat]       3.44  0.06      3.45      3.31      3.55       795       711  1.00 
#> mean[perf]      4.69  0.10      4.69      4.50      4.89       741       600  1.00 
#> sd[sat]         1.00  0.04      1.00      0.93      1.09      1126       830  1.01 
#> sd[perf]        1.77  0.08      1.78      1.62      1.92       970       778  1.00 

  bf_corr <- mcmc_corr$bayes_factor(fixed = list(corr = 0))
#> Calculating marginal likelihood for the full model...
#> Bridge Sampling Converged: LogML = -1013.161 (Error = 0.0104, ESS = 352.1)
#> 
#> --- Sampling from the comparison model ---
#> Starting sequential sampling (chains = 2)...
#> chain 1 started...
#> chain 1: iter 200/1500 (13%) warmup
#> chain 1: iter 400/1500 (26%) warmup
#> chain 1: iter 600/1500 (40%) warmup
#> chain 1: iter 800/1500 (53%) warmup
#> chain 1: iter 1000/1500 (66%) warmup
#> chain 1: iter 1200/1500 (80%) sampling
#> chain 1: iter 1400/1500 (93%) sampling
#> chain 1 done (100%)
#> chain 2 started...
#> chain 2: iter 200/1500 (13%) warmup
#> chain 2: iter 400/1500 (26%) warmup
#> chain 2: iter 600/1500 (40%) warmup
#> chain 2: iter 800/1500 (53%) warmup
#> chain 2: iter 1000/1500 (66%) warmup
#> chain 2: iter 1200/1500 (80%) sampling
#> chain 2: iter 1400/1500 (93%) sampling
#> chain 2 done (100%)
#> sampling: 100%
#> 
#> --- Calculating marginal likelihood for the comparison model ---
#> Bridge Sampling Converged: LogML = -1025.622 (Error = 0.0061, ESS = 285.4)
  print(bf_corr)
#> --- Bayes Factor Analysis (Bridge Sampling) ---
#> Bayes Factor (BF12) : 257903.3 
#> Log Bayes Factor    : 12.4603 (Approx. Error = 0.0121)
#> Evidence            : Decisive evidence for Model 1 
#> Comparison model    : Parameters fixed at list(corr = 0) 
#> 
  # }
```
