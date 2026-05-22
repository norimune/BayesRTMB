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
  \`prior_weak()\`. Hyperparameters can be specified within these
  functions (e.g., \`prior_normal(mean_sd = 10, sd_rate = 10)\`).
  Available parameters for FA: \`mean_sd\`, \`sd_rate\`,
  \`loadings_sd\`, and \`ssp_ratio\` (if \`rotate = "ssp"\`).

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
#> 
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



  # --- 2. Factor Analysis with Rotation and Factor Scores (2 Factors) ---
  # Extract 2 factors, apply Promax rotation during model fitting, and calculate factor scores
  fit_fa2 <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "promax", score = TRUE)
#> Pre-checking model code...
#> Checking RTMB setup...

  # MCMC sampling for the 2 factor-model (chains and iterations reduced for faster execution)
  if (FALSE) { # \dontrun{
  mcmc_fa2 <- fit_fa2$sample(sampling=500, warmup=500, chains=2)
  # The summary prioritizes rotated loadings (L_promax), standard deviations,
  # and factor correlations
  mcmc_fa2$summary()
  } # }

  # Setting se_method = "sampling" enables standard errors and 95% CIs
  # for transformed and generated quantities, such as factor scores and post-hoc rotations.
  # (It uses multivariate normal sampling from the unconstrained parameter space)
  map_fa2 <- fit_fa2$optimize(se_method = "sampling")
#> Starting RTMB optimization...
#> 
#> Using simulation-based error propagation (1000 samples)...
  map_fa2$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 2440.43
#> Approx. Log Marginal Likelihood (Laplace): -2503.48
#> 
#> Point Estimates and 95% Sampling-based CI:
#>               variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L_promax[BF1,Factor1]   -0.88412     0.19040   -0.99590   -0.29552 
#> L_promax[BF2,Factor1]   -0.07323     0.04808   -0.16583    0.02137 
#> L_promax[BF3,Factor1]   -0.11607     0.07986   -0.27315    0.04090 
#> L_promax[BF4,Factor1]   -0.07882     0.08795   -0.24846    0.09597 
#> L_promax[BF5,Factor1]    0.03564     0.08367   -0.12772    0.19089 
#> L_promax[BF6,Factor1]    0.68867     0.12762    0.39206    0.88820 
#> L_promax[BF7,Factor1]    0.01309     0.04747   -0.07218    0.09987 
#> L_promax[BF8,Factor1]   -0.15549     0.09734   -0.35668    0.02904 
#> L_promax[BF9,Factor1]    0.08481     0.09385   -0.10586    0.26402 
#> L_promax[BF10,Factor1]  -0.10916     0.07976   -0.26512    0.04661 
#> 

  # Post-hoc rotation using the fa_rotate() method
  # You can also apply other rotation methods (e.g., "varimax" from the stats package)
  # to the unrotated loading matrix ("L") after estimation.
  map_fa2$fa_rotate(target = "L", rotate = "varimax")
#> Applying varimax rotation to L (Saving to generate as _varimax)...
#> Generated quantities updated.

  # The post-hoc rotated loadings are automatically stored with the method's
  # suffix (e.g., "L_varimax")
  map_fa2$summary("L_varimax")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 2440.43
#> Approx. Log Marginal Likelihood (Laplace): -2503.48
#> 
#> Point Estimates and 95% Sampling-based CI:
#>        variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L_varimax[1,1]   -0.87739     0.18946   -0.98930   -0.27963 
#> L_varimax[2,1]   -0.00088     0.13712   -0.34153    0.19608 
#> L_varimax[3,1]   -0.12824     0.08587   -0.27944    0.04857 
#> L_varimax[4,1]   -0.09712     0.09498   -0.27465    0.09389 
#> L_varimax[5,1]    0.02261     0.08562   -0.13019    0.20885 
#> L_varimax[6,1]    0.68707     0.12478    0.38952    0.87753 
#> L_varimax[7,1]    0.08960     0.14207   -0.28461    0.28011 
#> L_varimax[8,1]   -0.17801     0.10885   -0.38106    0.03885 
#> L_varimax[9,1]    0.07460     0.09842   -0.11166    0.27545 
#> L_varimax[10,1]  -0.08975     0.08971   -0.28869    0.05730 
#> 



  # --- 3. Regularized Factor Analysis using Spike-and-Slab Prior (SSP) ---
  # Specifying 'rotate = "ssp"' enables sparse loading matrix estimation
  fit_ssp <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "ssp")
#> Pre-checking model code...
#> Checking RTMB setup...

  map_ssp <- fit_ssp$optimize()
#> Starting RTMB optimization...
#> 
#> Warning: Optimization did not converge ( convergence code = 1; message = function evaluation limit reached without convergence (9)). Estimates may be unreliable; consider increasing num_estimate, changing initial values, or adjusting optimizer control settings.
  map_ssp$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 2433.06
#> Approx. Log Marginal Likelihood (Laplace): -2692.30
#> 
#> Point Estimates and 95% Wald CI:
#>        variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L[BF1,Factor1]    0.00000     0.00000   -0.00000    0.00000 
#> L[BF2,Factor1]    0.78286     0.07926    0.62751    0.93821 
#> L[BF3,Factor1]   -0.11537     0.09039   -0.29252    0.06179 
#> L[BF4,Factor1]   -0.18356     0.08405   -0.34831   -0.01882 
#> L[BF5,Factor1]   -0.11921     0.08829   -0.29226    0.05383 
#> L[BF6,Factor1]    0.00000     0.00000   -0.00000    0.00000 
#> L[BF7,Factor1]    0.85423     0.08408    0.68944    1.01902 
#> L[BF8,Factor1]   -0.23762     0.08149   -0.39735   -0.07790 
#> L[BF9,Factor1]   -0.07829     0.09694   -0.26828    0.11170 
#> L[BF10,Factor1]   0.18712     0.09227    0.00628    0.36797 
#> 

  # MCMC sampling for the SSP model (chains and iterations reduced for faster execution)
  if (FALSE) { # \dontrun{
  mcmc_ssp <- fit_ssp$sample(sampling = 500, warmup = 500, chains = 2)

  # --- 4. Resolving Label Switching in MCMC ---
  # In Bayesian factor analysis, MCMC chains often suffer from label switching or sign flipping.
  # This can be resolved by applying post-hoc Procrustes rotation to the posterior samples.

  # Summary of unrotated loadings (may show poor convergence / large SE due to switching)
  mcmc_ssp$summary("L")
  mcmc_ssp$draws("L[BF1,1]") |> plot_dens()

  # Apply Procrustes rotation targeting the loading matrix "L"
  mcmc_ssp$rotate(target = "L")

  # Summary of the rotated loadings (L_rot) with stabilized estimates
  mcmc_ssp$summary("L_rot")
  mcmc_ssp$draws("L_rot[BF1,1]") |> plot_dens()
  } # }
```
