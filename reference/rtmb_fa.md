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
  missing = c("listwise", "fiml")
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

## Examples

``` r

  # Prepare a subset of the mtcars dataset for factor analysis
  # Scaling is recommended for variables with different units
  fa_data <- scale(mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")])

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
#> Negative Log-Posterior: 191.21
#> Approx. Log Marginal Likelihood (Laplace): -216.00
#> 
#> Point Estimates and 95% Wald CI:
#>        variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L[mpg,Factor1]    0.91547     0.03740    0.84216    0.98877 
#> L[disp,Factor1]  -0.95435     0.02609   -1.00548   -0.90322 
#> L[hp,Factor1]    -0.80119     0.06853   -0.93552   -0.66687 
#> L[drat,Factor1]   0.73233     0.08608    0.56363    0.90103 
#> L[wt,Factor1]    -0.92203     0.03253   -0.98578   -0.85828 
#> L[qsec,Factor1]   0.40864     0.15340    0.10798    0.70930 
#> sd[mpg]           0.39606     0.07319    0.27572    0.56893 
#> sd[disp]          0.29398     0.07497    0.17835    0.48460 
#> sd[hp]            0.58899     0.08135    0.44930    0.77210 
#> sd[drat]          0.67023     0.08827    0.51774    0.86762 
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
#> Negative Log-Posterior: 169.75
#> Approx. Log Marginal Likelihood (Laplace): -202.89
#> 
#> Point Estimates and 95% Sampling-based CI:
#>               variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L_promax[mpg,Factor1]   -0.84391     0.09234   -0.94105   -0.58754 
#> L_promax[disp,Factor1]   0.86488     0.09796    0.58671    0.95548 
#> L_promax[hp,Factor1]     0.48739     0.17207    0.06884    0.72643 
#> L_promax[drat,Factor1]  -0.81857     0.13502   -0.94834   -0.53021 
#> L_promax[wt,Factor1]     1.05698     0.10321    0.86654    1.11629 
#> L_promax[qsec,Factor1]   0.17103     0.11581    0.01298    0.26388 
#> L_promax[mpg,Factor2]    0.14304     0.12712   -0.02929    0.40829 
#> L_promax[disp,Factor2]  -0.15096     0.13212   -0.43374    0.01984 
#> L_promax[hp,Factor2]    -0.57624     0.18142   -0.91341   -0.26952 
#> L_promax[drat,Factor2]  -0.18382     0.13165   -0.38770    0.11398 
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
#> Negative Log-Posterior: 169.75
#> Approx. Log Marginal Likelihood (Laplace): -202.89
#> 
#> Point Estimates and 95% Sampling-based CI:
#>       variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L_varimax[1,1]  -0.84186     0.05720   -0.92282   -0.70183 
#> L_varimax[2,1]   0.86374     0.05952    0.70611    0.94049 
#> L_varimax[3,1]   0.59564     0.10327    0.37542    0.77071 
#> L_varimax[4,1]  -0.74507     0.10571   -0.86734   -0.51778 
#> L_varimax[5,1]   0.97372     0.06910    0.80969    0.99466 
#> L_varimax[6,1]  -0.06546     0.11314   -0.26331    0.12112 
#> L_varimax[1,2]   0.37581     0.13654    0.13184    0.56858 
#> L_varimax[2,2]  -0.38940     0.14306   -0.59211   -0.15393 
#> L_varimax[3,2]  -0.69839     0.14894   -0.86577   -0.39394 
#> L_varimax[4,2]   0.04998     0.13785   -0.24406    0.27966 
#> 



  # --- 3. Regularized Factor Analysis using Spike-and-Slab Prior (SSP) ---
  # Specifying 'rotate = "ssp"' enables sparse loading matrix estimation
  fit_ssp <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "ssp")
#> Pre-checking model code...
#> Checking RTMB setup...

  map_ssp <- fit_ssp$optimize()
#> Starting RTMB optimization...
#> 
#> Warning: Optimization did not converge ( convergence code = 1; message = singular convergence (7)). Estimates may be unreliable; consider increasing num_estimate, changing initial values, or adjusting optimizer control settings.
#> SE warning: sdreport() produced non-finite standard errors; Hessian-based fallback will be attempted.
#> SE warning: Hessian matrix was singular; using MASS::ginv() to approximate the covariance matrix.
  map_ssp$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 175.35
#> Approx. Log Marginal Likelihood (Laplace): NA
#> 
#> Point Estimates and 95% Wald CI:
#>        variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L[mpg,Factor1]    0.83826     0.05421    0.73201    0.94450 
#> L[disp,Factor1]  -0.86320     0.04963   -0.96047   -0.76594 
#> L[hp,Factor1]    -0.56058     0.09806   -0.75277   -0.36839 
#> L[drat,Factor1]   0.71346     0.09036    0.53636    0.89057 
#> L[wt,Factor1]    -0.97993     0.02361   -1.02620   -0.93365 
#> L[qsec,Factor1]   0.00000     0.00000   -0.00000    0.00000 
#> L[mpg,Factor2]    0.28214     0.09339    0.09909    0.46519 
#> L[disp,Factor2]  -0.29761     0.08450   -0.46323   -0.13198 
#> L[hp,Factor2]    -0.66648     0.11892   -0.89956   -0.43339 
#> L[drat,Factor2]   0.00000     0.00000    0.00000    0.00000 
#> 

  # MCMC sampling for the SSP model (chains and iterations reduced for faster execution)
  if (FALSE) { # \dontrun{
  mcmc_ssp <- fit_ssp$sample(sampling = 500, warmup = 500, chains = 2)

  # --- 4. Resolving Label Switching in MCMC ---
  # In Bayesian factor analysis, MCMC chains often suffer from label switching or sign flipping.
  # This can be resolved by applying post-hoc Procrustes rotation to the posterior samples.

  # Summary of unrotated loadings (may show poor convergence / large SE due to switching)
  mcmc_ssp$summary("L")
  mcmc_ssp$draws("L[mpg,1]") |> plot_dens()

  # Apply Procrustes rotation targeting the loading matrix "L"
  mcmc_ssp$rotate(target = "L")

  # Summary of the rotated loadings (L_rot) with stabilized estimates
  mcmc_ssp$summary("L_rot")
  mcmc_ssp$draws("L_rot[mpg,1]") |> plot_dens()
  mcmc_ssp$draws("L_rot[mpg,1]") |> plot_dens()
  } # }
```
