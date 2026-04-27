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
  prior = list(mean_sd = 10, loadings_sd = 1, sd_rate = 10, ssp_ratio = 0.25),
  init = NULL
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

  List of hyperparameters for prior distributions. \`ssp_ratio\`
  represents the proportion of non-zero loadings per factor when "ssp"
  is specified.

- init:

  List of initial values. If not provided, initial values are
  automatically generated based on PCA or the psych package.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Prepare a subset of the mtcars dataset for factor analysis
  # Scaling is recommended for variables with different units
  fa_data <- scale(mtcars[, c("mpg", "disp", "hp", "drat", "wt", "qsec")])

  # --- 1. Standard Exploratory Factor Analysis (1 Factor) ---
  fit_fa1 <- rtmb_fa(data = fa_data, nfactors = 1)

  # Maximum A Posteriori (MAP) estimation
  map_fa1 <- fit_fa1$optimize()
  map_fa1$summary()



  # --- 2. Factor Analysis with Rotation and Factor Scores (2 Factors) ---
  # Extract 2 factors, apply Promax rotation during model fitting, and calculate factor scores
  fit_fa2 <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "promax", score = TRUE)

  # MCMC sampling for the 2 factor-model (chains and iterations reduced for faster execution)
  mcmc_fa2 <- fit_fa2$sample(sampling=500, warmup=500, chains=2)
  # The summary prioritizes rotated loadings (L_promax), standard deviations,
  # and factor correlations
  mcmc_fa2$summary()

  # Setting 'se_sampling = TRUE' enables the calculation of standard errors and 95% CIs
  # for transformed and generated quantities, such as factor scores and post-hoc rotations.
  # (It uses multivariate normal sampling from the unconstrained parameter space)
  map_fa2 <- fit_fa2$optimize(se_sampling = TRUE)
  map_fa2$summary()

  # Post-hoc rotation using the fa_rotate() method
  # You can also apply other rotation methods (e.g., "varimax" from the stats package)
  # to the unrotated loading matrix ("L") after estimation.
  map_fa2$fa_rotate(target = "L", rotate = "varimax")

  # The post-hoc rotated loadings are automatically stored with the method's
  # suffix (e.g., "L_varimax")
  map_fa2$summary("L_varimax")



  # --- 3. Regularized Factor Analysis using Spike-and-Slab Prior (SSP) ---
  # Specifying 'rotate = "ssp"' enables sparse loading matrix estimation
  fit_ssp <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "ssp")

  map_ssp <- fit_ssp$optimize()
  map_ssp$summary()

  # MCMC sampling for the SSP model (chains and iterations reduced for faster execution)
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
} # }
```
