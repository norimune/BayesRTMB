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
  \`prior_normal(mean_sd = 10, sd_rate = 10)\`). Available parameters
  for FA: \`mean_sd\`, \`sd_rate\`, \`loadings_sd\`, and \`ssp_ratio\`
  (if \`rotate = "ssp"\`).

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
#> L_promax[BF1,Factor1]   -0.88412     0.19045   -0.99579   -0.29916 
#> L_promax[BF2,Factor1]   -0.07323     0.05383   -0.16711    0.01892 
#> L_promax[BF3,Factor1]   -0.11607     0.08270   -0.27486    0.04811 
#> L_promax[BF4,Factor1]   -0.07882     0.08821   -0.24792    0.09878 
#> L_promax[BF5,Factor1]    0.03564     0.08337   -0.12598    0.19492 
#> L_promax[BF6,Factor1]    0.68867     0.12798    0.38591    0.88749 
#> L_promax[BF7,Factor1]    0.01309     0.05165   -0.07151    0.10295 
#> L_promax[BF8,Factor1]   -0.15549     0.09684   -0.34706    0.02727 
#> L_promax[BF9,Factor1]    0.08481     0.09297   -0.09648    0.26530 
#> L_promax[BF10,Factor1]  -0.10916     0.07964   -0.26156    0.05252 
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
#> L_varimax[1,1]   -0.87739     0.18907   -0.98820   -0.29034 
#> L_varimax[2,1]   -0.00088     0.13695   -0.34629    0.19880 
#> L_varimax[3,1]   -0.12824     0.08883   -0.28335    0.04773 
#> L_varimax[4,1]   -0.09712     0.09463   -0.26598    0.10075 
#> L_varimax[5,1]    0.02261     0.08509   -0.13838    0.21208 
#> L_varimax[6,1]    0.68707     0.12595    0.38894    0.87870 
#> L_varimax[7,1]    0.08960     0.14184   -0.27815    0.29630 
#> L_varimax[8,1]   -0.17801     0.10725   -0.37010    0.04470 
#> L_varimax[9,1]    0.07460     0.09745   -0.09915    0.28085 
#> L_varimax[10,1]  -0.08975     0.08811   -0.28202    0.05605 
#> 



  # --- 3. Regularized Factor Analysis using Spike-and-Slab Prior (SSP) ---
  # Specifying 'rotate = "ssp"' enables sparse loading matrix estimation
  fit_ssp <- rtmb_fa(data = fa_data, nfactors = 2, rotate = "ssp")
#> Pre-checking model code...
#> Checking RTMB setup...

  map_ssp <- fit_ssp$optimize()
#> Starting RTMB optimization...
#> 
#> Warning: No optimization run reached an acceptable optimizer status; BEST is selected by the lowest objective among not converged runs ( convergence code = 1; message = function evaluation limit reached without convergence (9)). Estimates may be unreliable; consider increasing num_estimate, changing initial values, or adjusting optimizer control settings.
  map_ssp$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 2433.05
#> Approx. Log Marginal Likelihood (Laplace): -2716.78
#> 
#> Point Estimates and 95% Wald CI:
#>        variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> L[BF1,Factor1]    0.00000     0.00000   -0.00000    0.00000 
#> L[BF2,Factor1]    0.78201     0.08335    0.61865    0.94538 
#> L[BF3,Factor1]   -0.11589     0.09062   -0.29350    0.06172 
#> L[BF4,Factor1]   -0.18487     0.08384   -0.34920   -0.02054 
#> L[BF5,Factor1]   -0.11914     0.08861   -0.29282    0.05454 
#> L[BF6,Factor1]    0.00000     0.00000   -0.00000    0.00000 
#> L[BF7,Factor1]    0.85548     0.08895    0.68114    1.02982 
#> L[BF8,Factor1]   -0.23703     0.08306   -0.39982   -0.07424 
#> L[BF9,Factor1]   -0.07986     0.09650   -0.26900    0.10928 
#> L[BF10,Factor1]   0.18645     0.09386    0.00248    0.37042 
#> 
```
