# RTMB-based Mediation Analysis Wrapper

\`rtmb_mediation\` performs mediation analysis by simultaneously
estimating multiple GLM regression equations. It automatically
identifies mediation paths and calculates indirect, direct, and total
effects.

## Usage

``` r
rtmb_mediation(
  formula,
  data,
  family = "gaussian",
  prior = prior_flat(),
  y_range = NULL,
  fixed = NULL,
  view = NULL,
  ...
)
```

## Arguments

- formula:

  A list of formulas defining the regression paths (e.g., \`list(M ~ X,
  Y ~ X + M)\`).

- data:

  A data frame containing the variables.

- family:

  A single character string or a list of character strings specifying
  the error distribution for each equation (e.g., \`family =
  list("gaussian", "binomial")\`). Default is "gaussian".

- prior:

  An object of class "rtmb_prior" specifying the prior distribution.

- y_range:

  Theoretical minimum and maximum values of the response variable.

- fixed:

  A named list of parameter values to fix (optional).

- view:

  Character vector of parameter names to prioritize in summary.

- ...:

  Additional arguments passed to the model construction.

## Value

An \`RTMB_Model\` object.

## Details

The function identifies mediation paths by looking for variables that
are responses in one equation and predictors in another. Indirect
effects are calculated as the product of coefficients along these paths
(\\a \* b\\).

**Uncertainty Estimation**: When using \`\$optimize(ci_method =
"sampling")\`, the function provides asymmetric confidence intervals for
indirect effects based on the distribution of products, which is more
accurate than the standard Sobel test (Delta Method).

## Examples

``` r

  # Simulate mediation data
  set.seed(123)
  N <- 100
  X <- rnorm(N)
  # M is influenced by X
  M <- 0.5 * X + rnorm(N, 0, 0.5)
  # Y is influenced by both X and M
  Y <- 0.3 * X + 0.8 * M + rnorm(N, 0, 0.5)
  dat <- data.frame(X, M, Y)

  # Fit a mediation model
  # The formula list specifies the two regression equations
  fit_med <- rtmb_mediation(list(M ~ X, Y ~ X + M), data = dat)
#> Pre-checking model code...
#> Checking RTMB setup...
  
  # Maximum A Posteriori (MAP) estimation
  map_med <- fit_med$optimize()
#> Starting RTMB optimization...
#> 
  # The summary automatically calculates indirect, direct, and total effects
  map_med$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 134.66
#> Approx. Log Marginal Likelihood (Laplace): -149.35
#> 
#> Point Estimates and 95% Wald CI:
#>      variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b1[Intercept]  -0.05140     0.04829   -0.14604    0.04324 
#> b1[X]           0.47376     0.05290    0.37008    0.57745 
#> b2[Intercept]   0.06753     0.04734   -0.02526    0.16032 
#> b2[X]           0.22151     0.06924    0.08580    0.35721 
#> b2[M]           0.82381     0.09750    0.63272    1.01490 
#> IE_X_M_Y        0.39029     0.06351    0.26582    0.51476 
#> DE_X_Y          0.22151     0.06924    0.08580    0.35721 
#> TE_X_M_Y        0.61180     0.06753    0.47945    0.74415 
#> sigma1          0.48048     0.03398    0.41830    0.55191 
#> sigma2          0.46846     0.03313    0.40783    0.53810 
#> 
  
  # MCMC sampling for more reliable confidence intervals of the indirect effect
  # (chains and iterations reduced for faster execution)
  # \donttest{
  mcmc_med <- fit_med$sample(sampling = 500, warmup = 500, chains = 2)
#> Starting sequential sampling (chains = 2)...
#> chain 1 started... 
#> chain 1: iter 100 warmup 
#> chain 1: iter 200 warmup 
#> chain 1: iter 300 warmup 
#> chain 1: iter 400 warmup 
#> chain 1: iter 500 warmup 
#> chain 1: iter 600 sampling 
#> chain 1: iter 700 sampling 
#> chain 1: iter 800 sampling 
#> chain 1: iter 900 sampling 
#> chain 1: iter 1000 sampling 
#> chain 2 started... 
#> chain 2: iter 100 warmup 
#> chain 2: iter 200 warmup 
#> chain 2: iter 300 warmup 
#> chain 2: iter 400 warmup 
#> chain 2: iter 500 warmup 
#> chain 2: iter 600 sampling 
#> chain 2: iter 700 sampling 
#> chain 2: iter 800 sampling 
#> chain 2: iter 900 sampling 
#> chain 2: iter 1000 sampling 
#> Calculating transformed parameters...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================| 100%
  mcmc_med$summary()
#>      variable     mean    sd      map     q2.5    q97.5  ess_bulk  ess_tail  rhat 
#> lp             -139.63  1.88  -138.55  -144.17  -136.94       400       571  1.01 
#> b1[Intercept]    -0.05  0.05    -0.03    -0.15     0.05      1686       869  1.00 
#> b1[X]             0.47  0.05     0.46     0.36     0.58      1462       731  1.00 
#> b2[Intercept]     0.06  0.05     0.06    -0.02     0.16      1274       682  1.00 
#> b2[X]             0.22  0.07     0.22     0.07     0.36      1042       753  1.00 
#> b2[M]             0.83  0.10     0.84     0.63     1.02       928       840  1.00 
#> IE_X_M_Y          0.39  0.07     0.38     0.27     0.54      1022       635  1.00 
#> DE_X_Y            0.22  0.07     0.22     0.07     0.36      1042       753  1.00 
#> TE_X_M_Y          0.61  0.07     0.62     0.47     0.75      1361       760  1.00 
#> sigma1            0.49  0.04     0.48     0.43     0.57      1349       780  1.01 
  # }
```
