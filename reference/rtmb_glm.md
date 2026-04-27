# RTMB-based GLM wrapper function (no random effects)

RTMB-based GLM wrapper function (no random effects)

## Usage

``` r
rtmb_glm(
  formula,
  data,
  family = "gaussian",
  penalty = c("none", "rhs", "ssp"),
  y_range = NULL,
  use_weak_info = FALSE,
  prior = list(Intercept_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, nu_rate = 0.1,
    cutpoint_sd = 2.5, shape_rate = 1, phi_rate = 1, lkj_eta = 1),
  weak_info_prior = list(max_beta = 1, sd_ratio = 0.5, expected_vars = 3, slab_scale = 2,
    slab_df = 4, ssp_ratio = 0.25),
  init = NULL,
  null = NULL
)
```

## Arguments

- formula:

  Formula

- data:

  Data frame

- family:

  Character string of the distribution family (e.g., "gaussian",
  "binomial", "poisson")

- penalty:

  Type of regularization for fixed effects: "none", "rhs" (Regularized
  Horseshoe), or "ssp" (Spike and Slab Prior). Default is "none".

- y_range:

  Theoretical minimum and maximum values of the response variable as a
  vector c(min, max). Specifying this automatically enables weakly
  informative priors.

- use_weak_info:

  Logical; whether to explicitly use weakly informative priors (requires
  y_range for continuous models).

- prior:

  List of hyperparameters for the default fixed priors.

- weak_info_prior:

  List of hyperparameters for the weakly informative priors and
  regularization.

- init:

  List of initial values

- null:

  Character string specifying the target parameter for the null model.

## Examples

``` r
# \donttest{
  # --- 1. Linear Regression (rtmb_lm) ---
  # Fit a linear regression model using the mtcars dataset
  fit_lm <- rtmb_lm(mpg ~ wt + cyl, data = mtcars)
#> Pre-checking model code...
#> Checking RTMB setup...
  map_lm <- fit_lm$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 95.45
  map_lm$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 95.45
#> Approx. Log Marginal Likelihood (Laplace): -96.31
#> 
#> Point Estimates and 95% Wald CI:
#>    variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept    39.64829     1.41486   36.87522   42.42136 
#> b[wt]        -3.18109     0.62316   -4.40247   -1.95971 
#> b[cyl]       -1.51135     0.34159   -2.18086   -0.84183 
#> sigma         2.11659     0.21630    1.73241    2.58596 
#> Intercept_c  20.06248     0.37432   19.32883   20.79613 
#> 

  # --- 2. Generalized Linear Model (rtmb_glm) ---
  # Fit a logistic regression model
  fit_glm <- rtmb_glm(am ~ mpg + hp, data = mtcars, family = "binomial")
#> Pre-checking model code...
#> Checking RTMB setup...
  map_glm <- fit_glm$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 19.29
  map_glm$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 19.29
#> Approx. Log Marginal Likelihood (Laplace): -22.33
#> 
#> Point Estimates and 95% Wald CI:
#>    variable   Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept    -33.35406    14.99362  -62.74101   -3.96710 
#> b[mpg]         1.25569     0.56433    0.14962    2.36175 
#> b[hp]          0.05387     0.02678    0.00138    0.10636 
#> Intercept_c   -0.22426     0.57097   -1.34333    0.89482 
#> 

  # --- 3. Generalized Linear Mixed Model (rtmb_glmer) ---
  # Fit a linear mixed-effects model with random intercepts
  fit_glmer <- rtmb_glmer(weight ~ Time + (1 | Chick), data = ChickWeight, family = "gaussian")
#> Pre-checking model code...
#> Checking RTMB setup...

  # MAP estimation using Laplace approximation for random effects
  map_glmer <- fit_glmer$optimize(laplace = TRUE)
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 3193.13
  map_glmer$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 3193.13
#> Approx. Log Marginal Likelihood (Laplace): -3193.34
#> Note: Random effects are stored in $random_effects
#> 
#> Point Estimates and 95% Wald CI:
#>    variable   Estimate  Std. Error  Lower 95%  Upper 95% 
#> Intercept     25.09932     2.66286   19.88022   30.31843 
#> b[Time]        8.79909     0.21060    8.38633    9.21185 
#> sigma         34.16198     0.83919   32.55618   35.84699 
#> sd[Int]        0.03258     0.08259    0.00023    4.68872 
#> Intercept_c  119.40786     1.41298  116.63847  122.17725 
#> 

  # MCMC sampling (chains and iterations reduced for faster execution)
  mcmc_glmer <- fit_glmer$sample(sampling = 500, warmup = 500, chains = 2)
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
  mcmc_glmer$summary()
#>    variable      mean    sd       map      q2.5     q97.5  ess_bulk  ess_tail  rhat 
#> lp           -3138.51  5.74  -3138.33  -3150.94  -3128.10       305       460  1.00 
#> Intercept       21.49  2.86     21.53     15.72     26.98       593       708  1.00 
#> b[Time]          8.74  0.16      8.74      8.43      9.05      1387       577  1.00 
#> sigma           25.75  0.71     25.73     24.36     27.15      1058       892  1.00 
#> sd[Int]         15.12  1.07     15.14     13.04     17.32       538       476  1.00 
#> Intercept_c    115.12  2.29    115.72    110.38    119.46       494       639  1.01 
#> r_re[1]          0.17  0.78      0.20     -1.28      1.69      1314       752  1.00 
#> r_re[2]         -1.11  0.55     -1.10     -2.23     -0.03      1236       841  1.00 
#> r_re[3]         -1.10  0.53     -0.91     -2.20     -0.14       995       714  1.00 
#> r_re[4]         -2.61  0.50     -2.52     -3.62     -1.67      1124       810  1.00 

  # --- 4. Regularized Regression (Variable Selection) ---
  # You can apply regularization to the fixed effects to shrink noise variables towards zero.
  # Use penalty = "rhs" for the Regularized Horseshoe prior, or "ssp" for the Spike-and-Slab prior.
  # Note: When using regularization, you must specify 'y_range' (the theoretical minimum and maximum
  # values of the response variable) to automatically set up the required weakly informative priors.

  # Fit a linear regression using all predictors in mtcars with the Horseshoe prior
  # 'mpg' theoretically ranges roughly between 0 and 40
  fit_rhs <- rtmb_lm(mpg ~ ., data = mtcars, penalty = "rhs", y_range = c(0, 40))
#> Pre-checking model code...
#> Checking RTMB setup...
  map_rhs <- fit_rhs$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 101.68
  # Summarize only the fixed effects (slopes)
  map_rhs$summary("b")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 101.68
#> Approx. Log Marginal Likelihood (Laplace): -126.01
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[cyl]    -0.29045     0.75452   -1.76928    1.18839 
#> b[disp]    0.00404     0.01353   -0.02248    0.03057 
#> b[hp]     -0.01783     0.01719   -0.05152    0.01585 
#> b[drat]    0.82136     1.12172   -1.37717    3.01988 
#> b[wt]     -2.60936     1.32742   -5.21105   -0.00766 
#> b[qsec]    0.43886     0.55258   -0.64417    1.52189 
#> b[vs]      0.19642     1.32811   -2.40663    2.79948 
#> b[am]      1.74822     1.34119   -0.88046    4.37691 
#> b[gear]    0.82990     1.02696   -1.18289    2.84270 
#> b[carb]   -0.51563     0.59116   -1.67429    0.64303 
#> 

  # Fit a linear regression with the Spike-and-Slab prior
  fit_ssp <- rtmb_lm(mpg ~ ., data = mtcars, penalty = "ssp", y_range = c(0, 40))
#> Pre-checking model code...
#> Checking RTMB setup...
  map_ssp <- fit_ssp$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 79.81
  map_ssp$summary("b")
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 79.81
#> Approx. Log Marginal Likelihood (Laplace): NA
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[cyl]     0.00000     0.00000    0.00000    0.00000 
#> b[disp]    0.00000     0.00000    0.00000    0.00000 
#> b[hp]      0.00000     0.00000    0.00000    0.00000 
#> b[drat]    0.00000     0.00000    0.00000    0.00000 
#> b[wt]     -5.28858     0.00000   -5.28858   -5.28858 
#> b[qsec]    0.00000     0.00000    0.00000    0.00000 
#> b[vs]      0.00000     0.00000    0.00000    0.00000 
#> b[am]      0.00000     0.00000    0.00000    0.00000 
#> b[gear]   -0.00000     0.00000   -0.00000   -0.00000 
#> b[carb]    0.00000     0.00000    0.00000    0.00000 
#> 

  # For models with complex penalties, MCMC often provides more reliable credible intervals
  mcmc_ssp <- fit_ssp$sample(sampling = 500, warmup = 500, chains = 2)
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
  mcmc_ssp$summary("b")
#> variable   mean    sd    map   q2.5  q97.5  ess_bulk  ess_tail  rhat 
#> b[cyl]    -0.39  0.60  -0.02  -1.89   0.36       539       841  1.00 
#> b[disp]   -0.00  0.01   0.00  -0.02   0.01       770       748  1.00 
#> b[hp]     -0.01  0.01  -0.00  -0.05   0.01       514       753  1.00 
#> b[drat]    0.40  0.94   0.01  -0.99   3.01       787       843  1.00 
#> b[wt]     -2.90  1.33  -3.39  -5.30   0.01       547       369  1.00 
#> b[qsec]    0.26  0.42   0.00  -0.26   1.37       466       744  1.00 
#> b[vs]      0.33  1.04   0.01  -1.32   3.04       709       904  1.00 
#> b[am]      1.09  1.62   0.01  -0.92   5.10       731       864  1.00 
#> b[gear]    0.27  0.77   0.01  -0.84   2.27       750       761  1.00 
#> b[carb]   -0.29  0.48  -0.00  -1.63   0.26       627       692  1.00 
# }
```
