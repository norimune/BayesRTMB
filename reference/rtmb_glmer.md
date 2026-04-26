# RTMB-based GLMM wrapper function

RTMB-based GLMM wrapper function

## Usage

``` r
rtmb_glmer(
  formula,
  data,
  family = "gaussian",
  laplace = FALSE,
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

  lme4-style formula (e.g., Y ~ X + (1 \| GID))

- data:

  Data frame

- family:

  Character string of the distribution family (e.g., "gaussian",
  "binomial", "poisson")

- laplace:

  Logical; whether to marginalize random effects using Laplace
  approximation

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

  List of initial values (generated automatically based on glm if
  omitted)

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
#> lp           -3138.91  5.92  -3137.07  -3151.78  -3129.03       250       548  1.00 
#> Intercept       21.08  2.80     21.53     15.61     26.54       608       781  1.00 
#> b[Time]          8.74  0.15      8.75      8.44      9.02      1267       685  1.00 
#> sigma           25.78  0.71     25.58     24.41     27.20      1076       793  1.00 
#> sd[Int]         15.10  1.08     15.44     13.15     17.36       505       710  1.01 
#> Intercept_c    114.71  2.31    114.86    110.11    119.19       469       693  1.00 
#> r_re[1]          0.22  0.79      0.15     -1.29      1.76      1730       345  1.00 
#> r_re[2]         -1.10  0.54     -1.10     -2.16     -0.09      1406       686  1.00 
#> r_re[3]         -1.09  0.56     -1.12     -2.23     -0.02      1243       636  1.00 
#> r_re[4]         -2.59  0.49     -2.54     -3.55     -1.63      1034       381  1.00 

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
#> Approx. Log Marginal Likelihood (Laplace): -125.78
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> b[cyl]    -0.29030     0.75457   -1.76922    1.18862 
#> b[disp]    0.00404     0.01351   -0.02244    0.03052 
#> b[hp]     -0.01784     0.01718   -0.05151    0.01583 
#> b[drat]    0.82133     1.12171   -1.37718    3.01985 
#> b[wt]     -2.60905     1.32613   -5.20822   -0.00989 
#> b[qsec]    0.43880     0.55250   -0.64408    1.52168 
#> b[vs]      0.19646     1.32807   -2.40651    2.79942 
#> b[am]      1.74820     1.34113   -0.88036    4.37676 
#> b[gear]    0.83000     1.02695   -1.18278    2.84278 
#> b[carb]   -0.51571     0.59074   -1.67354    0.64212 
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
#> b[drat]   -0.00000     0.00000   -0.00000   -0.00000 
#> b[wt]     -5.28853     0.00000   -5.28853   -5.28853 
#> b[qsec]    0.00000     0.00000    0.00000    0.00000 
#> b[vs]      0.00000     0.00000    0.00000    0.00000 
#> b[am]      0.00000     0.00000    0.00000    0.00000 
#> b[gear]    0.00000     0.00000    0.00000    0.00000 
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
#> b[cyl]    -0.46  0.66  -0.03  -2.02   0.35       493       662  1.00 
#> b[disp]   -0.00  0.01  -0.00  -0.02   0.01       419       376  1.00 
#> b[hp]     -0.01  0.01  -0.00  -0.05   0.01       451       842  1.00 
#> b[drat]    0.38  0.91   0.01  -0.94   2.82       674       862  1.00 
#> b[wt]     -2.78  1.32  -3.04  -5.07   0.01       387       335  1.01 
#> b[qsec]    0.25  0.44   0.00  -0.33   1.33       605       861  1.00 
#> b[vs]      0.32  1.01   0.00  -1.42   3.01       737       729  1.00 
#> b[am]      1.19  1.62   0.04  -0.76   5.01       505       904  1.00 
#> b[gear]    0.21  0.71  -0.01  -1.02   2.26       781       782  1.01 
#> b[carb]   -0.31  0.45  -0.01  -1.42   0.21       507       893  1.01 
# }
```
