# Wrapper for estimating correlation matrix (multivariate normal distribution)

Estimates a correlation matrix (along with means and standard
deviations) assuming a multivariate normal distribution from observation
data. If there are 2 observed variables, it automatically switches to
estimate the scalar correlation coefficient (\`corr\`) directly.

## Usage

``` r
rtmb_corr(
  data,
  prior = list(lkj_eta = 1, mu_sd = 10, sigma_rate = 1),
  init = NULL,
  null = NULL
)
```

## Arguments

- data:

  Observation data frame or matrix (N x P).

- prior:

  List of hyperparameters for prior distributions. Default is
  `list(lkj_eta = 1.0, mu_sd = 10, sigma_rate = 1.0)`.

- init:

  List of initial values (optional).

- null:

  Target when creating a null model (e.g., `"corr"`). Optional.

## Value

An instance of the `RTMB_Model` class.

## Examples

``` r
# \donttest{
  # Simulate bivariate normal data with a true correlation of 0.5
  set.seed(123)
  N <- 50
  rho <- 0.5
  cov_mat <- matrix(c(1, rho, rho, 1), nrow = 2)

  if (requireNamespace("MASS", quietly = TRUE)) {
    data_corr <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = cov_mat)
    colnames(data_corr) <- c("X1", "X2")

    fit_corr <- rtmb_corr(data = data_corr)

    mcmc_corr <- fit_corr$sample(sampling = 500, warmup = 500, chains = 2)
    mcmc_corr$summary()

    bf_corr <- mcmc_corr$bayes_factor(null_model = "corr")
    print(bf_corr)
  }
#> Pre-checking model code...
#> Checking RTMB setup...
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
#> Calculating marginal likelihood for the full model...
#> Bridge Sampling Converged: LogML = -140.625 (Error = 0.0144, ESS = 243.1)
#> 
#> --- Preparing and Sampling Null Model (corr) ---
#> Auto-completed target: corr ~ lkj_corr(prior_lkj_eta)
#> Null model created. Fixed parameters:
#>   corr -> 0 (Dimensions: 1)
#>   Prior correction applied: -1.386294
#> Starting sequential sampling (chains = 2)...
#> chain 1 started... 
#> chain 1: iter 100 warmup 
#> chain 1: iter 200 warmup 
#> chain 1: iter 300 warmup 
#> chain 1: iter 400 warmup 
#> chain 1: iter 500 warmup 
#> chain 1: iter 600 warmup 
#> chain 1: iter 700 warmup 
#> chain 1: iter 800 warmup 
#> chain 1: iter 900 warmup 
#> chain 1: iter 1000 warmup 
#> chain 1: iter 1100 sampling 
#> chain 1: iter 1200 sampling 
#> chain 1: iter 1300 sampling 
#> chain 1: iter 1400 sampling 
#> chain 1: iter 1500 sampling 
#> chain 2 started... 
#> chain 2: iter 100 warmup 
#> chain 2: iter 200 warmup 
#> chain 2: iter 300 warmup 
#> chain 2: iter 400 warmup 
#> chain 2: iter 500 warmup 
#> chain 2: iter 600 warmup 
#> chain 2: iter 700 warmup 
#> chain 2: iter 800 warmup 
#> chain 2: iter 900 warmup 
#> chain 2: iter 1000 warmup 
#> chain 2: iter 1100 sampling 
#> chain 2: iter 1200 sampling 
#> chain 2: iter 1300 sampling 
#> chain 2: iter 1400 sampling 
#> chain 2: iter 1500 sampling 
#> 
#> --- Calculating marginal likelihood for the null model ---
#> Bridge Sampling Converged: LogML = -146.075 (Error = 0.0107, ESS = 375.4)
#> Bayes Factor (BF12) : 232.7278 
#> Log Bayes Factor    : 5.4499 (Approx. Error = 0.0179)
#> Interpretation      : Decisive evidence for Model 1 
# }
```
