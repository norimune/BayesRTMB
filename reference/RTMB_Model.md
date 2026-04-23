# Wrapper function to create an RTMB_Model instance

Combines user-defined data and model code (likelihood and priors) to
create an \`RTMB_Model\` (R6 class) instance for Bayesian inference
(MCMC, VB, MAP estimation).

## Usage

``` r
rtmb_model(
  data,
  code,
  par_names = list(),
  init = NULL,
  view = NULL,
  null_target = NULL
)
```

## Arguments

- data:

  A named list containing observation data and constants (e.g., sample
  size) used in the model.

- code:

  A model definition block described by \`rtmb_code(...)\` (including
  data, parameters, model, transform, generate).

- par_names:

  A list of specific variable names corresponding to the dimensions of
  each parameter (optional).

- init:

  A list or numeric vector of initial values for parameters (optional).
  If not specified, initialized randomly.

- view:

  Character vector of parameter names to be displayed preferentially at
  the top when outputting results like \`summary()\` (optional).

- null_target:

  Character string. To simultaneously create a null model, specify the
  target parameter to fix and the prior distribution to disable as a
  formula string (e.g., `"delta ~ cauchy(0, r)"`). Default is `NULL`.

## Value

An \`RTMB_Model\` class instance with compiled and pre-tested model.

## Examples

``` r
# \donttest{
# Simulate data for 3 groups
set.seed(123)
N <- 60
group_idx <- sample(1:3, N, replace = TRUE)
group_names <- c("Control", "Treatment_A", "Treatment_B")

# True group means: Control = 0, Treatment_A = 2, Treatment_B = -1
true_means <- c(0, 2, -1)
y <- true_means[group_idx] + rnorm(N, mean = 0, sd = 0.5)

data_list <- list(N = N, K = 3, group_idx = group_idx, y = y)

# Define the model using rtmb_code
model_code <- rtmb_code(
  parameters = {
    mu    = Dim(K)             # Vector of length K (group means)
    sigma = Dim(1, lower = 0)  # Scalar (residual standard deviation)
  },
  model = {
    # Priors
    for (k in 1:K) mu[k] ~ normal(0, 10)
    sigma ~ exponential(1)

    # Likelihood
    for (i in 1:N) {
      y[i] ~ normal(mu[group_idx[i]], sigma)
    }
  }
)

# --- 1. Basic Model Creation ---
# Create the RTMB_Model object
mod_basic <- rtmb_model(
  data = data_list,
  code = model_code
)
#> Pre-checking model code...
#> Checking RTMB setup...

# Perform Maximum A Posteriori (MAP) estimation
map_basic <- mod_basic$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 46.99
# The summary displays default parameter names: mu[1], mu[2], mu[3]
map_basic$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 46.99
#> Approx. Log Marginal Likelihood (Laplace): -52.63
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> mu[1]     -0.12319     0.09499   -0.30936    0.06298 
#> mu[2]      2.15919     0.10806    1.94741    2.37098 
#> mu[3]     -1.11887     0.09722   -1.30942   -0.92832 
#> sigma      0.44511     0.04041    0.37255    0.53179 
#> 

# Perform MCMC sampling using the named model (chains/iters reduced for speed)
mcmc_basic <- mod_basic$sample(sampling = 500, warmup = 500, chains = 2)
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
mcmc_basic$summary()
#> variable    mean    sd     map    q2.5   q97.5  ess_bulk  ess_tail  rhat 
#> lp        -49.80  1.41  -48.84  -53.55  -48.05       416       652  1.00 
#> mu[1]      -0.12  0.10   -0.11   -0.31    0.06      1106       715  1.00 
#> mu[2]       2.16  0.11    2.16    1.94    2.38       851       584  1.01 
#> mu[3]      -1.12  0.10   -1.12   -1.33   -0.91       857       700  1.00 
#> sigma       0.47  0.04    0.47    0.39    0.56       705       756  1.00 

# --- 2. Optional: Adding Custom Parameter Names and initial values ---
# You can optionally use 'par_names' to assign meaningful labels
# to vector or matrix elements for easier interpretation.
mod_named <- rtmb_model(
  data = data_list,
  code = model_code,
  init = list(mu = rep(0, 3), sigma = 1),
  par_names = list(mu = group_names)
)
#> Pre-checking model code...
#> Checking RTMB setup...

map_named <- mod_named$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 46.99
# The summary now displays: mu[Control], mu[Treatment_A], mu[Treatment_B]
map_named$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 46.99
#> Approx. Log Marginal Likelihood (Laplace): -52.63
#> 
#> Point Estimates and 95% Wald CI:
#>        variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> mu[Control]      -0.12319     0.09499   -0.30936    0.06298 
#> mu[Treatment_A]   2.15921     0.10806    1.94742    2.37099 
#> mu[Treatment_B]  -1.11887     0.09722   -1.30942   -0.92832 
#> sigma             0.44510     0.04041    0.37255    0.53178 
#> 

# Perform MCMC sampling using the named model (chains/iters reduced for speed)
mcmc_named <- mod_named$sample(sampling = 500, warmup = 500, chains = 2)
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
mcmc_named$summary()
#>        variable    mean    sd     map    q2.5   q97.5  ess_bulk  ess_tail  rhat 
#> lp               -49.84  1.43  -49.11  -53.41  -48.11       495       654  1.00 
#> mu[Control]       -0.12  0.10   -0.12   -0.32    0.08      1285       710  1.00 
#> mu[Treatment_A]    2.16  0.11    2.16    1.95    2.37      1014       733  1.00 
#> mu[Treatment_B]   -1.13  0.10   -1.13   -1.31   -0.92      1194       633  1.00 
#> sigma              0.47  0.04    0.47    0.39    0.56       964       576  1.00 
# }
```
