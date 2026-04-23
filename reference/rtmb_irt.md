# RTMB-based IRT (Item Response Theory) Wrapper

Performs Item Response Theory modeling (1PL, 2PL, 3PL, and Graded
Response Model) using RTMB. Missing values (NA) in the data are
automatically removed internally for efficient computation.

## Usage

``` r
rtmb_irt(
  data,
  model = c("2PL", "1PL", "3PL"),
  type = c("binary", "ordered"),
  prior = list(a_log_mean = 0, a_log_sd = 0.5, b_mean = 0, b_sd = 2.5, c_alpha = 1,
    c_beta = 4, theta_sd = 1),
  init = NULL
)
```

## Arguments

- data:

  A data frame or matrix of item responses (N persons x J items).

- model:

  Character string for the model type: "1PL", "2PL", or "3PL".

- type:

  Character string for the data type: "binary" or "ordered".

- prior:

  List of hyperparameters for prior distributions.

- init:

  List of initial values.

## Examples

``` r
# \donttest{
# --- 1. Binary Data (e.g., correct/incorrect answers) ---
# Simulate binary response data for 100 persons and 5 items
set.seed(123)
bin_data <- matrix(rbinom(500, size = 1, prob = 0.6), nrow = 100, ncol = 5)
colnames(bin_data) <- paste0("Item", 1:5)

# Introduce some missing values (NA) to demonstrate automatic handling
bin_data[sample(1:500, 10)] <- NA

# Fit a 2-Parameter Logistic (2PL) model
fit_2pl <- rtmb_irt(data = bin_data, model = "2PL", type = "binary")
#> Error: [Data error] Missing values (NA) are present in data 'Y'.
#>   * Solution: Please exclude NAs or impute missing values before passing to the model.

# Maximum A Posteriori (MAP) estimation
map_2pl <- fit_2pl$optimize()
#> Error: object 'fit_2pl' not found
map_2pl$summary()
#> Error: object 'map_2pl' not found

# --- 2. Ordered Data (e.g., Likert scale) ---
# Simulate ordered response data (categories 1 to 5)
ord_data <- matrix(sample(1:5, 500, replace = TRUE), nrow = 100, ncol = 5)
colnames(ord_data) <- paste0("Item", 1:5)

# Fit a Graded Response Model (2PL for ordered data)
fit_ord <- rtmb_irt(data = ord_data, model = "2PL", type = "ordered")
#> Pre-checking model code...
#> Checking RTMB setup...

map_ord <- fit_ord$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 844.69
map_ord$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 844.69
#> Approx. Log Marginal Likelihood (Laplace): -740.48
#> Note: Random effects are stored in $random_effects
#> 
#> Point Estimates and 95% Wald CI:
#>            variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> a[Item1]              0.56295     0.20237    0.27828    1.13882 
#> a[Item2]              0.45029     0.17162    0.21334    0.95043 
#> a[Item3]              0.49884     0.17987    0.24606    1.01129 
#> a[Item4]              0.53336     0.18806    0.26724    1.06449 
#> a[Item5]              0.61652     0.23564    0.29148    1.30402 
#> b[Item1,Threshold1]  -1.76683     0.28308   -2.32166   -1.21201 
#> b[Item2,Threshold1]  -1.40623     0.13643   -2.25247    0.66742 
#> b[Item3,Threshold1]  -1.40623     0.13643   -2.25247        Inf 
#> b[Item4,Threshold1]  -1.40623     0.13643   -2.25247        Inf 
#> b[Item5,Threshold1]  -1.40620     0.13656   -2.25247        Inf 
#> 

# MCMC sampling for the ordered model (chains and iterations reduced)
mcmc_ord <- fit_ord$sample(sampling = 500, warmup = 500, chains = 2)
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
mcmc_ord$summary()
#>            variable      mean     sd       map      q2.5    q97.5  ess_bulk  ess_tail  rhat 
#> lp                   -1014.56  10.09  -1012.97  -1034.60  -994.86       262       545  1.00 
#> a[Item1]                 0.62   0.22      0.52      0.29     1.11       904       893  1.00 
#> a[Item2]                 0.53   0.21      0.40      0.24     1.03       782       796  1.00 
#> a[Item3]                 0.58   0.21      0.52      0.26     1.06       778       475  1.00 
#> a[Item4]                 0.62   0.21      0.53      0.30     1.10       905       751  1.00 
#> a[Item5]                 0.78   0.32      0.71      0.31     1.56       489       739  1.00 
#> b[Item1,Threshold1]     -1.98   0.24     -1.94     -2.52    -1.56       702       471  1.00 
#> b[Item2,Threshold1]     -1.60   0.16     -1.57     -1.93    -1.29      1178       751  1.00 
#> b[Item3,Threshold1]     -1.48   0.15     -1.43     -1.78    -1.19      1113       804  1.00 
#> b[Item4,Threshold1]     -1.36   0.15     -1.33     -1.66    -1.08      1034       864  1.00 
# }
```
