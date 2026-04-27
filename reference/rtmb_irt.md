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
  #bin_data[sample(1:500, 10)] <- NA

  # Fit a 2-Parameter Logistic (2PL) model
  fit_2pl <- rtmb_irt(data = bin_data, model = "2PL", type = "binary")
#> Pre-checking model code...
#> Checking RTMB setup...

  # Maximum A Posteriori (MAP) estimation
  map_2pl <- fit_2pl$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 348.19
  map_2pl$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 348.19
#> Approx. Log Marginal Likelihood (Laplace): -351.42
#> Note: Random effects are stored in $random_effects
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> a[Item1]   0.50626     0.19539    0.23760    1.07868 
#> a[Item2]   0.53070     0.20473    0.24916    1.13036 
#> a[Item3]   0.53650     0.20799    0.25095    1.14698 
#> a[Item4]   0.68479     0.26912    0.31698    1.47936 
#> a[Item5]   0.53477     0.20714    0.25030    1.14253 
#> b[Item1]  -0.81807     0.49804   -1.79421    0.15808 
#> b[Item2]  -0.94996     0.50981   -1.94917    0.04926 
#> b[Item3]  -0.86002     0.48909   -1.81862    0.09857 
#> b[Item4]  -0.70531     0.39355   -1.47665    0.06604 
#> b[Item5]  -0.86159     0.48990   -1.82177    0.09859 
#> 

  # --- 2. Ordered Data (e.g., Likert scale) ---
  # Simulate ordered response data (categories 1 to 5) with a true latent structure
  set.seed(123)
  N <- 100
  J <- 5
  K <- 4

  # True parameters
  theta <- rnorm(N, 0, 1)          # Person abilities (latent traits)
  a <- runif(J, 0.8, 2.0)          # Item discriminations

  # Item-specific thresholds (must be strictly increasing)
  b <- matrix(NA, nrow = J, ncol = K - 1)
  for (j in 1:J) {
    b[j, ] <- sort(c(-1.5, 0, 1.5) + rnorm(3, 0, 0.2))
  }

  ord_data <- matrix(NA, nrow = N, ncol = J)
  colnames(ord_data) <- paste0("Item", 1:J)

  # Generate responses based on the Graded Response Model
  for (i in 1:N) {
    for (j in 1:J) {
      # Latent continuous response (eta + standard logistic noise)
      eta <- a[j] * theta[i]
      y_star <- eta + rlogis(1)

      # Apply thresholds to determine the observed category (1 to 4)
      category <- 1
      for (k in 1:(K - 1)) {
        if (y_star > b[j, k]) category <- k + 1
      }
      ord_data[i, j] <- category
    }
  }

  # Fit a Graded Response Model (2PL for ordered data)
  fit_ord <- rtmb_irt(data = ord_data, model = "2PL", type = "ordered")
#> Pre-checking model code...
#> Checking RTMB setup...
  fit_ord$print_code()
#> === RTMB Model Code ===
#> 
#> rtmb_code(
#>   setup = {
#>     obs_data <- which(!is.na(Y), arr.ind = TRUE)
#>     person_idx <- as.integer(obs_data[, "row"])
#>     item_idx <- as.integer(obs_data[, "col"])
#>     Y_obs <- Y[obs_data]
#>     N_persons <- nrow(Y)
#>     N_items <- ncol(Y)
#>     N_obs <- length(Y_obs)
#>     if (min(Y_obs) == 0) 
#>       Y_obs <- Y_obs + 1
#>     K_cat <- max(Y_obs)
#>   }, 
#>   parameters = {
#>     theta <- Dim(N_persons, random = TRUE)
#>     b <- Dim(c(N_items, K_cat - 1), type = "ordered")
#>     a <- Dim(N_items, lower = 0)
#>   }, 
#>   model = {
#>     # Likelihood
#>     for (i in 1:N_obs) {
#>       p <- person_idx[i]
#>       j <- item_idx[i]
#>       y <- Y_obs[i]
#>       eta <- a[j] * theta[p]
#>       y ~ ordered_logistic(eta, b[j, ])
#>     }
#>     # Priors
#>     a ~ lognormal(prior_a_log_mean, prior_a_log_sd)
#>     for (j in 1:N_items) b[j, ] ~ normal(prior_b_mean, prior_b_sd)
#>     theta ~ normal(0, prior_theta_sd)
#>   }
#> )

  map_ord <- fit_ord$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 690.99
  map_ord$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 690.99
#> Approx. Log Marginal Likelihood (Laplace): -731.73
#> Note: Random effects are stored in $random_effects
#> 
#> Point Estimates and 95% Wald CI:
#>            variable  Estimate  Std. Error  Lower 95%   Upper 95% 
#> a[Item1]              0.96615     0.26759    0.56143     1.66264 
#> a[Item2]              1.18876     0.31684    0.70506     2.00432 
#> a[Item3]              0.75313     0.20830    0.43797     1.29505 
#> a[Item4]              1.00441     0.24271    0.62550     1.61287 
#> a[Item5]              1.02515     0.26812    0.61399     1.71162 
#> b[Item1,Threshold1]  -1.69598     0.30728   -2.29823    -1.09373 
#> b[Item2,Threshold1]  -1.62229     0.32430   -2.29823  5924.72682 
#> b[Item3,Threshold1]  -1.51011     0.24906   -2.29809  6011.05137 
#> b[Item4,Threshold1]  -1.11941     0.20116   -2.20051  6012.61580 
#> b[Item5,Threshold1]  -1.11941     0.20116   -2.20051         Inf 
#> 

  # Note: For complex models like the Graded Response Model, the Wald CI from optimize()
  # may become extremely wide due to parameter transformations.
  # MCMC sampling is recommended for reliable interval estimation.

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
#>            variable     mean     sd      map     q2.5    q97.5  ess_bulk  ess_tail  rhat 
#> lp                   -815.26  11.12  -816.85  -839.59  -794.88       200       451  1.00 
#> a[Item1]                1.12   0.30     1.03     0.60     1.74       410       515  1.00 
#> a[Item2]                1.40   0.35     1.31     0.77     2.15       285       484  1.00 
#> a[Item3]                0.77   0.21     0.72     0.40     1.22       509       692  1.00 
#> a[Item4]                1.17   0.28     1.23     0.67     1.70       359       726  1.00 
#> a[Item5]                1.15   0.32     1.02     0.61     1.81       237       757  1.00 
#> b[Item1,Threshold1]    -1.97   0.25    -1.93    -2.45    -1.49       323       544  1.00 
#> b[Item2,Threshold1]    -1.72   0.22    -1.80    -2.17    -1.29       406       690  1.00 
#> b[Item3,Threshold1]    -1.47   0.20    -1.47    -1.86    -1.07       380       690  1.01 
#> b[Item4,Threshold1]    -1.17   0.19    -1.21    -1.55    -0.80       534       840  1.00 
 # }
```
