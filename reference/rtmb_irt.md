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
  prior = prior_flat(),
  init = NULL,
  fixed = NULL,
  view = NULL,
  missing = c("listwise", "fiml")
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

  Prior configuration: \`prior_flat()\`, \`prior_normal()\`, or
  \`prior_weak()\`. Hyperparameters can be specified within these
  functions (e.g., \`prior_normal(b_sd = 5)\`). Available parameters for
  IRT: \`a_rate\` (discrimination), \`b_sd\` (difficulty),
  \`c_alpha\`/\`c_beta\` (guessing).

- init:

  List of initial values.

- fixed:

  A named list of parameter values to fix (optional).

- view:

  Character vector of parameter names to prioritize in summary.

- missing:

  Missing value handling strategy: "listwise" (default) or "fiml" (Full
  Information Maximum Likelihood).

## Examples

``` r

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
#> Starting RTMB optimization...
#> 
#> Warning: Optimization did not converge ( convergence code = 1; message = function evaluation limit reached without convergence (9)). Estimates may be unreliable; consider increasing num_estimate, changing initial values, or adjusting optimizer control settings.
  map_2pl$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 334.06
#> Approx. Log Marginal Likelihood (Laplace): -326.51
#> Note: Random effects are stored in $random_effects (use ranef = TRUE to show them)
#> 
#> Point Estimates and 95% Wald CI:
#> variable     Estimate   Std. Error     Lower 95%       Upper 95% 
#> a[Item1]      0.00017      0.00222       0.00000  14803403.02410 
#> a[Item2]      0.00024      0.00268       0.00000    524214.81653 
#> a[Item3]      0.33853      0.63664       0.00849        13.50088 
#> a[Item4]      0.52727      0.55116       0.06796         4.09068 
#> a[Item5]      0.33846      0.63717       0.00845        13.55001 
#> b[Item1]  -2340.60797  30040.17736  -61218.27250     56537.05656 
#> b[Item2]  -2001.47440  21926.82445  -44977.26074     40974.31193 
#> b[Item3]     -1.35900      2.49182      -6.24287         3.52488 
#> b[Item4]     -0.90809      0.92095      -2.71311         0.89693 
#> b[Item5]     -1.35941      2.49533      -6.25018         3.53135 
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
#>     theta ~ normal(0, 1)
#>   }
#> )

  map_ord <- fit_ord$optimize()
#> Starting RTMB optimization...
#> 
  map_ord$summary()
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 652.00
#> Approx. Log Marginal Likelihood (Laplace): -661.02
#> Note: Random effects are stored in $random_effects (use ranef = TRUE to show them)
#> 
#> Point Estimates and 95% Wald CI:
#>            variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> a[Item1]              1.08980     0.33157    0.60030    1.97843 
#> a[Item2]              1.35449     0.39466    0.76518    2.39767 
#> a[Item3]              0.84779     0.28248    0.44124    1.62893 
#> a[Item4]              1.29604     0.37514    0.73491    2.28560 
#> a[Item5]              0.98734     0.31499    0.52834    1.84510 
#> b[Item1,Threshold1]  -1.80041     0.33328   -2.45362   -1.14720 
#> b[Item2,Threshold1]  -1.74260     0.36428   -2.45658   -1.02862 
#> b[Item3,Threshold1]  -1.31729     0.27357   -1.85348   -0.78110 
#> b[Item4,Threshold1]  -1.20024     0.30906   -1.80599   -0.59450 
#> b[Item5,Threshold1]  -1.53616     0.30305   -2.13012   -0.94220 
#> 

  # Note: For complex models like the Graded Response Model, the Wald CI from optimize()
  # may become extremely wide due to parameter transformations.
  # MCMC sampling is recommended for reliable interval estimation.

  # MCMC sampling for the ordered model (chains and iterations reduced)
  if (FALSE) { # \dontrun{
  mcmc_ord <- fit_ord$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_ord$summary()
  mcmc_ord$summary()
  } # }
```
