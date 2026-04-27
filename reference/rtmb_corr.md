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
if (FALSE) { # \dontrun{
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
} # }
```
