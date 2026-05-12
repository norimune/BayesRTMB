# Specify a uniform or manual prior

Specify a uniform or manual prior

## Usage

``` r
prior_uniform(
  Intercept_sd = NULL,
  b_sd = NULL,
  mu_sd = NULL,
  sigma_rate = NULL,
  tau_rate = NULL,
  ...
)
```

## Arguments

- Intercept_sd:

  Standard deviation for the intercept prior (Normal). Default is NULL
  (flat).

- b_sd:

  Standard deviation for the coefficients prior (Normal). Default is
  NULL (flat).

- mu_sd:

  Standard deviation for the mean/intercept prior (Normal). Default is
  NULL (flat).

- sigma_rate:

  Rate for the residual standard deviation prior (Exponential). Default
  is NULL (flat).

- tau_rate:

  Rate for the random effects standard deviation prior (Exponential).
  Default is NULL (flat).

- ...:

  Optional hyperparameters

## Value

A list with class "rtmb_prior"
