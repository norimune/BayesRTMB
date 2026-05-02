# Specify a uniform or manual prior

Specify a uniform or manual prior

## Usage

``` r
prior_uniform(
  Int_sd = NULL,
  beta_sd = NULL,
  sigma_rate = NULL,
  tau_rate = NULL,
  ...
)
```

## Arguments

- Int_sd:

  Standard deviation for the intercept prior (Normal). Default is NULL
  (flat).

- beta_sd:

  Standard deviation for the coefficients prior (Normal). Default is
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
