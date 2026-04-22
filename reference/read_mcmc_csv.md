# Restore MCMC Fit from CSV

Restore MCMC Fit from CSV

## Usage

``` r
read_mcmc_csv(model, name, dir = "BayesRTMB_mcmc", chains = 4, laplace = FALSE)
```

## Arguments

- model:

  An RTMB_Model object.

- name:

  Base name of the saved CSVs.

- dir:

  Directory where CSVs are saved. Default is "BayesRTMB_mcmc".

- chains:

  Number of chains. Default is 4.

- laplace:

  Logical; whether Laplace approximation was used. Default is FALSE.

## Value

An MCMC_Fit object.
