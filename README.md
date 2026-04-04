
# BayesRTMB

BayesRTMB is an R package for Bayesian modeling with **RTMB**.

The package provides:

- model objects built with `RTMB_Model`
- MAP estimation with `optimize()`
- MCMC sampling with `sample()`
- helper functions such as `lpdf` and `math`
- plotting functions for posterior samples

## Installation

You can install the development version from GitHub with:

``` r
remotes::install_github("YOUR_GITHUB_NAME/BayesRTMB")
```

## Basic idea

A model is defined by:

- a data list
- a parameter specification list
- a user-defined `log_prob()` function

Then you can optimize the posterior or draw samples.

## Example

``` r
library(BayesRTMB)

Trial <- 10
Y <- 6

dat <- list(Trial = Trial, Y = Y)

par <- list(
  theta = Dim(lower=0,upper=1)
)


log_prob <- function(dat,par){
  getAll(dat, par)
  lp <- binomial_lpmf(Y,Trial,par$theta)
  lp <- lp + beta_lpdf(par$theta, 1, 1)
  return(lp)
}

mdl <- rtmb_model(dat, par, log_prob)

fit.mcmc <- mdl$sample()
fit.mcmc

fit.mcmc$bridgesampling()

fit.map <- mdl$optimize()
fit.map

fit.mcmc$draws("theta") |> plot_dens()
fit.mcmc$draws("theta") |> plot_trace()
fit.mcmc$draws("theta") |> plot_acf()
```

## Posterior plots

``` r
fit.mcmc$draws("theta") |> plot_dens()
fit.mcmc$draws("theta") |> plot_trace()
fit.mcmc$draws("theta") |> plot_acf()
```

## Main functions

- `RTMB_Model`
- `MAP_Fit`
- `MCMC_Fit`
- `lpdf`
- `math`
- `plot_trace()`
- `plot_dens()`
- `plot_ac()`

## Development status

BayesRTMB is under active development. The interface and function names
may change in future versions.
