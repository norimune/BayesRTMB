
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

log_prob <- function(data, param) {
  lp <- 0
  lp <- lp + lpdf$normal(param$mu, 0, 1)
  lp <- lp + lpdf$normal(data$y, param$mu, 1)
  lp
}

data <- list(
  y = c(0.3, 0.7, -0.1, 0.5)
)

par_list <- list(
  mu = Dim()
)

model <- RTMB_Model$new(
  data = data,
  par_list = par_list,
  log_prob = log_prob
)

fit_map <- model$optimize()
fit_mcmc <- model$sample(chains = 2, warmup = 500)
```

## Posterior plots

``` r
draws <- fit_mcmc$draws()
plot_trace(draws)
plot_dens(draws)
plot_ac(draws, 1)
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
