# Specify normal/exponential priors for MAP and Bayesian inference

\`prior_normal()\` specifies normal priors for location parameters and
exponential priors for scale parameters. It is intended for MAP and
Bayesian inference, not for \`classic()\`.

## Usage

``` r
prior_normal(
  Intercept_sd = 10,
  mu_sd = 10,
  b_sd = 10,
  sigma_rate = 5,
  tau_rate = 5,
  ...
)
```

## Arguments

- Intercept_sd:

  Standard deviation for the intercept prior. If \`NULL\`, no intercept
  prior is added.

- mu_sd:

  Standard deviation for mean/intercept priors. If \`NULL\`, no mean
  prior is added.

- b_sd:

  Standard deviation for coefficient priors. If \`NULL\`, no coefficient
  prior is added.

- sigma_rate:

  Rate for residual standard deviation priors. If \`NULL\`, no sigma
  prior is added.

- tau_rate:

  Rate for random-effect standard deviation priors. If \`NULL\`, no tau
  prior is added.

- ...:

  Optional wrapper-specific hyperparameters.

## Value

A list with class \`"rtmb_prior"\` and \`type = "normal"\`.
