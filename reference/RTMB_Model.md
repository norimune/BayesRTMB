# Wrapper function to create an RTMB_Model instance

Combines user-defined data and model code (likelihood and priors) to
create an \`RTMB_Model\` (R6 class) instance for Bayesian inference
(MCMC, VB, MAP estimation).

## Usage

``` r
rtmb_model(
  data,
  code,
  par_names = list(),
  init = NULL,
  view = NULL,
  null_target = NULL
)
```

## Arguments

- data:

  A named list containing observation data and constants (e.g., sample
  size) used in the model.

- code:

  A model definition block described by \`rtmb_code(...)\` (including
  data, parameters, model, transform, generate).

- par_names:

  A list of specific variable names corresponding to the dimensions of
  each parameter (optional).

- init:

  A list or numeric vector of initial values for parameters (optional).
  If not specified, initialized randomly.

- view:

  Character vector of parameter names to be displayed preferentially at
  the top when outputting results like \`summary()\` (optional).

- null_target:

  Character string. To simultaneously create a null model, specify the
  target parameter to fix and the prior distribution to disable as a
  formula string (e.g., `"delta ~ cauchy(0, r)"`). Default is `NULL`.

## Value

An \`RTMB_Model\` class instance with compiled and pre-tested model.
