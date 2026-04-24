# VB fit object

VB fit object

VB fit object

## Details

An R6 class storing posterior samples and related information from
Automatic Differentiation Variational Inference (ADVI).

## Super class

[`BayesRTMB::RTMB_Fit_Base`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.md)
-\> `advi_fit`

## Public fields

- `model`:

  An \`RTMB_Model\` object used for estimation.

- `fit`:

  A 3D array of posterior draws for fixed model parameters.

- `random_fit`:

  A 3D array of posterior draws for random effects, if available.

- `transform_fit`:

  A 3D array of posterior draws for transformed parameters, if
  available.

- `generate_fit`:

  A 3D array of posterior draws for generated quantities, if available.

- `transform_dims`:

  A list storing dimension information for transformed parameters.

- `generate_dims`:

  A list storing dimension information for generated quantities.

- `elbo_history`:

  A list of numeric vectors storing the Evidence Lower Bound (ELBO)
  history during optimization for each chain.

- `laplace`:

  Logical; whether Laplace approximation was used to marginalize random
  effects.

- `posterior_mean`:

  A named numeric vector of posterior mean estimates.

- `ELBO`:

  A numeric vector of final ELBO values for each chain.

- `rel_obj_vals`:

  A numeric vector of final relative objective tolerance values for each
  chain.

- `best_chain`:

  Integer; the index of the chain with the maximum ELBO.

- `mu_history`:

  Matrix of the parameter trajectory from the final window.

## Methods

### Public methods

- [`VB_Fit$get_point_estimate()`](#method-advi_fit-get_point_estimate)

- [`VB_Fit$new()`](#method-advi_fit-new)

- [`VB_Fit$print()`](#method-advi_fit-print)

- [`VB_Fit$draws()`](#method-advi_fit-draws)

- [`VB_Fit$summary()`](#method-advi_fit-summary)

- [`VB_Fit$plot_elbo()`](#method-advi_fit-plot_elbo)

- [`VB_Fit$plot_trajectory()`](#method-advi_fit-plot_trajectory)

- [`VB_Fit$transformed_draws()`](#method-advi_fit-transformed_draws)

- [`VB_Fit$generated_quantities()`](#method-advi_fit-generated_quantities)

- [`VB_Fit$clone()`](#method-advi_fit-clone)

Inherited methods

- [`BayesRTMB::RTMB_Fit_Base$EAP()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-EAP)
- [`BayesRTMB::RTMB_Fit_Base$MAP()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-MAP)
- [`BayesRTMB::RTMB_Fit_Base$fa_rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-fa_rotate)
- [`BayesRTMB::RTMB_Fit_Base$rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-rotate)

------------------------------------------------------------------------

### Method `get_point_estimate()`

Get point estimate for a target parameter (internal use).

#### Usage

    VB_Fit$get_point_estimate(target)

#### Arguments

- `target`:

  Target parameter name.

#### Returns

Matrix or array of point estimate.

------------------------------------------------------------------------

### Method `new()`

Create a new \`VB_Fit\` object.

#### Usage

    VB_Fit$new(
      model,
      fit,
      random_fit,
      elbo_history,
      laplace,
      posterior_mean,
      ELBO,
      rel_obj_vals,
      best_chain,
      mu_history
    )

#### Arguments

- `model`:

  An \`RTMB_Model\` object.

- `fit`:

  A 3D array of parameter draws.

- `random_fit`:

  A 3D array of random effect draws.

- `elbo_history`:

  A list of numeric vectors of ELBO values for each chain.

- `laplace`:

  Logical; indicates if Laplace approximation was used.

- `posterior_mean`:

  A named numeric vector of posterior means.

- `ELBO`:

  A numeric vector of final ELBO values for each chain.

- `rel_obj_vals`:

  A numeric vector of final relative objective tolerance values for each
  chain.

- `best_chain`:

  Integer; the index of the chain with the maximum ELBO.

- `mu_history`:

  Matrix of the parameter trajectory from the final window.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a brief summary of the fitted object.

#### Usage

    VB_Fit$print(...)

#### Arguments

- `...`:

  Additional arguments passed to the \`summary\` method.

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### Method `draws()`

Extract posterior draws for selected parameters.

#### Usage

    VB_Fit$draws(
      pars = NULL,
      chains = NULL,
      best_chains = NULL,
      inc_random = FALSE,
      inc_transform = TRUE,
      inc_generate = TRUE,
      best_only = FALSE
    )

#### Arguments

- `pars`:

  Character or numeric vector specifying the names or indices of
  parameters to extract. If NULL, all available parameters are
  extracted.

- `chains`:

  Numeric vector specifying the chains to extract. If NULL, draws from
  all chains are returned.

- `best_chains`:

  Integer; number of best chains to retain based on ELBO. If supplied,
  chains with the highest ELBO are retained.

- `inc_random`:

  Logical; whether to include random effects in the output. Default is
  FALSE.

- `inc_transform`:

  Logical; whether to include transformed parameters in the output.
  Default is TRUE.

- `inc_generate`:

  Logical; whether to include generated quantities in the output.
  Default is TRUE.

- `best_only`:

  Logical; whether to extract only from the chain with the maximum ELBO.
  Default is FALSE unless explicitly requested.

#### Returns

A 3D array of posterior draws \`\[iterations, chains, parameters\]\`.

------------------------------------------------------------------------

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Summarize posterior draws.

#### Usage

    VB_Fit$summary(
      pars = NULL,
      max_rows = 10,
      digits = 2,
      inc_random = FALSE,
      inc_transform = TRUE,
      inc_generate = TRUE
    )

#### Arguments

- `pars`:

  Character or numeric vector specifying the names or indices of
  parameters to summarize. If NULL, all available parameters are
  summarized.

- `max_rows`:

  Integer; maximum number of rows to print in the summary table. Default
  is 10.

- `digits`:

  Integer; number of decimal places to print. Default is 2.

- `inc_random`:

  Logical; whether to include random effects in the summary. Default is
  FALSE.

- `inc_transform`:

  Logical; whether to include transformed parameters in the summary.
  Default is TRUE.

- `inc_generate`:

  Logical; whether to include generated quantities in the summary.
  Default is TRUE.

#### Returns

A data frame containing the summarized posterior statistics.

------------------------------------------------------------------------

### Method `plot_elbo()`

Plot the ELBO history to diagnose convergence.

#### Usage

    VB_Fit$plot_elbo(tail_n = 1000, ests = NULL, type = "l", ...)

#### Arguments

- `tail_n`:

  Integer; the number of recent iterations to plot. If NULL, plots the
  entire history. Default is 2000.

- `ests`:

  Character string \`"best"\`, numeric vector of estimate indices (e.g.,
  \`c(1, 3)\`), or \`NULL\` to plot all. Default is \`NULL\`.

- `type`:

  Character string; the type of plot. Default is "l" (lines).

- `...`:

  Additional arguments passed to the \`plot\` function.

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### Method `plot_trajectory()`

Plot the parameter trajectory from the final optimization window.

#### Usage

    VB_Fit$plot_trajectory(pars = NULL, type = "l", ...)

#### Arguments

- `pars`:

  Character vector specifying the names of parameters to plot. If NULL,
  plots all available parameters.

- `type`:

  Character string; the type of plot. Default is "l" (lines).

- `...`:

  Additional arguments passed to the \`matplot\` function.

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### Method `transformed_draws()`

Compute transformed parameters from posterior draws.

#### Usage

    VB_Fit$transformed_draws(tran_fn = NULL)

#### Arguments

- `tran_fn`:

  An optional user-supplied function that takes data and parameter lists
  to return transformed quantities.

#### Returns

The \`VB_Fit\` object itself, invisibly.

------------------------------------------------------------------------

### Method `generated_quantities()`

Compute generated quantities from posterior draws.

#### Usage

    VB_Fit$generated_quantities(code)

#### Arguments

- `code`:

  An \`rtmb_code( ... )\` or \` ... \` block containing the logic to be
  calculated using posterior samples.

#### Returns

The \`VB_Fit\` object itself (invisibly). Results are appended to the
\`generate_fit\` field.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    VB_Fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
