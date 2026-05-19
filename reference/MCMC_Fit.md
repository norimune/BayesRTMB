# MCMC fit object

An R6 class storing posterior samples and related information from MCMC
estimation.

Bayes factors are computed as ratios of marginal likelihoods estimated
by bridge sampling. The comparison model can be constructed by fixing
parameters with \`fixed = list(...)\`, or supplied as an already fitted
\`MCMC_Fit\` object via \`comparison_fit\`. Bayes factors are computed
only by marginal-likelihood model comparison.

## Super class

[`RTMB_Fit_Base`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.md)
-\> `mcmc_fit`

## Public fields

- `model`:

  An \`RTMB_Model\` object used for estimation.

- `fit`:

  Posterior draws for model parameters.

- `random_fit`:

  Posterior draws for random effects.

- `transform_fit`:

  Posterior draws for transformed parameters.

- `transform_dims`:

  Dimension information for transformed parameters.

- `generate_fit`:

  Posterior draws for generated quantities.

- `generate_dims`:

  Dimension information for generated quantities.

- `eps`:

  Step size used by the sampler.

- `accept`:

  Acceptance statistics from sampling.

- `treedepth`:

  Tree depth used in HMC/NUTS sampling.

- `laplace`:

  Logical; whether Laplace approximation was used.

- `posterior_mean`:

  Posterior mean estimates.

- `log_ml`:

  Numeric value storing the calculated log marginal likelihood from
  bridge sampling.

- `comparison_fit`:

  An `MCMC_Fit` object containing the fitted comparison model.

## Methods

### Public methods

- [`mcmc_fit$get_point_estimate()`](#method-mcmc_fit-get_point_estimate)

- [`mcmc_fit$new()`](#method-mcmc_fit-initialize)

- [`mcmc_fit$print()`](#method-mcmc_fit-print)

- [`mcmc_fit$draws()`](#method-mcmc_fit-draws)

- [`mcmc_fit$summary()`](#method-mcmc_fit-summary)

- [`mcmc_fit$unconstrain_draws()`](#method-mcmc_fit-unconstrain_draws)

- [`mcmc_fit$log_prob()`](#method-mcmc_fit-log_prob)

- [`mcmc_fit$bridgesampling()`](#method-mcmc_fit-bridgesampling)

- [`mcmc_fit$bayes_factor()`](#method-mcmc_fit-bayes_factor)

- [`mcmc_fit$transformed_draws()`](#method-mcmc_fit-transformed_draws)

- [`mcmc_fit$generated_quantities()`](#method-mcmc_fit-generated_quantities)

- [`mcmc_fit$resolve_switching()`](#method-mcmc_fit-resolve_switching)

- [`mcmc_fit$clone()`](#method-mcmc_fit-clone)

Inherited methods

- [`RTMB_Fit_Base$EAP()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-EAP)
- [`RTMB_Fit_Base$MAP()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-MAP)
- [`RTMB_Fit_Base$estimate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-estimate)
- [`RTMB_Fit_Base$fa_rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-fa_rotate)
- [`RTMB_Fit_Base$rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-rotate)

------------------------------------------------------------------------

### `mcmc_fit$get_point_estimate()`

Get point estimate for a target parameter.

#### Usage

    mcmc_fit$get_point_estimate(target, chains = NULL, best_chains = NULL)

#### Arguments

- `target`:

  Target parameter name.

- `chains`:

  Numeric vector of chains to include. If NULL, all chains are used.

- `best_chains`:

  Integer; number of best chains to retain based on mean log-posterior.

#### Returns

Matrix, array, vector, or scalar point estimate.

------------------------------------------------------------------------

### `mcmc_fit$new()`

Create a new \`MCMC_Fit\` object.

#### Usage

    mcmc_fit$new(
      model,
      fit,
      random_fit,
      eps,
      accept,
      treedepth,
      laplace,
      posterior_mean
    )

#### Arguments

- `model`:

  An \`RTMB_Model\` object used for estimation.

- `fit`:

  Posterior draws for model parameters.

- `random_fit`:

  Posterior draws for random effects, if available.

- `eps`:

  Step size used by the sampler.

- `accept`:

  Acceptance statistics from sampling.

- `treedepth`:

  Tree depth used in HMC/NUTS sampling.

- `laplace`:

  Logical; whether Laplace approximation was used.

- `posterior_mean`:

  Posterior mean estimates.

------------------------------------------------------------------------

### `mcmc_fit$print()`

Print a brief summary of the fitted object.

#### Usage

    mcmc_fit$print(...)

#### Arguments

- `...`:

  Additional arguments.

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### `mcmc_fit$draws()`

Extract posterior draws for selected parameters.

#### Usage

    mcmc_fit$draws(
      pars = NULL,
      chains = NULL,
      best_chains = NULL,
      inc_random = FALSE,
      inc_transform = TRUE,
      inc_generate = TRUE
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

  Integer; number of best chains to retain based on mean log-posterior
  (lp).

- `inc_random`:

  Logical; whether to include random effects in the output. Default is
  FALSE.

- `inc_transform`:

  Logical; whether to include transformed parameters in the output.
  Default is TRUE.

- `inc_generate`:

  Logical; whether to include generated quantities in the output.
  Default is TRUE.

#### Returns

Posterior draws.

------------------------------------------------------------------------

### `mcmc_fit$summary()`

Summarize posterior draws.

#### Usage

    mcmc_fit$summary(
      pars = NULL,
      chains = NULL,
      best_chains = NULL,
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

- `chains`:

  Numeric vector specifying the chains to extract. If NULL, draws from
  all chains are used.

- `best_chains`:

  Integer; number of best chains to retain based on mean log-posterior
  (lp).

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

A summary object.

------------------------------------------------------------------------

### `mcmc_fit$unconstrain_draws()`

Transform posterior draws to the unconstrained scale.

#### Usage

    mcmc_fit$unconstrain_draws()

#### Returns

Posterior draws on the unconstrained scale.

------------------------------------------------------------------------

### `mcmc_fit$log_prob()`

Evaluate log-probability values.

#### Usage

    mcmc_fit$log_prob(safe = FALSE)

#### Arguments

- `safe`:

  Logical; whether to wrap the evaluation in a tryCatch block. Default
  is FALSE for speed.

#### Returns

Numeric vector of log-probability values.

------------------------------------------------------------------------

### `mcmc_fit$bridgesampling()`

Estimate the marginal likelihood by bridge sampling.

#### Usage

    mcmc_fit$bridgesampling(
      method = "normal",
      use_neff = TRUE,
      seed = NULL,
      max_iter = 100
    )

#### Arguments

- `method`:

  Character; the method to use for bridge sampling (e.g., "warp3",
  "normal"). Default is "warp3".

- `use_neff`:

  Logical; whether to use the effective sample size (ESS) to adjust for
  autocorrelation. Default is TRUE.

- `seed`:

  Integer; random seed for reproducibility. Default is NULL.

- `max_iter`:

  Integer; maximum number of iterations for the estimation algorithm.
  Default is 100.

#### Returns

Bridge sampling result.

------------------------------------------------------------------------

### `mcmc_fit$bayes_factor()`

Calculate the Bayes factor by marginal-likelihood model comparison.

#### Usage

    mcmc_fit$bayes_factor(
      fixed = NULL,
      comparison_fit = NULL,
      bs_method = "normal",
      error_threshold = 0.2,
      ...
    )

#### Arguments

- `fixed`:

  Named list of parameter values used to construct the comparison model.
  For example, `fixed = list(delta = 0)` or `fixed = list("b[x]" = 0)`.

- `comparison_fit`:

  Optional \`MCMC_Fit\` object for an already fitted comparison model.

- `bs_method`:

  Character; the method to use for bridge sampling ("normal" or
  "warp3").

- `error_threshold`:

  Numeric; threshold for the approximate error warning.

- `...`:

  Additional arguments passed to \`sample()\` when fitting the
  comparison model (e.g., `chains = 4`, `sampling = 4000`).

#### Returns

A list of class `bayes_factor_rtmb` containing Bayes factor results.

------------------------------------------------------------------------

### `mcmc_fit$transformed_draws()`

Compute transformed parameters from posterior draws.

#### Usage

    mcmc_fit$transformed_draws(tran_fn = NULL)

#### Arguments

- `tran_fn`:

  A function for transformed parameters.

#### Returns

Transformed parameter draws.

------------------------------------------------------------------------

### `mcmc_fit$generated_quantities()`

Compute generated quantities from posterior draws.

#### Usage

    mcmc_fit$generated_quantities(code)

#### Arguments

- `code`:

  An \`rtmb_code({ ... })\` or \`{ ... }\` block containing the logic to
  be calculated using posterior samples.

#### Returns

The \`MCMC_Fit\` object itself (invisibly). Results are appended to the
\`generate_fit\` field.

------------------------------------------------------------------------

### `mcmc_fit$resolve_switching()`

Resolve label switching in posterior draws.

#### Usage

    mcmc_fit$resolve_switching(
      target,
      linked = NULL,
      overwrite = TRUE,
      scalar_fns = list()
    )

#### Arguments

- `target`:

  Character string specifying the target variable to base the relabeling
  on.

- `linked`:

  Character vector of variable names to be relabeled in the same order
  as the target. Default is NULL.

- `overwrite`:

  Logical; whether to overwrite the stored draws in the current object.
  Default is TRUE.

- `scalar_fns`:

  A named list of functions to apply to scalar variables for relabeling.
  Default is an empty list.

#### Returns

Relabeled draws or updated object.

------------------------------------------------------------------------

### `mcmc_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    mcmc_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
