# MCMC fit object

MCMC fit object

MCMC fit object

## Details

An R6 class storing posterior samples and related information from MCMC
estimation.

## Super class

[`BayesRTMB::RTMB_Fit_Base`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.md)
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

- `null_fit`:

  An `MCMC_Fit` object containing the fitted null model. This is
  automatically cached when calculating a Bayes factor using a target
  string.

## Methods

### Public methods

- [`MCMC_Fit$get_point_estimate()`](#method-mcmc_fit-get_point_estimate)

- [`MCMC_Fit$new()`](#method-mcmc_fit-new)

- [`MCMC_Fit$print()`](#method-mcmc_fit-print)

- [`MCMC_Fit$draws()`](#method-mcmc_fit-draws)

- [`MCMC_Fit$summary()`](#method-mcmc_fit-summary)

- [`MCMC_Fit$unconstrain_draws()`](#method-mcmc_fit-unconstrain_draws)

- [`MCMC_Fit$log_prob()`](#method-mcmc_fit-log_prob)

- [`MCMC_Fit$bridgesampling()`](#method-mcmc_fit-bridgesampling)

- [`MCMC_Fit$bayes_factor()`](#method-mcmc_fit-bayes_factor)

- [`MCMC_Fit$transformed_draws()`](#method-mcmc_fit-transformed_draws)

- [`MCMC_Fit$generated_quantities()`](#method-mcmc_fit-generated_quantities)

- [`MCMC_Fit$resolve_switching()`](#method-mcmc_fit-resolve_switching)

- [`MCMC_Fit$clone()`](#method-mcmc_fit-clone)

Inherited methods

- [`BayesRTMB::RTMB_Fit_Base$EAP()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-EAP)
- [`BayesRTMB::RTMB_Fit_Base$MAP()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-MAP)
- [`BayesRTMB::RTMB_Fit_Base$fa_rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-fa_rotate)
- [`BayesRTMB::RTMB_Fit_Base$rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-rotate)

------------------------------------------------------------------------

### Method `get_point_estimate()`

Get point estimate for a target parameter (internal use).

#### Usage

    MCMC_Fit$get_point_estimate(target)

#### Arguments

- `target`:

  Target parameter name.

#### Returns

Matrix or array of point estimate.

------------------------------------------------------------------------

### Method `new()`

Create a new \`MCMC_Fit\` object.

#### Usage

    MCMC_Fit$new(
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

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a brief summary of the fitted object.

#### Usage

    MCMC_Fit$print(...)

#### Arguments

- `...`:

  Additional arguments.

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### Method `draws()`

Extract posterior draws for selected parameters.

#### Usage

    MCMC_Fit$draws(
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

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Summarize posterior draws.

#### Usage

    MCMC_Fit$summary(
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

### Method `unconstrain_draws()`

Transform posterior draws to the unconstrained scale.

#### Usage

    MCMC_Fit$unconstrain_draws()

#### Returns

Posterior draws on the unconstrained scale.

------------------------------------------------------------------------

### Method `log_prob()`

Evaluate log-probability values.

#### Usage

    MCMC_Fit$log_prob(safe = FALSE)

#### Arguments

- `safe`:

  Logical; whether to wrap the evaluation in a tryCatch block. Default
  is FALSE for speed.

#### Returns

Numeric vector of log-probability values.

------------------------------------------------------------------------

### Method `bridgesampling()`

Estimate the marginal likelihood by bridge sampling.

#### Usage

    MCMC_Fit$bridgesampling(
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

### Method [`bayes_factor()`](https://norimune.github.io/BayesRTMB/reference/bayes_factor.md)

Calculate the Bayes Factor against a null model or another fit object.

#### Usage

    MCMC_Fit$bayes_factor(
      null_model,
      bs_method = "normal",
      error_threshold = 0.2,
      ...
    )

#### Arguments

- `null_model`:

  Either a character string specifying the null target (e.g., "rho ~
  uniform(-1, 1)") or another MCMC_Fit object.

- `bs_method`:

  Character; the method to use for bridge sampling ("normal" or
  "warp3"). Default is "normal".

- `error_threshold`:

  Numeric; threshold for the approximate error warning. Default is 0.2.

- `...`:

  Additional arguments passed to the sample() method when fitting a null
  model (e.g., `chains = 4`, `sampling = 4000`).

------------------------------------------------------------------------

### Method `transformed_draws()`

Compute transformed parameters from posterior draws.

#### Usage

    MCMC_Fit$transformed_draws(tran_fn = NULL)

#### Arguments

- `tran_fn`:

  A function for transformed parameters.

#### Returns

Transformed parameter draws.

------------------------------------------------------------------------

### Method `generated_quantities()`

Compute generated quantities from posterior draws.

#### Usage

    MCMC_Fit$generated_quantities(code)

#### Arguments

- `code`:

  An \`rtmb_code( ... )\` or \` ... \` block containing the logic to be
  calculated using posterior samples.

#### Returns

The \`MCMC_Fit\` object itself (invisibly). Results are appended to the
\`generate_fit\` field.

------------------------------------------------------------------------

### Method `resolve_switching()`

Resolve label switching in posterior draws.

#### Usage

    MCMC_Fit$resolve_switching(
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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MCMC_Fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
