# MAP fit object

MAP fit object

MAP fit object

## Details

An R6 class storing optimization results from maximum a posteriori (MAP)
estimation.

## Super class

[`BayesRTMB::RTMB_Fit_Base`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.md)
-\> `map_fit`

## Public fields

- `model`:

  The \`RTMB_Model\` object used for estimation.

- `par_vec`:

  Parameter vector on the unconstrained scale (constrained values
  unlisted).

- `par`:

  Parameter list on the constrained scale.

- `par_unc`:

  Parameter vector on the unconstrained scale (raw unconstrained
  values).

- `ci_method`:

  Method used for CI estimation ("wald", "profile", or "sampling").

- `objective`:

  RTMB objective function object.

- `log_ml`:

  Log marginal likelihood or related model criterion.

- `convergence`:

  Optimizer convergence code.

- `sd_rep`:

  Standard deviation report object.

- `df_fixed`:

  Summary table for fixed-effect parameters.

- `random_effects`:

  Random effect estimates.

- `df_transform`:

  Summary table for transformed parameter estimates.

- `df_generate`:

  Summary table for generated quantity estimates.

- `opt_history`:

  A vector of optimize objective history.

- `transform`:

  List of transformed parameters maintaining their original dimensions.

- `generate`:

  List of generated quantities maintaining their original dimensions.

- `se_samples`:

  List of simulated samples for standard error estimation.

- `laplace`:

  Logical; whether Laplace approximation was used.

- `vcov_unc`:

  Variance-covariance matrix of parameters in unconstrained space.

- `map`:

  List; the parameter mapping used.

- `marginal_vars`:

  Character vector of parameter names that were marginalized (integrated
  out).

## Methods

### Public methods

- [`MAP_Fit$get_point_estimate()`](#method-map_fit-get_point_estimate)

- [`MAP_Fit$EAP()`](#method-map_fit-EAP)

- [`MAP_Fit$MAP()`](#method-map_fit-MAP)

- [`MAP_Fit$new()`](#method-map_fit-new)

- [`MAP_Fit$ranef()`](#method-map_fit-ranef)

- [`MAP_Fit$draws()`](#method-map_fit-draws)

- [`MAP_Fit$summary()`](#method-map_fit-summary)

- [`MAP_Fit$savage_dickey()`](#method-map_fit-savage_dickey)

- [`MAP_Fit$print()`](#method-map_fit-print)

- [`MAP_Fit$generated_quantities()`](#method-map_fit-generated_quantities)

- [`MAP_Fit$profile()`](#method-map_fit-profile)

- [`MAP_Fit$clone()`](#method-map_fit-clone)

Inherited methods

- [`BayesRTMB::RTMB_Fit_Base$fa_rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-fa_rotate)
- [`BayesRTMB::RTMB_Fit_Base$rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-rotate)

------------------------------------------------------------------------

### Method `get_point_estimate()`

Get point estimate for a target parameter (internal use).

#### Usage

    MAP_Fit$get_point_estimate(target)

#### Arguments

- `target`:

  Target parameter name.

#### Returns

Matrix or array of point estimate.

------------------------------------------------------------------------

### Method `EAP()`

Return point estimates (EAP is not applicable for MAP).

#### Usage

    MAP_Fit$EAP(pars = NULL)

#### Arguments

- `pars`:

  Optional character vector of parameter names to extract.

#### Returns

A named list of point estimates.

------------------------------------------------------------------------

### Method `MAP()`

Return point estimates (MAP sampling method is not applicable).

#### Usage

    MAP_Fit$MAP(pars = NULL)

#### Arguments

- `pars`:

  Optional character vector of parameter names to extract.

#### Returns

A named list of point estimates.

------------------------------------------------------------------------

### Method `new()`

Create a new \`MAP_Fit\` object.

#### Usage

    MAP_Fit$new(
      model,
      par_vec,
      par,
      objective,
      log_ml,
      convergence,
      sd_rep,
      df_fixed,
      random_effects,
      df_transform = NULL,
      df_generate = NULL,
      opt_history = NULL,
      transform = NULL,
      generate = NULL,
      se_samples = NULL,
      par_unc = NULL,
      vcov_unc = NULL,
      ci_method = "wald",
      laplace = TRUE,
      map = NULL,
      marginal_vars = NULL
    )

#### Arguments

- `model`:

  The \`RTMB_Model\` object used for estimation.

- `par_vec`:

  Parameter vector on the unconstrained scale (constrained values
  unlisted).

- `par`:

  Parameter list on the constrained scale.

- `objective`:

  The objective function value at the optimum.

- `log_ml`:

  Log marginal likelihood.

- `convergence`:

  Optimizer convergence code.

- `sd_rep`:

  The \`sdreport\` object from TMB.

- `df_fixed`:

  Data frame of fixed effects estimates and CIs.

- `random_effects`:

  Data frame of random effects estimates and CIs.

- `df_transform`:

  Data frame of transformed parameters.

- `df_generate`:

  Data frame of generated quantities.

- `opt_history`:

  Data frame of optimization history.

- `transform`:

  List of transformed parameters maintaining their original dimensions.

- `generate`:

  List of generated quantities maintaining their original dimensions.

- `se_samples`:

  List of simulated samples for standard error estimation.

- `par_unc`:

  Parameter vector on the unconstrained scale (raw values).

- `vcov_unc`:

  Variance-covariance matrix of parameters in unconstrained space.

- `ci_method`:

  Method used for CI estimation ("wald" or "sampling").

- `laplace`:

  Logical; whether Laplace approximation was used.

- `map`:

  List; the parameter mapping used.

- `marginal_vars`:

  Character vector of parameter names that were marginalized (integrated
  out).

------------------------------------------------------------------------

### Method `ranef()`

Return random effect estimates as a named list.

#### Usage

    MAP_Fit$ranef()

#### Returns

A named list of random effect estimates.

------------------------------------------------------------------------

### Method `draws()`

Extract samples from the asymptotic posterior distribution.

#### Usage

    MAP_Fit$draws(
      pars = NULL,
      inc_random = FALSE,
      inc_transform = TRUE,
      inc_generate = TRUE,
      ...
    )

#### Arguments

- `pars`:

  Character or numeric vector specifying the names or indices of
  parameters to extract.

- `inc_random`:

  Logical; whether to include random effects.

- `inc_transform`:

  Logical; whether to include transformed parameters.

- `inc_generate`:

  Logical; whether to include generated quantities.

- `...`:

  Ignored.

#### Returns

An array of samples \[iterations, 1, parameters\].

------------------------------------------------------------------------

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Summarize MAP estimates.

#### Usage

    MAP_Fit$summary(pars = NULL, max_rows = 10, digits = 5, ranef = FALSE)

#### Arguments

- `pars`:

  Character vector specifying the names of parameters to summarize. If
  NULL, all available parameters are summarized.

- `max_rows`:

  Maximum number of rows to print in summaries. Default is 10.

- `digits`:

  Number of digits to print.

- `ranef`:

  Logical; whether to also display random effect estimates. Default is
  FALSE.

#### Returns

A summary object, typically a data frame.

------------------------------------------------------------------------

### Method `savage_dickey()`

Calculate Bayes Factors using the Savage-Dickey density ratio.

#### Usage

    MAP_Fit$savage_dickey(pars = NULL, null = 0, digits = 3)

#### Arguments

- `pars`:

  Optional character vector of parameter names to test. If NULL, tests
  all fixed effects.

- `null`:

  The null hypothesis value (on the constrained scale). Default is 0.

- `digits`:

  Number of decimal places in the output.

#### Returns

A data frame containing Bayes Factors and evidence descriptions.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a brief summary of the fitted object.

#### Usage

    MAP_Fit$print(pars = NULL, max_rows = 10, digits = 5, ...)

#### Arguments

- `pars`:

  Character vector specifying the names of parameters to summarize.

- `max_rows`:

  Maximum number of rows to print in summaries.

- `digits`:

  Number of digits to print.

- `...`:

  Additional arguments passed to the \`summary\` method.

#### Returns

The object itself, invisibly.

------------------------------------------------------------------------

### Method `generated_quantities()`

Compute generated quantities from the MAP estimate.

#### Usage

    MAP_Fit$generated_quantities(code)

#### Arguments

- `code`:

  An \`rtmb_code({ ... })\` or \`{ ... }\` block containing the logic to
  be calculated using the MAP estimate.

#### Returns

The \`MAP_Fit\` object itself (invisibly). Results are added or updated
in the \`generate\` list and \`df_generate\`.

------------------------------------------------------------------------

### Method [`profile()`](https://rdrr.io/r/stats/profile.html)

Calculate Profile Likelihood confidence intervals for specific
parameters.

#### Usage

    MAP_Fit$profile(
      pars = NULL,
      level = 0.95,
      trace = FALSE,
      digits = 5,
      show_plot = FALSE,
      quiet = FALSE,
      jacobian = "none",
      ...
    )

#### Arguments

- `pars`:

  Character vector of parameter names to profile. If NULL, all fixed
  parameters are profiled.

- `level`:

  Confidence level (default is 0.95).

- `trace`:

  Logical; whether to print profiling progress. Default is FALSE.

- `digits`:

  Integer; number of decimal places to print. Default is 5.

- `show_plot`:

  Logical; whether to plot the profile likelihood curves. Default is
  FALSE.

- `quiet`:

  Logical; whether to suppress text output. Default is FALSE.

- `jacobian`:

  Character; "none" (default), "random", or "all". Whether to include
  Jacobian adjustments for transformations.

- `...`:

  Additional arguments passed to TMB::tmbprofile (e.g., ytol).

#### Returns

A data frame containing the profile-based confidence intervals, with the
raw profile objects stored in the "profiles" attribute.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MAP_Fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
