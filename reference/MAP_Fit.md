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

  Parameter vector on the unconstrained scale.

- `par`:

  Parameter list on the constrained scale.

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

## Methods

### Public methods

- [`MAP_Fit$get_point_estimate()`](#method-map_fit-get_point_estimate)

- [`MAP_Fit$EAP()`](#method-map_fit-EAP)

- [`MAP_Fit$MAP()`](#method-map_fit-MAP)

- [`MAP_Fit$new()`](#method-map_fit-new)

- [`MAP_Fit$summary()`](#method-map_fit-summary)

- [`MAP_Fit$print()`](#method-map_fit-print)

- [`MAP_Fit$generated_quantities()`](#method-map_fit-generated_quantities)

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

Return point estimates (EAP is not applicable).

#### Usage

    MAP_Fit$EAP(...)

#### Arguments

- `...`:

  Ignored.

#### Returns

A named list of point estimates.

------------------------------------------------------------------------

### Method `MAP()`

Return point estimates (MAP sampling method is not applicable).

#### Usage

    MAP_Fit$MAP(...)

#### Arguments

- `...`:

  Ignored.

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
      se_samples = NULL
    )

#### Arguments

- `model`:

  The \`RTMB_Model\` object used for estimation.

- `par_vec`:

  Parameter vector on the unconstrained scale.

- `par`:

  Parameter list on the constrained scale.

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

------------------------------------------------------------------------

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Summarize MAP estimates.

#### Usage

    MAP_Fit$summary(pars = NULL, max_rows = 10, digits = 5)

#### Arguments

- `pars`:

  Character vector specifying the names of parameters to summarize. If
  NULL, all available parameters are summarized.

- `max_rows`:

  Maximum number of rows to print in summaries. Default is 10.

- `digits`:

  Number of digits to print.

#### Returns

A summary object, typically a data frame.

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

  An \`rtmb_code( ... )\` or \` ... \` block containing the logic to be
  calculated using the MAP estimate.

#### Returns

The \`MAP_Fit\` object itself (invisibly). Results are added or updated
in the \`generate\` list and \`df_generate\`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    MAP_Fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
