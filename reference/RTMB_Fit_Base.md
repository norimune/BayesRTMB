# Base class for RTMB Fit objects

An R6 base class providing common methods for Bayesian and MAP inference
results.

## Public fields

- `model`:

  The \`RTMB_Model\` object used for estimation.

## Methods

### Public methods

- [`RTMB_Fit_Base$get_point_estimate()`](#method-RTMB_Fit_Base-get_point_estimate)

- [`RTMB_Fit_Base$estimate()`](#method-RTMB_Fit_Base-estimate)

- [`RTMB_Fit_Base$EAP()`](#method-RTMB_Fit_Base-EAP)

- [`RTMB_Fit_Base$MAP()`](#method-RTMB_Fit_Base-MAP)

- [`RTMB_Fit_Base$rotate()`](#method-RTMB_Fit_Base-rotate)

- [`RTMB_Fit_Base$fa_rotate()`](#method-RTMB_Fit_Base-fa_rotate)

- [`RTMB_Fit_Base$clone()`](#method-RTMB_Fit_Base-clone)

------------------------------------------------------------------------

### `RTMB_Fit_Base$get_point_estimate()`

Abstract method to get a point estimate for a target parameter.

#### Usage

    RTMB_Fit_Base$get_point_estimate(target, ...)

#### Arguments

- `target`:

  Character string specifying the target parameter name.

- `...`:

  Additional arguments.

#### Returns

A numeric array or matrix of the point estimate.

------------------------------------------------------------------------

### `RTMB_Fit_Base$estimate()`

Get point estimates for parameters, transformed parameters, and
generated quantities.

#### Usage

    RTMB_Fit_Base$estimate(
      pars = NULL,
      type = c("mean", "EAP", "marginal_map", "joint_map", "MAP"),
      component = c("all", "parameters", "transform", "generate"),
      chains = NULL,
      best_chains = NULL,
      drop = TRUE,
      ...
    )

#### Arguments

- `pars`:

  Optional character or numeric vector of parameter names or indices to
  extract. Supports special keywords: "parameters", "transform",
  "generate", and "all".

- `type`:

  Character string specifying the estimation type.

- `component`:

  Character string specifying the component to filter by.

- `chains`:

  Numeric vector of chains to include.

- `best_chains`:

  Number of best chains to include.

- `drop`:

  Logical; if TRUE and only one parameter is selected, return the value
  directly instead of a list.

- `...`:

  Additional arguments passed to draws().

#### Returns

A named list of point estimates, or a single value if \`drop = TRUE\`.

------------------------------------------------------------------------

### `RTMB_Fit_Base$EAP()`

Calculate Expected A Posteriori (EAP) estimates from posterior samples.

#### Usage

    RTMB_Fit_Base$EAP(
      pars = "parameters",
      chains = NULL,
      best_chains = NULL,
      drop = FALSE,
      ...
    )

#### Arguments

- `pars`:

  Optional character vector of parameter names to extract.

- `chains`:

  Numeric vector of chains to include.

- `best_chains`:

  Number of best chains to include.

- `drop`:

  Logical; whether to drop the list if only one parameter is selected.

- `...`:

  Additional arguments passed to \`estimate()\`.

#### Returns

A named list of EAP estimates.

------------------------------------------------------------------------

### `RTMB_Fit_Base$MAP()`

Calculate Maximum A Posteriori (MAP) estimates.

#### Usage

    RTMB_Fit_Base$MAP(
      pars = "parameters",
      chains = NULL,
      best_chains = NULL,
      type = c("marginal", "joint"),
      drop = FALSE,
      ...
    )

#### Arguments

- `pars`:

  Optional character vector of parameter names to extract.

- `chains`:

  Numeric vector of chains to include.

- `best_chains`:

  Number of best chains to include.

- `type`:

  Character string; "marginal" or "joint" MAP.

- `drop`:

  Logical; whether to drop the list if only one parameter is selected.

- `...`:

  Additional arguments passed to \`estimate()\`.

#### Returns

A named list of MAP estimates.

------------------------------------------------------------------------

### `RTMB_Fit_Base$rotate()`

Rotate sampled parameters.

#### Usage

    RTMB_Fit_Base$rotate(target, reference = NULL, linked = NULL)

#### Arguments

- `target`:

  Character string specifying the target variable to base the rotation
  on.

- `reference`:

  Matrix to rotate towards. If NULL, the target's point estimate is
  used.

- `linked`:

  Character vector of variable names to be rotated in the same
  direction.

- `overwrite`:

  Logical; whether to overwrite the stored draws. If FALSE, adds to
  generated quantities. Default is FALSE.

- `suffix`:

  Character string to append to the rotated variable names when
  overwrite is FALSE. Default is "rot".

- `...`:

  Additional arguments passed to the rotation function.

#### Returns

The updated object invisibly.

------------------------------------------------------------------------

### `RTMB_Fit_Base$fa_rotate()`

Rotate factor loadings and optional factor scores.

#### Usage

    RTMB_Fit_Base$fa_rotate(
      target = "L",
      linked = NULL,
      scores = NULL,
      rotate = "promax",
      ...
    )

#### Arguments

- `target`:

  Character string specifying the target variable to base the rotation
  on.

- `linked`:

  Character vector of variable names to be rotated in the same
  direction.

- `scores`:

  Character vector of variable names to be rotated as factor scores
  (inverse direction).

- `rotate`:

  Character string specifying the rotation method.

- `...`:

  Additional arguments passed to the rotation function.

#### Returns

The updated object invisibly.

------------------------------------------------------------------------

### `RTMB_Fit_Base$clone()`

The objects of this class are cloneable with this method.

#### Usage

    RTMB_Fit_Base$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
