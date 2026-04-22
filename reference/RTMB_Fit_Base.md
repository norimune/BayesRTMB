# Base class for RTMB Fit objects

Base class for RTMB Fit objects

Base class for RTMB Fit objects

## Details

An R6 base class providing common methods for Bayesian and MAP inference
results.

## Public fields

- `model`:

  The \`RTMB_Model\` object used for estimation.

## Methods

### Public methods

- [`RTMB_Fit_Base$get_point_estimate()`](#method-RTMB_Fit_Base-get_point_estimate)

- [`RTMB_Fit_Base$EAP()`](#method-RTMB_Fit_Base-EAP)

- [`RTMB_Fit_Base$MAP()`](#method-RTMB_Fit_Base-MAP)

- [`RTMB_Fit_Base$rotate()`](#method-RTMB_Fit_Base-rotate)

- [`RTMB_Fit_Base$fa_rotate()`](#method-RTMB_Fit_Base-fa_rotate)

- [`RTMB_Fit_Base$clone()`](#method-RTMB_Fit_Base-clone)

------------------------------------------------------------------------

### Method `get_point_estimate()`

Abstract method to get a point estimate for a target parameter.

#### Usage

    RTMB_Fit_Base$get_point_estimate(target)

#### Arguments

- `target`:

  Character string specifying the target parameter name.

#### Returns

A numeric array or matrix of the point estimate.

------------------------------------------------------------------------

### Method `EAP()`

Calculate Expected A Posteriori (EAP) estimates from posterior samples.

#### Usage

    RTMB_Fit_Base$EAP(pars = "parameters", chains = NULL, best_chains = NULL)

#### Arguments

- `pars`:

  Character vector specifying the names of parameters to extract. Use
  "parameters" for only model parameters, "all" for all variables
  including transformed and generated quantities, or a character vector
  of specific variable names. Default is "parameters".

- `chains`:

  Numeric vector specifying the chains to use. Default is NULL (all
  chains).

- `best_chains`:

  Integer; number of best chains to retain based on mean log-posterior
  (lp) or ELBO. Default is NULL.

#### Returns

A named list of EAP estimates structured for use as \`init\`.

------------------------------------------------------------------------

### Method `MAP()`

Calculate Maximum A Posteriori (MAP) estimates from posterior samples.

#### Usage

    RTMB_Fit_Base$MAP(
      pars = "parameters",
      chains = NULL,
      best_chains = NULL,
      type = c("marginal", "joint")
    )

#### Arguments

- `pars`:

  Character vector specifying the names of parameters to extract. Use
  "parameters" for only model parameters, "all" for all variables
  including transformed and generated quantities, or a character vector
  of specific variable names. Default is "parameters".

- `chains`:

  Numeric vector specifying the chains to use. Default is NULL (all
  chains).

- `best_chains`:

  Integer; number of best chains to retain based on mean log-posterior
  (lp) or ELBO. Default is NULL.

- `type`:

  Character string specifying the type of MAP estimate. "marginal"
  (default) calculates the peak of the marginal posterior density for
  each parameter. "joint" returns the parameter values from the
  iteration with the highest log-posterior.

#### Returns

A named list of MAP estimates structured for use as \`init\`.

------------------------------------------------------------------------

### Method `rotate()`

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

### Method `fa_rotate()`

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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    RTMB_Fit_Base$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
