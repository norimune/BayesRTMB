# Classic Model Class for Frequentist Estimation

An R6 class representing a classical (frequentist) statistical model.
This class is used when \`classic = TRUE\` is specified in wrapper
functions like \`rtmb_lm\`.

## Public fields

- `type`:

  Character string specifying the model type (e.g., "lm").

- `formula`:

  The formula used for the model.

- `data`:

  The data used for estimation.

- `family`:

  The distribution family.

- `view`:

  Parameter names to prioritize in summary display.

- `obj`:

  Optional underlying model object (e.g., RTMB_Model for lmer).

- `refit_fn`:

  Optional function to re-fit the model on new data.

- `extra`:

  Optional list of additional parameters for specific model types.

## Methods

### Public methods

- [`Classic_Model$new()`](#method-Classic_Model-new)

- [`Classic_Model$estimate()`](#method-Classic_Model-estimate)

- [`Classic_Model$.resample_data()`](#method-Classic_Model-.resample_data)

- [`Classic_Model$.perform_fit()`](#method-Classic_Model-.perform_fit)

- [`Classic_Model$clone()`](#method-Classic_Model-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \`Classic_Model\` object.

#### Usage

    Classic_Model$new(
      type,
      formula,
      data,
      family = "gaussian",
      view = NULL,
      obj = NULL,
      refit_fn = NULL,
      extra = list()
    )

#### Arguments

- `type`:

  Character string specifying the model type (e.g., "lm", "lmer").

- `formula`:

  The formula used for the model.

- `data`:

  The data used for estimation.

- `family`:

  The distribution family.

- `view`:

  Parameter names to prioritize in summary display.

- `obj`:

  Optional underlying model object (e.g., RTMB_Model for lmer).

- `refit_fn`:

  Optional function to re-fit the model on new data.

- `extra`:

  Optional list of additional parameters for specific model types.

------------------------------------------------------------------------

### Method `estimate()`

Perform estimation using classical methods.

#### Usage

    Classic_Model$estimate(bootstrap = FALSE, n_boot = 1000)

#### Arguments

- `bootstrap`:

  Logical; whether to perform non-parametric bootstrap.

- `n_boot`:

  Integer; number of bootstrap samples.

#### Returns

A \`Classic_Fit\` object containing the results.

------------------------------------------------------------------------

### Method `.resample_data()`

(Internal) Resample data for bootstrap.

#### Usage

    Classic_Model$.resample_data()

#### Returns

Resampled data frame.

------------------------------------------------------------------------

### Method `.perform_fit()`

(Internal) Perform a single fit.

#### Usage

    Classic_Model$.perform_fit(data)

#### Arguments

- `data`:

  The data to fit.

#### Returns

A list containing the fit result and covariance matrix.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Classic_Model$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
