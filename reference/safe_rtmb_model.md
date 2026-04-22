# Safe RTMB model construction (with error message translation)

Wraps \`rtmb_model\` to translate cryptic error messages generated
during \`MakeADFun\` execution (originating from C++/RTMB) into
user-friendly hints.

## Usage

``` r
safe_rtmb_model(
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

  A list of data that has passed \`validate_data\`.

- code:

  Code blocks for likelihood and priors defined with \`model_code()\`.

- par_names:

  A list of variable names corresponding to the dimensions of each
  parameter (optional).

- init:

  A list or numeric vector of initial values (optional).

- view:

  Character vector of parameter names to be displayed preferentially in
  summary outputs (optional).

- null_target:

  A character string specifying the parameter to fix and the prior to
  disable for creating a null model (optional).

## Value

Returns an \`rtmb_model\` object (R6 class instance) upon successful
compilation.
