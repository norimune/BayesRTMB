# Prepare GLMM Formula Components

Builds the response, fixed-effect model matrix, and lme4-style
random-effect design components used by \`rtmb_glmer()\`. This is the
user-facing helper that reproduces what the wrapper places in the
generated \`setup\` block.

## Usage

``` r
make_glmer_re_terms(
  formula,
  data,
  family = "gaussian",
  resid_group = NULL,
  resid_time = NULL,
  within = NULL,
  factors = NULL,
  missing = "listwise"
)
```

## Arguments

- formula:

  lme4-style formula.

- data:

  Data frame.

- family:

  Character string of the distribution family.

- resid_group:

  Optional residual-correlation grouping variable.

- resid_time:

  Optional residual-correlation time variable.

- within:

  Optional list for wide-to-long conversion when the response is written
  as \`cbind(...)\`.

- factors:

  Optional character vector of variables to treat as factors.

## Value

A list containing \`Y\`, \`X\`, \`trials\`, \`offset\`, \`N\`,
fixed-effect metadata, and random-effect terms.
