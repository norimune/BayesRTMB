# RTMB-based IRT (Item Response Theory) Wrapper

Performs Item Response Theory modeling (1PL, 2PL, 3PL, and Graded
Response Model) using RTMB. Missing values (NA) in the data are
automatically removed internally for efficient computation.

## Usage

``` r
rtmb_irt(
  data,
  model = c("2PL", "1PL", "3PL"),
  type = c("binary", "ordered"),
  prior = prior_uniform(),
  init = NULL,
  fixed = NULL,
  view = NULL
)
```

## Arguments

- data:

  A data frame or matrix of item responses (N persons x J items).

- model:

  Character string for the model type: "1PL", "2PL", or "3PL".

- type:

  Character string for the data type: "binary" or "ordered".

- prior:

  Prior configuration: \`prior_uniform()\` (default) or
  \`prior_weak()\`. Hyperparameters can be specified within these
  functions (e.g., \`prior_weak(b_sd = 5)\`). Available parameters for
  IRT: \`a_rate\` (discrimination), \`b_sd\` (difficulty),
  \`c_alpha\`/\`c_beta\` (guessing).

- init:

  List of initial values.

- fixed:

  A named list of parameter values to fix (optional).

- view:

  Character vector of parameter names to prioritize in summary.
