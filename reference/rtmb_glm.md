# RTMB-based GLM wrapper function (no random effects)

RTMB-based GLM wrapper function (no random effects)

## Usage

``` r
rtmb_glm(
  formula,
  data,
  family = "gaussian",
  penalty = c("none", "rhs", "ssp"),
  y_range = NULL,
  use_weak_info = FALSE,
  prior = list(Intercept_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, nu_rate = 0.1,
    cutpoint_sd = 2.5, shape_rate = 1, phi_rate = 1, lkj_eta = 1),
  weak_info_prior = list(max_beta = 1, sd_ratio = 0.5, expected_vars = 3, slab_scale = 2,
    slab_df = 4, ssp_ratio = 0.25),
  init = NULL,
  null = NULL
)
```

## Arguments

- formula:

  Formula

- data:

  Data frame

- family:

  Character string of the distribution family (e.g., "gaussian",
  "binomial", "poisson")

- penalty:

  Type of regularization for fixed effects: "none", "rhs" (Regularized
  Horseshoe), or "ssp" (Spike and Slab Prior). Default is "none".

- y_range:

  Theoretical minimum and maximum values of the response variable as a
  vector c(min, max). Specifying this automatically enables weakly
  informative priors.

- use_weak_info:

  Logical; whether to explicitly use weakly informative priors (requires
  y_range for continuous models).

- prior:

  List of hyperparameters for the default fixed priors.

- weak_info_prior:

  List of hyperparameters for the weakly informative priors and
  regularization.

- init:

  List of initial values

- null:

  Character string specifying the target parameter for the null model.
