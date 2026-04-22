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
  prior = list(a_log_mean = 0, a_log_sd = 0.5, b_mean = 0, b_sd = 2.5, c_alpha = 1,
    c_beta = 4, theta_sd = 1),
  init = NULL
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

  List of hyperparameters for prior distributions.

- init:

  List of initial values.
