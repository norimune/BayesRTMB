# RTMB-based Bayesian two-sample t-test wrapper function

Performs a Bayesian two-sample t-test using RTMB. It estimates the
effect size (delta) with a Cauchy prior, allowing for robust inference
and calculation of Bayes factors.

## Usage

``` r
rtmb_ttest(
  y1,
  y2,
  r = 0.707,
  y_range = NULL,
  use_weak_info = FALSE,
  prior = list(mean_sd = 10, sd_rate = 0.1),
  weak_info_prior = list(sd_ratio = 0.5),
  init = NULL,
  null = NULL
)
```

## Arguments

- y1:

  Numeric vector of responses for group 1.

- y2:

  Numeric vector of responses for group 2.

- r:

  Numeric; Cauchy prior scale for the effect size (delta). Default is
  0.707.

- y_range:

  Theoretical minimum and maximum values of the response variable as a
  vector c(min, max). Specifying this automatically enables weakly
  informative priors.

- use_weak_info:

  Logical; whether to explicitly use weakly informative priors.

- prior:

  List of hyperparameters for the default fixed priors.

- weak_info_prior:

  List of hyperparameters for the weakly informative priors.

- init:

  List of initial values.

- null:

  Character string specifying the target parameter for the null model
  (e.g., "delta" or "delta ~ cauchy(0, r)").

## Value

An `RTMB_Model` object.

## Examples

``` r
inst/examples/ex_ttest.R
#> Error: object 'inst' not found
```
