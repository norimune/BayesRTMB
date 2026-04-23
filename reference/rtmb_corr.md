# Wrapper for estimating correlation matrix (multivariate normal distribution)

Estimates a correlation matrix (along with means and standard
deviations) assuming a multivariate normal distribution from observation
data. If there are 2 observed variables, it automatically switches to
estimate the scalar correlation coefficient (\`corr\`) directly.

## Usage

``` r
rtmb_corr(
  data,
  prior = list(lkj_eta = 1, mu_sd = 10, sigma_rate = 1),
  init = NULL,
  null = NULL
)
```

## Arguments

- data:

  Observation data frame or matrix (N x P).

- prior:

  List of hyperparameters for prior distributions. Default is
  `list(lkj_eta = 1.0, mu_sd = 10, sigma_rate = 1.0)`.

- init:

  List of initial values (optional).

- null:

  Target when creating a null model (e.g., `"corr"`). Optional.

## Value

An instance of the `RTMB_Model` class.

## Examples

``` r
inst/examples/ex_corr.R
#> Error: object 'inst' not found
```
