# RTMB-based Factor Analysis Wrapper

Performs exploratory factor analysis (EFA) using RTMB. It supports
standard rotation methods (e.g., varimax, promax) as well as regularized
factor analysis using a Spike-and-Slab Prior (SSP) for estimating sparse
loading matrices.

## Usage

``` r
rtmb_fa(
  data,
  nfactors = 1,
  rotate = NULL,
  score = FALSE,
  prior = list(mean_sd = 10, loadings_sd = 1, sd_rate = 10, ssp_ratio = 0.25),
  init = NULL
)
```

## Arguments

- data:

  Observation data frame or matrix (N x J).

- nfactors:

  Number of factors (K).

- rotate:

  String specifying the rotation method (e.g., "varimax", "promax",
  "ssp"). If NULL, no rotation is applied. Specifying "ssp" performs
  regularized factor analysis.

- score:

  Logical; if TRUE, factor scores are calculated in the generate block
  (default is FALSE).

- prior:

  List of hyperparameters for prior distributions. \`ssp_ratio\`
  represents the proportion of non-zero loadings per factor when "ssp"
  is specified.

- init:

  List of initial values. If not provided, initial values are
  automatically generated based on PCA or the psych package.
