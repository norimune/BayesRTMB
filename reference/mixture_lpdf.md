# Generic Mixture log-probability density function

Generic Mixture log-probability density function

## Usage

``` r
mixture_lpdf(x, pi_w, lpdf_list, sum = TRUE)
```

## Arguments

- x:

  Vector of quantiles.

- pi_w:

  Vector or matrix of mixing proportions.

- lpdf_list:

  A list of log-density vectors (one for each cluster).

## Value

The sum of the log-density.
