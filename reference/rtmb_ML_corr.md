# Multilevel Correlation Analysis (Kenny's Model)

This function fits a multilevel correlation model to estimate
between-group and within-group covariance/correlation structures, along
with Intraclass Correlation Coefficients (ICC).

## Usage

``` r
rtmb_ML_corr(
  formula,
  ID,
  data = NULL,
  prior_type = c("weakly_informative", "uniform"),
  ...
)
```

## Arguments

- formula:

  A formula or a character vector of variable names.

- ID:

  A character string or expression specifying the group ID variable.

- data:

  A data frame.

- prior_type:

  Prior type: "weakly_informative" or "uniform".

- ...:

  Additional arguments passed to `rtmb_model`.

## Value

A `RTMB_Model` object.
