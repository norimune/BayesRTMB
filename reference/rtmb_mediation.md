# RTMB-based Mediation Analysis Wrapper

\`rtmb_mediation\` performs mediation analysis by simultaneously
estimating multiple GLM regression equations. It automatically
identifies mediation paths and calculates indirect, direct, and total
effects.

## Usage

``` r
rtmb_mediation(
  formula,
  data,
  family = "gaussian",
  prior = prior_flat(),
  y_range = NULL,
  fixed = NULL,
  view = NULL,
  ...
)
```

## Arguments

- formula:

  A list of formulas defining the regression paths (e.g., \`list(M ~ X,
  Y ~ X + M)\`).

- data:

  A data frame containing the variables.

- family:

  A single character string or a list of character strings specifying
  the error distribution for each equation (e.g., \`family =
  list("gaussian", "binomial")\`). Default is "gaussian".

- prior:

  An object of class "rtmb_prior" specifying the prior distribution.

- y_range:

  Theoretical minimum and maximum values of the response variable.

- fixed:

  A named list of parameter values to fix (optional).

- view:

  Character vector of parameter names to prioritize in summary.

- ...:

  Additional arguments passed to the model construction.

## Value

An \`RTMB_Model\` object.

## Details

The function identifies mediation paths by looking for variables that
are responses in one equation and predictors in another. Indirect
effects are calculated as the product of coefficients along these paths
(\\a \* b\\).

**Uncertainty Estimation**: When using \`\$optimize(ci_method =
"sampling")\`, the function provides asymmetric confidence intervals for
indirect effects based on the distribution of products, which is more
accurate than the standard Sobel test (Delta Method).
