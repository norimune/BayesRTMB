# Plot parameter estimates and credible intervals (Forest Plot)

Plot parameter estimates and credible intervals (Forest Plot)

## Usage

``` r
plot_forest(
  x,
  prob = 0.95,
  point_estimate = c("median", "EAP", "mean", "MAP", "marginal_map", "joint_map")
)
```

## Arguments

- x:

  A 3D array of posterior samples with dimensions \`(iterations, chains,
  variables)\`.

- prob:

  Numeric. Probability mass for the credible interval (default is 0.95).

- point_estimate:

  Character. Point estimate shown as dots. One of \`"median"\`
  (default), \`"EAP"\`/\`"mean"\`, \`"MAP"\`/\`"marginal_map"\`, or
  \`"joint_map"\`. \`"joint_map"\` uses the draw with the largest \`lp\`
  variable.

## Value

No return value.
