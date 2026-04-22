# Plot MCMC trace plots

Draw trace plots for each parameter across chains.

## Usage

``` r
plot_trace(x, mono = FALSE)
```

## Arguments

- x:

  A 3D array of posterior samples with dimensions \`(iterations, chains,
  variables)\`. A 2D matrix \`(iterations, chains)\` is also allowed and
  is treated as a single variable.

- mono:

  Logical. If \`TRUE\`, plot in monochrome. If \`FALSE\`, use blue
  shades.

## Value

No return value. This function is called for its side effect of
plotting.
