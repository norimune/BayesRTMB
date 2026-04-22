# Plot autocorrelation for one variable across chains

Draw autocorrelation plots for a selected parameter across all chains.

## Usage

``` r
plot_acf(x, var_idx = 1)
```

## Arguments

- x:

  A 3D array of posterior samples with dimensions \`(iterations, chains,
  variables)\`.

- var_idx:

  Integer. Index of the variable to plot.

## Value

No return value. This function is called for its side effect of
plotting.
