# Plot pairs for posterior samples

Draw a scatterplot matrix to examine posterior correlations between
parameters.

## Usage

``` r
plot_pairs(x, pars = NULL)
```

## Arguments

- x:

  A 3D array of posterior samples with dimensions \`(iterations, chains,
  variables)\`.

- pars:

  Character vector or integer vector specifying which parameters to
  plot. If NULL (default), up to the first 10 parameters are plotted to
  prevent overplotting.

## Value

No return value.
