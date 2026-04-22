# Calculate conditional effects for MCMC fit objects

Calculate conditional effects for MCMC fit objects

## Usage

``` r
# S3 method for class 'mcmc_fit'
conditional_effects(fit, effect, resolution = 100, prob = 0.95, ...)
```

## Arguments

- fit:

  An object of class \`MCMC_Fit\`.

- effect:

  Name of the explanatory variable to visualize (e.g., "X1" or "X1:X2").

- resolution:

  Grid resolution to calculate for continuous variables (default is
  100).

- prob:

  Probability for the credible interval (default is 0.95).

- ...:

  Additional arguments.
