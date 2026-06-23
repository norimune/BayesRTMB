# Simple effects for MCMC fit objects

Simple effects for MCMC fit objects

## Usage

``` r
# S3 method for class 'mcmc_fit'
simple_effects(
  fit,
  effect,
  prob = 0.95,
  sd_multiplier = 1,
  sd_slice = NULL,
  ...
)
```

## Arguments

- fit:

  An object of class \`MCMC_Fit\`.

- effect:

  Interaction term (e.g., "A:B").

- prob:

  Probability for credible intervals.

- sd_multiplier:

  Multiplier for SD for continuous moderators.

- sd_slice:

  Logical or NULL; controls whether continuous moderators are evaluated
  at mean - SD, mean, and mean + SD.

- ...:

  Additional arguments.
