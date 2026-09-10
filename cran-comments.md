## Test environments

* local Windows 11 x86_64-w64-mingw32, R 4.6.0:
  0 errors, 0 warnings, 0 notes

## Reverse dependencies

There are no known CRAN reverse dependencies.

## Submission notes

This is a feature update from BayesRTMB 0.3.0 to 0.4.0.

Changes since 0.3.0:

* added `MCMC_Fit$continue_sampling()` to extend existing NUTS chains while
  reusing their adapted sampler state;
* corrected fixed-effect standard errors for Gaussian mixed models fitted with
  `classic()` and REML;
* corrected marginal-likelihood normalization for correlation models using a
  flat prior;
* fixed simulation-based standard errors for one-parameter models;
* improved error reporting for posterior prediction from `classic()` fits;
* changed the fallback distance used by `plot_mdu()` to Euclidean; and
* updated the related tests and documentation.
