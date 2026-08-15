## Test environments

* local Windows 11 x86_64-w64-mingw32, R 4.5.3:
  0 errors, 0 warnings, 0 notes
* win-builder, R-release 4.6.1 x86_64-w64-mingw32:
  Status: OK
* win-builder, R-devel 2026-08-12 (r90394) x86_64-w64-mingw32:
  Status: OK

## Reverse dependencies

There are no known CRAN reverse dependencies.

## Submission notes

This is an update from BayesRTMB 0.2.4 to 0.3.0.

Changes since 0.2.4:

* added the setup-only `.data` binding and made wrapper-generated setup code
  reproduce the original data interface more directly;
* added standardized regression coefficients through `std = TRUE`;
* added `summary_mcmc()` for ordinary numeric MCMC arrays;
* added grand-mean and within-cluster centering helpers;
* improved generated wrapper code and `rtmb_table(..., data = ...)` handling;
* corrected Jacobian adjustments when constrained parameters are fixed; and
* updated the English and Japanese documentation.
