## Test environments

* local Windows x86_64-w64-mingw32, R 4.5.3:
  0 errors, 0 warnings, 0 notes
* win-builder R-release x86_64-w64-mingw32, R 4.6.1:
  Status: OK
* win-builder R-devel x86_64-w64-mingw32,
  R-devel 2026-07-23 (r90295): Status: OK

## Reverse dependencies

There are no known CRAN reverse dependencies.

## Submission notes

This is an update from BayesRTMB 0.2.3 to 0.2.4.

Changes since 0.2.3:

* allowed regression wrappers to resolve variables from the formula environment
  when `data` is omitted;
* added the Jeffreys scale prior to JZS t-tests;
* corrected exponential-prior rates and the default IRT discrimination prior;
* improved weak-prior calibration for mixture and latent-rank models;
* improved multivariate normal log-density performance; and
* updated English and Japanese vignettes.
