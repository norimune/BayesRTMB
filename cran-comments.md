## Test environments

* local Windows 10 x86_64-w64-mingw32, R 4.5.3
* win-builder R-devel x86_64-w64-mingw32, Windows Server 2022 x64,
  R Under development (2026-06-21 r90185 ucrt): OK

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no known CRAN reverse dependencies.

## Submission notes

This is a resubmission.

In the previous submission, CRAN reported that the following README URLs
returned 404:

* https://norimune.github.io/BayesRTMB/articles/analysis_reference.html
* https://norimune.github.io/BayesRTMB/articles/ja-analysis_reference.html

These URLs were invalid because the pkgdown GitHub Pages deployment had failed.
The pkgdown configuration has been fixed, the site has been redeployed, and both
URLs now return HTTP 200 OK.

This release updates BayesRTMB to 0.2.1.

The update includes performance and robustness improvements for MCMC,
variational inference, transformed parameters, generated quantities, and
wrapper model code. It also improves AD-compatible helper containers and
softmax/log-sum-exp handling used in model code.
