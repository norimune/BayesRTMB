## Test environments

* local Windows x86_64-w64-mingw32, R 4.5.3
* win-builder R-devel x86_64-w64-mingw32, Windows Server 2022 x64: OK
* GitHub Actions macOS 15 arm64, R-devel: OK

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no known CRAN reverse dependencies.

## Submission notes

This release updates BayesRTMB to 0.2.2.

This is a maintenance release submitted in response to the CRAN team's
macOS arm64 check report for the previous 0.2.1 release. The factor-analysis
wrapper code has been updated to avoid an automatic-differentiation container
construction issue that could occur on the M1 macOS R-devel check platform.

Additional changes in this release include:

* response-time distributions, including exponentially modified normal and
  diffusion model log-density functions;
* improved `obs(...)` sampling syntax for distributions with multiple observed
  variables;
* recursive capture of helper functions and constants used by user-defined
  functions supplied through the `setup` block;
* improved generated-quantity/report handling and AD-compatible helper
  containers.
