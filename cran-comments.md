## Test environments

* local Windows x86_64-w64-mingw32, R 4.5.3: targeted `rtmb_fa()`
  examples OK
* win-builder R-devel x86_64-w64-mingw32, Windows Server 2022 x64: to be
  rerun before resubmission
* GitHub Actions macOS arm64, R-devel: to be rerun before resubmission

## R CMD check results

To be updated before resubmission.

## Reverse dependencies

There are no known CRAN reverse dependencies.

## Submission notes

This is a resubmission as BayesRTMB 0.2.3 with an increased version number.

In this resubmission, the `rtmb_fa()` example failure reported by CRAN on
the M1mac check platform has been addressed. The failure was narrowed to
factor-analysis code that used lower-triangular AD loading matrices as ordinary
dense matrices in transformed quantities. On the CRAN M1mac platform this could
touch structural-zero entries above the triangular part and produce an invalid
`advector`.

The correction has two parts:

* constrained matrix parameters are now initialized with `rtmb_array()` rather
  than `matrix(ad_zero, ...)`, so structural zeros are represented by a stable
  AD-compatible container;
* the standard `rtmb_fa()` transform and generated quantities now compute
  lower-triangular loading summaries by looping only over `k <= min(j, K)`,
  avoiding reads from the structural-zero upper-triangular entries.

The runnable `rtmb_fa()` example now demonstrates Promax rotation via post-hoc
`fa_rotate()` after model fitting, rather than including the rotation in the
model-construction step.

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
