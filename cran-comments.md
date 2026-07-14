## Test environments

* local Windows x86_64-w64-mingw32, R 4.5.3: `R CMD check` completed with
  0 errors, 0 warnings, and 0 notes
* win-builder R-devel x86_64-w64-mingw32, Windows Server 2022 x64,
  R-devel 2026-07-12 (r90242): OK
* external Apple M1 Mac mini, macOS 26.5.1, R 4.6.1 aarch64,
  Apple clang 21, RTMB 1.9, TMB 1.9.21: OK
  * source tarball build, tarball installation, `example(rtmb_fa)`, the
    targeted `rtmb_fa(data = fa_data, nfactors = 2, score = TRUE)`
    regression, post-hoc Promax rotation, and `R CMD check --as-cran
    --no-manual` all passed
* external Apple M1 Mac tests of the previous 0.2.2 source with RTMB 1.9:
  R 4.6.1 release and R 4.5.3 oldrel both passed build, installation,
  targeted `rtmb_fa()` execution, and `R CMD check --as-cran --run-donttest`
* external Apple Silicon R-devel check configured to closely match the BDR
  M1mac setup (macOS 26.5.x, Apple clang 21, CRAN gfortran 14.2 fork,
  stack size 20M, Europe/London timezone, en_GB.UTF-8 with LC_COLLATE=C,
  and the relevant `_R_CHECK_*` environment variables): the source containing
  the factor-analysis AD fixes passed `R CMD check --as-cran --run-donttest`
  without package-caused ERRORs or WARNINGs
* GitHub Actions macOS arm64, R-devel, commit bafad68:
  * macos-15: OK
  * macos-26: OK
  * source tarball build, tarball installation, `example(rtmb_fa)`, the
    targeted `rtmb_fa(data = fa_data, nfactors = 2, score = TRUE)`
    regression, post-hoc Promax rotation, and `R CMD check --as-cran
    --no-manual` all passed

## R CMD check results

win-builder R-devel (2026-07-12, r90242): OK.

GitHub Actions macOS arm64 checks passed on macos-15 and macos-26. External
Apple M1 checks also passed for the updated source tarball, and the originally
reported failure was not reproducible on an independent Apple M1 Mac using
R 4.6.1 release or R 4.5.3 oldrel with RTMB 1.9.
The failure was also not reproduced in an external Apple Silicon R-devel
environment configured to closely match the BDR M1mac setup; the remaining
known difference was the exact R-devel revision available on that machine.
We also attempted to use the macOS builder service for the CRAN M1mac setup,
but the service was unavailable during preparation for this resubmission
("The macOS package build services are currently not responding."). As an
alternative, we checked the package on GitHub Actions macOS arm64 and on
independent Apple M1 hardware.

## Reverse dependencies

There are no known CRAN reverse dependencies.

## Submission notes

This is a resubmission as BayesRTMB 0.2.3 with an increased version number.

In this resubmission, the `rtmb_fa()` example failure reported by CRAN on
the M1mac check platform has been addressed. The reported failure occurred in
the factor-analysis transform block. The failing code path combined
base-matrix construction and recycling involving automatic-differentiation
values with accesses to structural-zero upper-triangular entries. On the CRAN
M1mac platform this could result in an invalid `advector` during subsequent AD
expressions.

The correction has two parts:

* constrained matrix parameters are now initialized with `rtmb_array()` rather
  than base `matrix(...)` construction when AD values are involved, preserving
  an AD-compatible container;
* the standard `rtmb_fa()` transform and generated quantities now compute
  lower-triangular loading summaries by looping only over `k <= min(j, K)`,
  avoiding reads from structural-zero upper-triangular entries and avoiding
  base-matrix reconstruction of AD values.

The runnable `rtmb_fa()` example has also been simplified to a one-factor
example only. More advanced factor-analysis workflows, including factor scores,
rotations, and SSP regularization, remain documented outside the runnable Rd
example and are checked in CI regression tests rather than during CRAN example
execution.

This is a maintenance release prepared in response to the CRAN team's
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
