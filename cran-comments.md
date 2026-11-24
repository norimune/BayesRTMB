## Test environments

* local Windows 11 x86_64-w64-mingw32, R 4.5.3:
  0 errors, 0 warnings, 0 notes

## Reverse dependencies

There are no known CRAN reverse dependencies.

## Submission notes

This is a feature and bug-fix update from BayesRTMB 0.4.0 to 0.5.0.

Changes since 0.4.0:

* added AD-compatible log-CDF and log-CCDF functions for supported
  distributions;
* made repeated generated-quantity calculations replace variables with the
  same name instead of duplicating them;
* preserved fitted contrast coding in conditional and simple effects; and
* corrected the documented dimensions of the `BigFive` data set.
