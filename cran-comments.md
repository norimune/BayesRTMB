## R CMD check results

0 errors | 0 warnings | 0 notes

## check_win_devel() results

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Possibly misspelled words in DESCRIPTION:
    RTMB

RTMB is the name of an R package used as the automatic differentiation
engine by BayesRTMB.

## Resubmission

This is a resubmission. I have addressed the reviewer comments by:

* writing 'RTMB' in single quotes in the DESCRIPTION title;
* adding method references to the DESCRIPTION field in CRAN format;
* replacing \dontrun{} examples with \donttest{}.

Since the original submission on 2026-05-20, the package also includes
small bug fixes and documentation updates, including an AD-compatible
negative-binomial log-density implementation and a correction for the
unequal-variance JZS t-test example/documentation.
