## R CMD check results

0 errors | 0 warnings | 1 note

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  New submission
  Possibly misspelled words in DESCRIPTION:
    Gelman
    Kristensen
    al
    et

This is a new submission. The possibly misspelled words are author names
(Gelman, Kristensen) and the abbreviation "et al." used in references.

## Resubmission

This is a resubmission. I have addressed the reviewer comments by:

* writing 'RTMB' in single quotes in the DESCRIPTION title;
* adding method references to the DESCRIPTION field in CRAN format;
* replacing \dontrun{} examples with \donttest{}.

Since the original submission on 2026-05-20, the package also includes
small bug fixes and documentation updates, including an AD-compatible
negative-binomial log-density implementation and a correction for the
unequal-variance JZS t-test example/documentation.
