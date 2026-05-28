# BayesRTMB 0.1.1

* Addressed CRAN resubmission feedback by quoting 'RTMB' in the package title,
  adding method references to DESCRIPTION, and replacing `\dontrun{}` examples
  with `\donttest{}` where appropriate.
* Trimmed long-running examples for CRAN checks while retaining representative
  MCMC examples for correlation, t-test, and mixed-model workflows.
* Updated IRT post-estimation examples to use ordered response data explicitly.
* Fixed an AD-compatible negative-binomial log-density issue.
* Corrected the unequal-variance JZS t-test example and documentation.
