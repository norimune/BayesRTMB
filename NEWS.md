# BayesRTMB 0.2.1

* Added `upgrade_fit()` to rebuild saved MCMC, VB, MAP, and classic fit objects
  with the currently loaded class definitions, optionally upgrading their
  embedded model objects as well.
* Improved `rtmb_vector()` and `rtmb_array()` tape construction time by
  automatically reusing an AD seed from model parameters when available.
* Made `log_sum_exp()`, `softmax()`, and `log_softmax()` work more reliably
  with RTMB automatic-differentiation values, including baseline-category
  patterns such as `softmax(c(0, eta))` inside `rtmb_code()`.
* Improved wrapper-generated model code to use AD-compatible `rtmb_vector()`
  and `rtmb_array()` containers in loop-filled generated quantities and
  generated likelihood contributions where needed.
* Improved `report()` handling in transformed and generated quantities,
  including namespaced `BayesRTMB::report()` calls and wrapper-generated
  `print_code()` output.
* Changed VB point estimates to use only the best variational estimate by
  default, aligning `EAP()`, `MAP()`, and rotation references with the selected
  best ELBO run while still allowing explicit `chains` or `best_chains`
  selection.
* Made `EAP()` and `MAP()` drop their list wrapper by default when a single
  parameter is requested, matching the behavior of `estimate()`.
* Added optional taped evaluation for transformed parameters and generated
  quantities, with automatic fallback to R evaluation when taping is not
  possible.
* Improved parallel worker robustness by reducing exported globals and
  preserving wrapper setup environments needed by generated model code.
* Improved MCMC runtime behavior by caching metric calculations, speeding up
  retained draw conversion, and refining progress checks.
* Changed bootstrap progress reporting to use percentage-style progress output,
  consistent with other long-running workflows.
* Added diagnostic recommendations to help interpret common fitting warnings.
* Improved matrix-valued Gaussian process log-density evaluation.
* Updated MDU defaults and internals, including Euclidean distance as the
  default MDU distance and more explicit use of namespace-qualified factor
  rotations.

# BayesRTMB 0.2.0

* Reworked NUTS sampling internals with Stan-style multinomial tree expansion,
  warmup diagnostics, Stan-window metric adaptation, and support for diagonal,
  dense, hybrid, and automatic metric selection.
* Improved MCMC diagnostics by reporting divergence counts and percentages,
  metric auto-selection details, warmup summaries, metric condition numbers, and
  positive-definite fallback counts.
* Added configurable progress output for MCMC and VB workflows, including
  streamed message-style progress for parallel workers and percentage reporting.
* Added delta-method standard errors and confidence intervals for
  `conditional_effects()` and `simple_effects()` with optimized and classic fits;
  `simple_effects()` for classic fits now also reports `df`, `t value`, and `Pr`.
* Added `sd_slice` and `sd_multiplier` controls for conditional and simple
  effects, including automatic SD slicing for moderators with many observed
  values.
* Added `rhat_summary()` for MCMC fits, returning a numeric R-hat vector with a
  compact printed summary.
* Expanded data-reshaping helpers: `to_long()` now supports multiple value
  columns, list-based column groups, and preserves input row order by default
  while still allowing sorted output with `sort = TRUE`.
* Added AD-compatible helper constructors `rtmb_vector()` and `rtmb_array()` for
  model code that needs mutable RTMB-compatible containers.
* Improved RTMB model setup error messages for common AD and NA/NaN failures.
* Improved wrapper behavior, including `rtmb_glmer(cwc = list(ID, "all"))`,
  hierarchical `lambda` in `rtmb_mdu()`, stronger prior validation, and more
  robust handling of non-finite VB optimization attempts.
* Improved MDU plotting and initialization, including principal-axis reference
  rotation and clearer radius display controls.
* Fixed several model-specific issues, including AD-compatible negative-binomial
  log densities and unequal-variance JZS t-test examples.

# BayesRTMB 0.1.1

* Addressed CRAN resubmission feedback by quoting 'RTMB' in the package title,
  adding method references to DESCRIPTION, and replacing `\dontrun{}` examples
  with `\donttest{}` where appropriate.
* Trimmed long-running examples for CRAN checks while retaining representative
  MCMC examples for correlation, t-test, and mixed-model workflows.
* Updated IRT post-estimation examples to use ordered response data explicitly.
* Fixed an AD-compatible negative-binomial log-density issue.
* Corrected the unequal-variance JZS t-test example and documentation.
