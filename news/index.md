# Changelog

## BayesRTMB 0.2.4

- Allowed
  [`rtmb_lm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lm.md),
  [`rtmb_glm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glm.md),
  [`rtmb_lmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lmer.md),
  and
  [`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md)
  to resolve bare variable names from the formula environment when
  `data` is omitted. Formulas using `$`, `[[`, or `.` continue to
  require explicit data.
- Corrected JZS t-test priors to include the Jeffreys scale prior (p()
  /), including paired and unequal-variance models.
- Corrected exponential-prior scaling across wrappers. The default
  `sigma_rate` and `tau_rate` in
  [`prior_normal()`](https://norimune.github.io/BayesRTMB/reference/prior_normal.md)
  are now `1 / 5`, giving exponential priors with mean 5, and
  wrapper-specific aliases are applied consistently.
- Added the default discrimination prior `a ~ exponential(1 / 2)` when
  [`prior_normal()`](https://norimune.github.io/BayesRTMB/reference/prior_normal.md)
  is used with IRT models.
- Improved weak-prior calibration from `y_range` in mixture and
  latent-rank models, including response-specific ranges for
  multivariate outcomes.
- Improved multivariate normal log-density performance and removed an
  unintended blank line from optimization progress output.
- Updated English and Japanese vignettes for the new defaults and
  wrapper behavior.

## BayesRTMB 0.2.3

CRAN release: 2026-07-13

- Fixed factor-analysis model construction on platforms where RTMB is
  sensitive to AD matrix containers that include structural zeros from
  lower-triangular parameters. The
  [`rtmb_fa()`](https://norimune.github.io/BayesRTMB/reference/rtmb_fa.md)
  wrapper now avoids reading upper-triangular structural-zero entries of
  `lower_tri` loading matrices and constructs constrained AD matrices
  with
  [`rtmb_array()`](https://norimune.github.io/BayesRTMB/reference/rtmb_array.md).
- Simplified the runnable
  [`rtmb_fa()`](https://norimune.github.io/BayesRTMB/reference/rtmb_fa.md)
  example to a one-factor model. Advanced factor-analysis workflows
  remain covered by documentation and CI regression checks.

## BayesRTMB 0.2.2

- Added response-time distributions for model code:
  [`exp_mod_normal_lpdf()`](https://norimune.github.io/BayesRTMB/reference/exp_mod_normal_lpdf.md)
  and
  [`diffusion_lpdf()`](https://norimune.github.io/BayesRTMB/reference/diffusion_lpdf.md),
  with sampling syntax support via `exp_mod_normal(...)` and
  `obs(RT, Choice) ~ diffusion(...)`.
- Added `obs(...)` sampling syntax for multivariate observed values on
  the left side of `~`.
- Improved setup helper capture so functions referenced in `setup` are
  more reliably available when building models and running parallel
  workers.
- Improved factor-analysis AD robustness and rotation output naming.

## BayesRTMB 0.2.1

CRAN release: 2026-06-23

- Added
  [`upgrade_fit()`](https://norimune.github.io/BayesRTMB/reference/upgrade_fit.md)
  to rebuild saved MCMC, VB, MAP, and classic fit objects with the
  currently loaded class definitions, optionally upgrading their
  embedded model objects as well.
- Improved
  [`rtmb_vector()`](https://norimune.github.io/BayesRTMB/reference/rtmb_vector.md)
  and
  [`rtmb_array()`](https://norimune.github.io/BayesRTMB/reference/rtmb_array.md)
  tape construction time by automatically reusing an AD seed from model
  parameters when available.
- Made
  [`log_sum_exp()`](https://norimune.github.io/BayesRTMB/reference/log_sum_exp.md),
  [`softmax()`](https://norimune.github.io/BayesRTMB/reference/softmax.md),
  and
  [`log_softmax()`](https://norimune.github.io/BayesRTMB/reference/log_softmax.md)
  work more reliably with RTMB automatic-differentiation values,
  including baseline-category patterns such as `softmax(c(0, eta))`
  inside
  [`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md).
- Improved wrapper-generated model code to use AD-compatible
  [`rtmb_vector()`](https://norimune.github.io/BayesRTMB/reference/rtmb_vector.md)
  and
  [`rtmb_array()`](https://norimune.github.io/BayesRTMB/reference/rtmb_array.md)
  containers in loop-filled generated quantities and generated
  likelihood contributions where needed.
- Improved `report()` handling in transformed and generated quantities,
  including namespaced `BayesRTMB::report()` calls and wrapper-generated
  `print_code()` output.
- Changed VB point estimates to use only the best variational estimate
  by default, aligning `EAP()`, `MAP()`, and rotation references with
  the selected best ELBO run while still allowing explicit `chains` or
  `best_chains` selection.
- Made `EAP()` and `MAP()` drop their list wrapper by default when a
  single parameter is requested, matching the behavior of `estimate()`.
- Added optional taped evaluation for transformed parameters and
  generated quantities, with automatic fallback to R evaluation when
  taping is not possible.
- Improved parallel worker robustness by reducing exported globals and
  preserving wrapper setup environments needed by generated model code.
- Improved MCMC runtime behavior by caching metric calculations,
  speeding up retained draw conversion, and refining progress checks.
- Changed bootstrap progress reporting to use percentage-style progress
  output, consistent with other long-running workflows.
- Added diagnostic recommendations to help interpret common fitting
  warnings.
- Improved matrix-valued Gaussian process log-density evaluation.
- Updated MDU defaults and internals, including Euclidean distance as
  the default MDU distance and more explicit use of namespace-qualified
  factor rotations.

## BayesRTMB 0.2.0

- Reworked NUTS sampling internals with Stan-style multinomial tree
  expansion, warmup diagnostics, Stan-window metric adaptation, and
  support for diagonal, dense, hybrid, and automatic metric selection.
- Improved MCMC diagnostics by reporting divergence counts and
  percentages, metric auto-selection details, warmup summaries, metric
  condition numbers, and positive-definite fallback counts.
- Added configurable progress output for MCMC and VB workflows,
  including streamed message-style progress for parallel workers and
  percentage reporting.
- Added delta-method standard errors and confidence intervals for
  [`conditional_effects()`](https://norimune.github.io/BayesRTMB/reference/conditional_effects.md)
  and
  [`simple_effects()`](https://norimune.github.io/BayesRTMB/reference/simple_effects.md)
  with optimized and classic fits;
  [`simple_effects()`](https://norimune.github.io/BayesRTMB/reference/simple_effects.md)
  for classic fits now also reports `df`, `t value`, and `Pr`.
- Added `sd_slice` and `sd_multiplier` controls for conditional and
  simple effects, including automatic SD slicing for moderators with
  many observed values.
- Added
  [`rhat_summary()`](https://norimune.github.io/BayesRTMB/reference/rhat_summary.md)
  for MCMC fits, returning a numeric R-hat vector with a compact printed
  summary.
- Expanded data-reshaping helpers:
  [`to_long()`](https://norimune.github.io/BayesRTMB/reference/to_long.md)
  now supports multiple value columns, list-based column groups, and
  preserves input row order by default while still allowing sorted
  output with `sort = TRUE`.
- Added AD-compatible helper constructors
  [`rtmb_vector()`](https://norimune.github.io/BayesRTMB/reference/rtmb_vector.md)
  and
  [`rtmb_array()`](https://norimune.github.io/BayesRTMB/reference/rtmb_array.md)
  for model code that needs mutable RTMB-compatible containers.
- Improved RTMB model setup error messages for common AD and NA/NaN
  failures.
- Improved wrapper behavior, including
  `rtmb_glmer(cwc = list(ID, "all"))`, hierarchical `lambda` in
  [`rtmb_mdu()`](https://norimune.github.io/BayesRTMB/reference/rtmb_mdu.md),
  stronger prior validation, and more robust handling of non-finite VB
  optimization attempts.
- Improved MDU plotting and initialization, including principal-axis
  reference rotation and clearer radius display controls.
- Fixed several model-specific issues, including AD-compatible
  negative-binomial log densities and unequal-variance JZS t-test
  examples.

## BayesRTMB 0.1.1

CRAN release: 2026-06-01

- Addressed CRAN resubmission feedback by quoting ‘RTMB’ in the package
  title, adding method references to DESCRIPTION, and replacing
  `\dontrun{}` examples with `\donttest{}` where appropriate.
- Trimmed long-running examples for CRAN checks while retaining
  representative MCMC examples for correlation, t-test, and mixed-model
  workflows.
- Updated IRT post-estimation examples to use ordered response data
  explicitly.
- Fixed an AD-compatible negative-binomial log-density issue.
- Corrected the unequal-variance JZS t-test example and documentation.
