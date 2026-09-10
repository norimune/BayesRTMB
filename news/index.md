# Changelog

## BayesRTMB 0.4.0

- Added `MCMC_Fit$continue_sampling()` to extend existing NUTS chains
  without another warmup. It resumes from each chain’s final
  unconstrained state, reuses the adapted step size and mass matrix, and
  can update the fit in place or return an independently extended copy.
- Corrected fixed-effect standard errors for Gaussian
  [`rtmb_lmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lmer.md)
  models fitted with `classic()`. REML fits now use the generalized
  least-squares covariance at the estimated variance components,
  matching `lmerTest`/`lme4`.
- Corrected `rtmb_corr(prior_flat())` marginal-likelihood normalization
  so Bayes factors for correlations use the normalized LKJ(1)/uniform
  correlation prior without changing the printed flat-prior model code.
- Made
  [`posterior_predict()`](https://norimune.github.io/BayesRTMB/reference/posterior_predict.md)
  and
  [`pp_check()`](https://norimune.github.io/BayesRTMB/reference/pp_check.md)
  report a clear error for `classic()` results, which do not store
  posterior draws.
- Fixed simulation-based standard errors for one-parameter models,
  including `optimize(se_method = "sampling")` on scalar binomial
  models.
- Changed `plot_mdu(distance = "auto")` to use Euclidean distance when
  no distance metadata is available. Stored model distance settings
  continue to take precedence.
- Clarified posterior-prediction examples in the Japanese wrapper
  vignette.

## BayesRTMB 0.3.0

CRAN release: 2026-08-20

- Added posterior predictive simulation and checks through
  [`posterior_predict()`](https://norimune.github.io/BayesRTMB/reference/posterior_predict.md)
  and
  [`pp_check()`](https://norimune.github.io/BayesRTMB/reference/pp_check.md).
  Regression wrappers automatically support continuous density checks,
  discrete binned-bar checks, scalar test statistics, predictor-based
  and fitted-value calibration scatter checks, and conditional,
  population-level, or newly simulated random effects. Custom models can
  provide replicated outcomes in a `generate` block.
- Added random-intercept mediation models to
  [`rtmb_mediation()`](https://norimune.github.io/BayesRTMB/reference/rtmb_mediation.md).
  A random intercept can be included in one equation or in multiple
  equations; random intercepts sharing a grouping variable are modeled
  jointly with an estimated correlation matrix.
- Added grand-mean centering (`gmc`/`centering`) and centering within
  cluster (`cwc`) to
  [`rtmb_mediation()`](https://norimune.github.io/BayesRTMB/reference/rtmb_mediation.md).
  Centering is applied to predictor uses while preserving response
  variables on their original scale; cluster means remain user-specified
  model terms.
- Corrected degrees of freedom for classical mediation models. Fixed
  Gaussian equations now use the rank of their own design matrix,
  derived effects inherit degrees of freedom from their contributing
  coefficients, and Gaussian random-intercept mediation models use
  Satterthwaite degrees of freedom by default.
- Added the setup-only `.data` binding for accessing the original object
  passed to
  [`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/rtmb_model.md).
  Wrapper-generated setup code now reads matrices and data frames
  directly, or reads related inputs such as responses, IDs, covariates,
  and choice sets from one named data list. Structural options remain
  visible as fixed assignments in the generated code.
- Added `std = TRUE` to the regression wrappers. It reports post-hoc
  standardized fixed-effect coefficients as `b_std`, using every column
  of the fixed-effect design matrix, including factor and interaction
  columns.
- Added
  [`summary_mcmc()`](https://norimune.github.io/BayesRTMB/reference/summary_mcmc.md)
  for summarizing numeric MCMC arrays arranged as iterations by chains
  by variables. It provides the same summary columns and print format as
  `MCMC_Fit$summary()`, with parameter and chain selection.
- Added exported
  [`center_grand_mean()`](https://norimune.github.io/BayesRTMB/reference/center_grand_mean.md)
  and
  [`center_within_cluster()`](https://norimune.github.io/BayesRTMB/reference/center_within_cluster.md)
  helpers for generated wrapper code and hand-written models.
- Made regression-wrapper `print_code()` expose data-frame columns
  instead of a nested data/formula payload, with explicit formula
  preprocessing and readable grand-mean and within-cluster centering
  helpers.
- Fixed `rtmb_table(x, y, data = ...)` so unquoted column names are
  resolved from the supplied data frame while preserving reproducible
  generated code.
- Simplified factor-analysis and multidimensional-unfolding
  `print_code()` output by assigning `nfactors` and `ndim` values
  directly to the internal constants `K` and `D`.
- Updated
  [`rtmb_corr()`](https://norimune.github.io/BayesRTMB/reference/rtmb_corr.md)
  normal priors to use the model-specific `mean_sd` and `sd_rate`
  aliases and to add an LKJ prior with `lkj_eta = 1` by default.
- Corrected Jacobian adjustments for constrained parameters supplied
  through `fixed`, so elements removed from the free-parameter map no
  longer contribute to optimization or sampling target densities.
- Updated the English and Japanese vignettes to document `.data`,
  unified named data inputs, structural constants in `setup`, and
  reproducible wrapper-generated code.

## BayesRTMB 0.2.4

CRAN release: 2026-07-24

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
