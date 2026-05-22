# Classic fit object

An R6 class representing the results of a classical (frequentist)
estimation.

## Super class

[`BayesRTMB::RTMB_Fit_Base`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.md)
-\> `Classic_Fit`

## Public fields

- `model`:

  The \`RTMB_Model\` object used for estimation.

- `fit`:

  The result of the estimation (dataframe or lm object).

- `par`:

  a named list of parameter estimates.

- `vcov`:

  Variance-covariance matrix of fixed effects.

- `se_method`:

  Character string specifying the method used for standard errors.

- `cluster`:

  Character string specifying the cluster variable name, if any.

- `robust_type`:

  Character string specifying the robust standard error correction type,
  if used.

- `bootstrap_results`:

  A matrix containing bootstrap samples, if applicable.

- `test_results`:

  List of additional test results (e.g., chisq.test).

- `view`:

  Character vector of parameter names to prioritize in summary.

- `par_vec`:

  Numeric vector of parameter estimates.

- `objective`:

  Final objective value.

- `log_lik`:

  Log-likelihood.

- `restricted_log_lik`:

  Restricted log-likelihood, if classical estimation used REML-style
  marginalization.

- `rss`:

  Residual sum of squares for classical linear models, if available.

- `df_residual`:

  Residual degrees of freedom for classical linear models, if available.

- `convergence`:

  Convergence code.

- `sd_rep`:

  TMB sdreport object.

- `df_fixed`:

  Dataframe of fixed effects results.

- `random_effects`:

  Dataframe of random effects results.

- `df_transform`:

  Dataframe of transformed parameters.

- `df_generate`:

  Dataframe of generated quantities.

- `opt_history`:

  Dataframe of optimization history.

- `transform`:

  List of transformed parameters.

- `generate`:

  List of generated quantities.

- `se_samples`:

  List of simulated samples for SE estimation.

- `par_unc`:

  Numeric vector of unconstrained parameter estimates.

- `vcov_unc`:

  Variance-covariance matrix of parameters in unconstrained space.

- `ci_method`:

  Method used for CI estimation.

- `laplace`:

  Whether Laplace approximation was used.

- `map`:

  Parameter mapping used.

- `df_method`:

  Character string specifying the degrees of freedom calculation method.

- `idx_fix_active`:

  Numeric vector; mapping between active parameters and full
  unconstrained vector.

- `show_df`:

  Logical; whether to display degrees of freedom in the summary output.

## Methods

### Public methods

- [`Classic_Fit$get_point_estimate()`](#method-Classic_Fit-get_point_estimate)

- [`Classic_Fit$new()`](#method-Classic_Fit-new)

- [`Classic_Fit$robust_se()`](#method-Classic_Fit-robust_se)

- [`Classic_Fit$compute_robust()`](#method-Classic_Fit-compute_robust)

- [`Classic_Fit$bootstrap()`](#method-Classic_Fit-bootstrap)

- [`Classic_Fit$compute_bootstrap()`](#method-Classic_Fit-compute_bootstrap)

- [`Classic_Fit$AIC()`](#method-Classic_Fit-AIC)

- [`Classic_Fit$BIC()`](#method-Classic_Fit-BIC)

- [`Classic_Fit$WAIC()`](#method-Classic_Fit-WAIC)

- [`Classic_Fit$diagnose()`](#method-Classic_Fit-diagnose)

- [`Classic_Fit$print()`](#method-Classic_Fit-print)

- [`Classic_Fit$logLik()`](#method-Classic_Fit-logLik)

- [`Classic_Fit$EAP()`](#method-Classic_Fit-EAP)

- [`Classic_Fit$MAP()`](#method-Classic_Fit-MAP)

- [`Classic_Fit$summary()`](#method-Classic_Fit-summary)

- [`Classic_Fit$anova()`](#method-Classic_Fit-anova)

- [`Classic_Fit$lsmeans()`](#method-Classic_Fit-lsmeans)

- [`Classic_Fit$.calc_contrast()`](#method-Classic_Fit-.calc_contrast)

- [`Classic_Fit$.get_lsmeans_df()`](#method-Classic_Fit-.get_lsmeans_df)

- [`Classic_Fit$.construct_par_list()`](#method-Classic_Fit-.construct_par_list)

- [`Classic_Fit$clone()`](#method-Classic_Fit-clone)

Inherited methods

- [`BayesRTMB::RTMB_Fit_Base$estimate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-estimate)
- [`BayesRTMB::RTMB_Fit_Base$fa_rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-fa_rotate)
- [`BayesRTMB::RTMB_Fit_Base$rotate()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.html#method-rotate)

------------------------------------------------------------------------

### Method `get_point_estimate()`

Get point estimate for a target parameter.

#### Usage

    Classic_Fit$get_point_estimate(target, ...)

#### Arguments

- `target`:

  Target parameter name.

- `...`:

  Additional arguments, ignored for classic fits.

#### Returns

Matrix, array, vector, or scalar point estimate.

------------------------------------------------------------------------

### Method `new()`

Create a new \`Classic_Fit\` object.

#### Usage

    Classic_Fit$new(
      model,
      par_vec = NULL,
      par = NULL,
      objective = NULL,
      log_lik = NULL,
      restricted_log_lik = NULL,
      convergence = NULL,
      sd_rep = NULL,
      df_fixed = NULL,
      random_effects = NULL,
      df_transform = NULL,
      df_generate = NULL,
      opt_history = NULL,
      transform = NULL,
      generate = NULL,
      se_samples = NULL,
      par_unc = NULL,
      vcov_unc = NULL,
      ci_method = "wald",
      laplace = TRUE,
      map = NULL,
      test_results = list(),
      view = NULL,
      fit = NULL,
      vcov = NULL,
      df_method = "auto",
      idx_fix_active = NULL,
      show_df = TRUE,
      rss = NULL,
      df_residual = NULL,
      ...
    )

#### Arguments

- `model`:

  The \`RTMB_Model\` object.

- `par_vec`:

  Numeric vector of parameter estimates.

- `par`:

  List of parameter estimates.

- `objective`:

  Final objective value.

- `log_lik`:

  Full log-likelihood used for information criteria.

- `restricted_log_lik`:

  Restricted log-likelihood, if classical estimation used REML-style
  marginalization.

- `convergence`:

  Convergence code.

- `sd_rep`:

  TMB sdreport object.

- `df_fixed`:

  Dataframe of fixed effects results.

- `random_effects`:

  Dataframe of random effects results.

- `df_transform`:

  Dataframe of transformed parameters.

- `df_generate`:

  Dataframe of generated quantities.

- `opt_history`:

  Dataframe of optimization history.

- `transform`:

  List of transformed parameters.

- `generate`:

  List of generated quantities.

- `se_samples`:

  List of simulated samples for SE estimation.

- `par_unc`:

  Parameter vector on unconstrained scale.

- `vcov_unc`:

  Covariance matrix on unconstrained scale.

- `ci_method`:

  Method used for CI estimation.

- `laplace`:

  Whether Laplace approximation was used.

- `map`:

  Parameter mapping used.

- `test_results`:

  List of additional test results (e.g., chisq.test).

- `view`:

  Character vector of parameter names to prioritize in summary.

- `fit`:

  Legacy argument for backward compatibility (maps to df_fixed).

- `vcov`:

  Variance-covariance matrix of parameters.

- `df_method`:

  Method for degrees of freedom calculation.

- `idx_fix_active`:

  Numeric vector; mapping between active parameters and full
  unconstrained vector.

- `show_df`:

  Logical; whether to display degrees of freedom in the summary output.

- `rss`:

  Residual sum of squares for classical linear models, if available.

- `df_residual`:

  Residual degrees of freedom for classical linear models, if available.

- `...`:

  Additional arguments passed to the constructor.

------------------------------------------------------------------------

### Method `robust_se()`

Compute robust standard errors (sandwich estimator).

#### Usage

    Classic_Fit$robust_se(
      cluster = NULL,
      type = c("HC3", "HC0", "HC1", "CR1", "CR0"),
      inplace = FALSE,
      ...
    )

#### Arguments

- `cluster`:

  Character; variable name for clustering.

- `type`:

  Character; "HC3" (default), "HC0", or "HC1" for non-clustered robust
  SE; "CR1" or "CR0" for clustered robust SE.

- `inplace`:

  Logical; if FALSE, return a robust-SE copy without modifying the
  original fit. If TRUE, update the object in place.

- `...`:

  Additional arguments.

#### Returns

Self.

------------------------------------------------------------------------

### Method `compute_robust()`

(Deprecated) Use robust_se() instead.

#### Usage

    Classic_Fit$compute_robust(...)

#### Arguments

- `...`:

  Arguments passed to robust_se(), including inplace.

#### Returns

Self.

------------------------------------------------------------------------

### Method `bootstrap()`

Compute nonparametric bootstrap standard errors and confidence
intervals. Currently implemented only for mediation models. Bootstrap
refits mediation equations with base R lm.fit()/glm.fit() when possible,
and falls back to RTMB refits for unsupported families.

#### Usage

    Classic_Fit$bootstrap(n_boot = 1000, seed = NULL, inplace = FALSE, ...)

#### Arguments

- `n_boot`:

  Integer; number of bootstrap samples.

- `seed`:

  Optional integer random seed.

- `inplace`:

  Logical; if FALSE, return a bootstrap-SE copy without modifying the
  original fit. If TRUE, update the object in place.

- `...`:

  Additional arguments.

#### Returns

A Classic_Fit object.

------------------------------------------------------------------------

### Method `compute_bootstrap()`

Compute nonparametric bootstrap standard errors and confidence
intervals. Currently implemented only for mediation models. Bootstrap
refits mediation equations with base R lm.fit()/glm.fit() when possible,
and falls back to RTMB refits for unsupported families.

#### Usage

    Classic_Fit$compute_bootstrap(n_boot = 1000, seed = NULL, inplace = FALSE, ...)

#### Arguments

- `n_boot`:

  Integer; number of bootstrap samples.

- `seed`:

  Optional integer random seed.

- `inplace`:

  Logical; if FALSE, return a bootstrap-SE copy without modifying the
  original fit. If TRUE, update the object in place.

- `...`:

  Additional arguments.

#### Returns

A Classic_Fit object.

------------------------------------------------------------------------

### Method [`AIC()`](https://rdrr.io/r/stats/AIC.html)

Get the AIC of the fitted model.

#### Usage

    Classic_Fit$AIC()

------------------------------------------------------------------------

### Method [`BIC()`](https://rdrr.io/r/stats/AIC.html)

Get the BIC of the fitted model.

#### Usage

    Classic_Fit$BIC()

------------------------------------------------------------------------

### Method `WAIC()`

WAIC is not defined for classical fits.

#### Usage

    Classic_Fit$WAIC()

------------------------------------------------------------------------

### Method `diagnose()`

Run basic diagnostics for the classical fit.

#### Usage

    Classic_Fit$diagnose(...)

#### Arguments

- `...`:

  Additional arguments passed to \`diagnose_classic_fit()\`.

#### Returns

A \`diagnose_BayesRTMB\` object.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the fit results.

#### Usage

    Classic_Fit$print(...)

#### Arguments

- `...`:

  Additional arguments passed to \`summary()\`.

------------------------------------------------------------------------

### Method [`logLik()`](https://rdrr.io/r/stats/logLik.html)

Get the Log-Likelihood of the fitted model.

#### Usage

    Classic_Fit$logLik()

------------------------------------------------------------------------

### Method `EAP()`

EAP estimates are not available for Classic_Fit.

#### Usage

    Classic_Fit$EAP(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### Method `MAP()`

MAP estimates are not available for Classic_Fit.

#### Usage

    Classic_Fit$MAP(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Display a summary of the estimation results.

#### Usage

    Classic_Fit$summary(view = NULL, digits = 5, max_rows = 10, ranef = FALSE)

#### Arguments

- `view`:

  Character vector of parameter names to prioritize or filter by.

- `digits`:

  Number of digits to print for estimates.

- `max_rows`:

  Maximum number of rows to display in the coefficient table.

- `ranef`:

  Logical; whether to show random effects in the summary.

------------------------------------------------------------------------

### Method [`anova()`](https://rdrr.io/r/stats/anova.html)

Perform ANOVA (Wald F-tests / Chisq-tests) on the fitted model.

#### Usage

    Classic_Fit$anova(method = c("reml", "ls"), type = 3)

#### Arguments

- `method`:

  Character; "reml" (standard) or "ls" (experimental).

- `type`:

  Integer; Type of Sum of Squares (only Type III supported currently).

#### Returns

A data frame containing the ANOVA table.

------------------------------------------------------------------------

### Method [`lsmeans()`](https://norimune.github.io/BayesRTMB/reference/lsmeans.md)

Calculate Least Squares Means (Marginal Means) and contrasts.

#### Usage

    Classic_Fit$lsmeans(
      specs = NULL,
      pairwise = FALSE,
      simple = NULL,
      adjust = "holm",
      protect = FALSE
    )

#### Arguments

- `specs`:

  Character vector of factors to calculate means for.

- `pairwise`:

  Logical; whether to perform pairwise comparisons.

- `simple`:

  Character vector of factors to hold constant for simple main effects.

- `adjust`:

  Character; p-value adjustment method (e.g., "bonferroni", "holm",
  "none").

- `protect`:

  Logical; whether to use hierarchical (protected) testing.

#### Returns

A data frame containing the marginal means or contrasts.

------------------------------------------------------------------------

### Method `.calc_contrast()`

(Internal) Calculate metrics for a contrast.

#### Usage

    Classic_Fit$.calc_contrast(L, specs)

#### Arguments

- `L`:

  Contrast matrix.

- `specs`:

  Variable names for lsmeans.

#### Returns

A data frame with estimate, SE, df, t-value, and p-value.

------------------------------------------------------------------------

### Method `.get_lsmeans_df()`

(Internal) Get representative DF for lsmeans.

#### Usage

    Classic_Fit$.get_lsmeans_df(specs)

#### Arguments

- `specs`:

  Variable names for lsmeans.

#### Returns

Degrees of freedom.

------------------------------------------------------------------------

### Method `.construct_par_list()`

(Internal) Construct a list of parameters from the fit.

#### Usage

    Classic_Fit$.construct_par_list(fit)

#### Arguments

- `fit`:

  The fit result (dataframe or lm object).

#### Returns

A named list of parameters.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    Classic_Fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
