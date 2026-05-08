# Classic fit object

An R6 class representing the results of a classical (frequentist)
estimation.

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

- `bootstrap_results`:

  A matrix containing bootstrap samples, if applicable.

- `test_results`:

  List of additional test results (e.g., chisq.test).

- `view`:

  Character vector of parameter names to prioritize in summary.

- `par_unc`:

  Numeric vector of unconstrained parameter estimates.

- `vcov_unc`:

  Variance-covariance matrix of parameters in unconstrained space.

## Methods

### Public methods

- [`Classic_Fit$new()`](#method-Classic_Fit-new)

- [`Classic_Fit$compute_robust()`](#method-Classic_Fit-compute_robust)

- [`Classic_Fit$compute_bootstrap()`](#method-Classic_Fit-compute_bootstrap)

- [`Classic_Fit$.update_fit_with_vcov()`](#method-Classic_Fit-.update_fit_with_vcov)

- [`Classic_Fit$AIC()`](#method-Classic_Fit-AIC)

- [`Classic_Fit$BIC()`](#method-Classic_Fit-BIC)

- [`Classic_Fit$print()`](#method-Classic_Fit-print)

- [`Classic_Fit$logLik()`](#method-Classic_Fit-logLik)

- [`Classic_Fit$summary()`](#method-Classic_Fit-summary)

- [`Classic_Fit$savage_dickey()`](#method-Classic_Fit-savage_dickey)

- [`Classic_Fit$anova()`](#method-Classic_Fit-anova)

- [`Classic_Fit$lsmeans()`](#method-Classic_Fit-lsmeans)

- [`Classic_Fit$.calc_contrast()`](#method-Classic_Fit-.calc_contrast)

- [`Classic_Fit$.get_lsmeans_df()`](#method-Classic_Fit-.get_lsmeans_df)

- [`Classic_Fit$.construct_par_list()`](#method-Classic_Fit-.construct_par_list)

- [`Classic_Fit$clone()`](#method-Classic_Fit-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \`Classic_Fit\` object.

#### Usage

    Classic_Fit$new(
      model,
      fit,
      vcov = NULL,
      par_unc = NULL,
      vcov_unc = NULL,
      se_method = "wald",
      cluster = NULL,
      test_results = list(),
      view = NULL
    )

#### Arguments

- `model`:

  The \`RTMB_Model\` object.

- `fit`:

  The result of the estimation.

- `vcov`:

  Variance-covariance matrix of fixed effects.

- `par_unc`:

  Numeric vector of unconstrained parameter estimates.

- `vcov_unc`:

  Variance-covariance matrix of parameters in unconstrained space.

- `se_method`:

  Character; "wald", "robust", or "bootstrap".

- `cluster`:

  Character; cluster variable name.

- `test_results`:

  List of additional test results (e.g., chisq.test).

- `view`:

  Character vector of parameter names to prioritize in summary.

------------------------------------------------------------------------

### Method `compute_robust()`

Compute robust standard errors (sandwich estimator).

#### Usage

    Classic_Fit$compute_robust(cluster = NULL)

#### Arguments

- `cluster`:

  Character; variable name for clustering.

#### Returns

Self.

------------------------------------------------------------------------

### Method `compute_bootstrap()`

Compute bootstrap standard errors.

#### Usage

    Classic_Fit$compute_bootstrap(n_boot = 1000, cluster = NULL)

#### Arguments

- `n_boot`:

  Integer; number of samples.

- `cluster`:

  Character; clustering variable.

#### Returns

Self.

------------------------------------------------------------------------

### Method `.update_fit_with_vcov()`

(Internal) Update fit data frame with current vcov.

#### Usage

    Classic_Fit$.update_fit_with_vcov()

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

### Method [`summary()`](https://rdrr.io/r/base/summary.html)

Display a summary of the estimation results.

#### Usage

    Classic_Fit$summary(digits = 5, max_rows = 10)

#### Arguments

- `digits`:

  Number of digits to print for estimates.

- `max_rows`:

  Maximum number of rows to display in the coefficient table.

------------------------------------------------------------------------

### Method `savage_dickey()`

Calculate Bayes Factors using the Savage-Dickey density ratio method.

#### Usage

    Classic_Fit$savage_dickey(pars = NULL, null = 0, digits = 3)

#### Arguments

- `pars`:

  Character vector of parameter names to test.

- `null`:

  Numeric; the null value to test against (in constrained space).
  Default is 0.

- `digits`:

  Number of decimal places to round results.

#### Returns

A data frame containing the Bayes Factors and evidence descriptors.

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
      specs,
      pairwise = FALSE,
      simple = NULL,
      adjust = "bonferroni"
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
