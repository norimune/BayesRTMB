# Common Features and Arguments of RTMB Wrapper Functions

The RTMB wrapper functions (\`rtmb_lm\`, \`rtmb_glm\`, \`rtmb_glmer\`,
\`rtmb_fa\`, etc.) share a unified interface designed to make Bayesian
and Frequentist inference accessible through familiar R formulas and
standard model specifications.

## Details

All wrapper functions in this package are built upon the same core
engine. This ensures that regardless of the model type, the workflow for
estimation, summary, and expansion remains consistent.

**1. Unified Inference Methods:** Every model object returned by a
wrapper function provides the following methods:

- `$classic()`: Performs standard frequentist estimation (MLE/REML)
  using \`prior_flat()\`. Supports robust standard errors and
  Satterthwaite/Between-Within degrees of freedom approximations.

- `$optimize()`: Performs Maximum A Posteriori (MAP) estimation
  (Empirical Bayes).

- `$sample()`: Performs MCMC sampling (NUTS via tmbstan) for full
  Bayesian inference.

- `$variational()`: Performs Variational Inference (ADVI) for fast
  posterior approximation.

**2. Prior API and Regularization:** You can specify the prior
distribution using the \`prior\` argument:

- [`prior_flat()`](https://norimune.github.io/BayesRTMB/reference/prior_flat.md):
  No additional prior (default). Required for `classic()` inference.

- [`prior_normal()`](https://norimune.github.io/BayesRTMB/reference/prior_normal.md):
  Manual normal/exponential priors for parameters.

- [`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md):
  Weakly informative priors based on `y_range`.

- [`prior_rhs()`](https://norimune.github.io/BayesRTMB/reference/prior_rhs.md):
  Regularized Horseshoe prior for continuous shrinkage.

- [`prior_ssp()`](https://norimune.github.io/BayesRTMB/reference/prior_ssp.md):
  Spike-and-Slab prior for sparse variable selection.

*Note: When using \`prior_weak()\`, \`prior_rhs()\`, or \`prior_ssp()\`,
you must specify `y_range = c(min, max)` to let the model set
appropriate global scales for the priors.*

**3. Model Comparison and Bayes Factors:** Model comparison is performed
exclusively via Bridge Sampling. To compare a full model against a
nested null model, use the `$bayes_factor(fixed = list(...))` method on
an MCMC fit object. For example, `fit$bayes_factor(fixed = list(b = 0))`
tests the hypothesis that the fixed effect `b` is exactly zero.

**4. Standard Errors and Degrees of Freedom:** The `$classic()` and
`$optimize()` methods provide advanced summary options:

- `se_method`: Calculate \`"robust"\` (Huber-White) or
  \`"cluster-robust"\` standard errors by providing a cluster variable.

- `df_method`: Use \`"satterthwaite"\` for mixed models or \`"bw"\`
  (Between-Within) for multilevel correlation models to obtain accurate
  p-values and confidence intervals.

**5. Fixed parameters:** Use `fixed = list(...)` during model
construction to fix model parameters to specified values. For example,
`fixed = list(delta = 0)` fixes the parameter `delta` to zero. This same
syntax is used in `$bayes_factor()` to define the restricted model.
