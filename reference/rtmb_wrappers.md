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

- `$optimize()`: Performs MAP estimation (comparable to MLE).

- `$sample()`: Performs MCMC sampling (NUTS) for full Bayesian
  inference.

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

**3. Weakly Informative Priors (`y_range`):** By providing the
theoretical range of your response variable via `y_range` while keeping
the default `prior = prior_flat()`, the wrappers automatically switch to
[`prior_weak()`](https://norimune.github.io/BayesRTMB/reference/prior_weak.md).
These priors are designed to be broad enough to cover any reasonable
value but narrow enough to stabilize the estimation and prevent the
sampler from wandering into non-sensical parameter space.

**4. Fixed vs. Random Effects:** For mixed-effect models (e.g.,
`rtmb_glmer`), random effects are marginalized using the Laplace
approximation during `$optimize(laplace = TRUE)`. When using
`$sample()`, random effects are treated as unknown parameters and
sampled alongside fixed effects.

**5. Fixed parameters:** Use `fixed = list(...)` to fix model parameters
to specified values. For example, `fixed = list(delta = 0)` fixes the
parameter `delta` to zero, and `fixed = list("b[x]" = 0)` fixes one
element of a vector parameter.
