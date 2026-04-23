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

**2. Regularization and Variable Selection:** For rtmb_lm, rtmb_glm, and
rtmb_glmer, you can handle high-dimensional predictors (where the number
of variables is large relative to the sample size) using the penalty
argument:

- `"none"`: Standard flat or weakly informative priors.

- `"rhs"`: Regularized Horseshoe prior for continuous shrinkage.

- `"ssp"`: Spike-and-Slab prior for sparse variable selection.

*Note: When using regularization, you must specify
`y_range = c(min, max)` to let the model set appropriate global scales
for the priors.*

**3. Weakly Informative Priors (`y_range`):** By providing the
theoretical range of your response variable via `y_range`, the wrappers
automatically construct "Weakly Informative Priors". These priors are
designed to be broad enough to cover any reasonable value but narrow
enough to stabilize the estimation and prevent the sampler from
wandering into non-sensical parameter space.

**4. Fixed vs. Random Effects:** For mixed-effect models (e.g.,
`rtmb_glmer`), random effects are marginalized using the Laplace
approximation during `$optimize(laplace = TRUE)`. When using
`$sample()`, random effects are treated as unknown parameters and
sampled alongside fixed effects.

**5. Null Model Creation (`null`):** You can specify a `null` argument
(e.g., `null = "x1 ~ normal(0, 0.1)"`) in the wrappers to simultaneously
create a restricted version of your model. This is particularly useful
for computing Bayes Factors or performing model comparisons.
