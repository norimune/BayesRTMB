# Probability Distributions for RTMB Models

This page summarizes the probability distributions available within the
\`rtmb_code\` block. These functions are designed to be used with the
Stan-like sampling syntax (\`~\`), which internally adds the log-density
to the model's total log-posterior.

## Details

**Syntax Styles:** In \`rtmb_code\`, you can specify distributions in
two ways:

- **Sampling Syntax (Recommended):** `y ~ normal(mu, sigma)`

- **Explicit Function Call:**
  `target <- target + normal_lpdf(y, mu, sigma)`

**Continuous Distributions (LPDF):**

- `normal(mean, sd)`: Normal distribution.

- `lognormal(meanlog, sdlog)`: Lognormal distribution.

- `exponential(rate)`: Exponential distribution.

- `cauchy(location, scale)`: Cauchy distribution.

- `student_t(df, mu, sigma)`: Student's t-distribution.

- `gamma(shape, rate)`: Gamma distribution.

- `inverse_gamma(shape, scale)`: Inverse-gamma distribution.

- `beta(a, b)`: Beta distribution.

**Discrete Distributions (LPMF):**

- `bernoulli(prob)` / `bernoulli_logit(eta)`: Binary outcomes.

- `binomial(size, prob)` / `binomial_logit(size, eta)`: Binomial
  outcomes.

- `poisson(mean)`: Poisson count data.

- `neg_binomial_2(mu, size)`: Negative binomial (mean/dispersion
  parameterization).

- `ordered_logistic(eta, cutpoints)`: Ordered categorical outcomes.

**Multivariate and Matrix Distributions:**

- `multi_normal(mean, Sigma)`: Standard multivariate normal
  distribution.

- `lkj_corr(eta)`: LKJ prior for correlation matrices.

- `dirichlet(alpha)`: Dirichlet distribution for simplexes.

- `lower_tri_normal(mean, sd)`: Normal distribution for elements of a
  lower-triangular matrix.

- `centered_tri_multi_normal(sigma)`: Multivariate normal for centered
  triangular matrices (used in identification constraints).

- `sufficient_multi_normal_fa(S_mat, N, y_bar, mu, psi, Lambda)`: Factor
  analysis likelihood using sufficient statistics (highly efficient for
  large sample sizes).

**Vectorization:** Most univariate distributions are vectorized. If `y`
and `mu` are vectors, `y ~ normal(mu, sigma)` will calculate the sum of
log-densities for all elements efficiently.
