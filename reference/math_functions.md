# Mathematical and Matrix Utility Functions for RTMB Models

This page summarizes the utility functions provided to assist in writing
efficient and numerically stable model code within \`rtmb_code\`. These
functions are specifically designed to be compatible with RTMB's
Automatic Differentiation (AD).

## Details

**1. Link and Inverse Functions:**

- `logit(x)`: Computes the logit transformation `log(x/(1-x))`.

- `inv_logit(x)`: Computes the inverse logit (logistic) transformation
  `1/(1+exp(-x))`.

**2. Numerical Stability Utilities:** These functions use specialized
algorithms to prevent overflow/underflow in AD calculations.

- `log_sum_exp(x)`: Safely computes `log(sum(exp(x)))` using the
  log-sum-exp trick.

- `log1p_exp(x)`: Computes `log(1 + exp(x))` stably.

- `log1m_exp(x)`: Computes `log(1 - exp(x))` stably for `x < 0`.

- `log_softmax(x)`: Computes the log of the softmax function.

- `fabs(x)`: A smooth version of the absolute value function `abs(x)` to
  ensure differentiability at zero.

**3. Matrix and Vector Transformations:** Used to convert unconstrained
vectors into structured matrices (e.g., for factor analysis or
identification constraints).

- `sum_to_zero(x)`: Transforms a vector of length K-1 into a vector of
  length K that sums to zero.

- `to_lower_tri(x, M, D)`: Fills a matrix of size M x D with elements
  from vector `x` in a lower-triangular fashion.

- `to_centered_matrix(x, R, C)`: Creates an R x C matrix where each
  column sums to zero.

- `to_centered_tri(x, R, C)`: Creates an R x C matrix with column-wise
  sum-to-zero constraints on the lower elements (useful for
  identification in factor analysis).

**4. Linear Algebra for AD:**

- `log_det_chol(L)`: Calculates the log-determinant of a covariance
  matrix from its Cholesky factor `L`.

- `quad_form_chol(x, L)`: Computes the quadratic form x^T Sigma^(-1) x
  using the Cholesky factor `L`.

- `distance(x, y)`: Computes the Euclidean distance between two vectors
  with a small epsilon for stability.
