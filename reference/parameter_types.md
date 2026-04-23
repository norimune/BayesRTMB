# Parameter Types and Constraints in RTMB Models

When declaring parameters in the `parameters` block using
[`Dim()`](https://norimune.github.io/BayesRTMB/reference/Dim.md), you
can specify a `type` to automatically apply structural constraints.
These constraints ensure that parameters remain within their valid
mathematical space during estimation.

## Details

**Standard Types:**

- `"vector"`, `"matrix"`, `"array"`: Standard unconstrained containers.
  Use `lower` and `upper` for element-wise bounds.

**Ordering Constraints:**

- `"ordered"`: A vector where \$x_1 \< x_2 \< ... \< x_K\$.

- `"positive_ordered"`: A vector where \$0 \< x_1 \< x_2 \< ... \<
  x_K\$.

**Constrained Vectors:**

- `"simplex"`: A vector where all elements are \$ 0\$ and \$ x = 1\$.

- `"sum_to_zero"`: A vector of length \$K\$ where \$ x = 0\$. (Estimated
  using \$K-1\$ degrees of freedom).

**Correlation and Covariance Matrices:**

- `"corr_matrix"`: A symmetric, positive-definite correlation matrix
  (diagonal elements are 1).

- `"cov_matrix"`: A symmetric, positive-definite covariance matrix.

- `"CF_corr"`: The Cholesky factor of a correlation matrix.

- `"CF_cov"`: The Cholesky factor of a covariance matrix (diagonal
  elements are positive).

**Special Structural Matrices:**

- `"lower_tri"`: A lower-triangular matrix.

- `"positive_lower_tri"`: A lower-triangular matrix with positive
  diagonal elements.

- `"centered_matrix"`: A matrix where each column sums to zero.

- `"centered_tri"`: A triangular matrix with column-wise sum-to-zero
  constraints (often used for identification in factor analysis).

## See also

[`rtmb_code`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md),
[`math_functions`](https://norimune.github.io/BayesRTMB/reference/math_functions.md),
[`distributions`](https://norimune.github.io/BayesRTMB/reference/distributions.md)
