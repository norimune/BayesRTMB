# Guidelines for Writing RTMB-Compatible Code

Models defined in `rtmb_code` rely on Automatic Differentiation (AD) via
the `RTMB` package. To ensure the model is differentiable and
numerically stable, specific coding practices must be followed.

## Details

**1. Automatic Differentiation and `advector`:** Parameters and
intermediate calculations involving them are treated as `advector`
objects. RTMB "records" all mathematical operations to compute
derivatives automatically. If a function strips these attributes or is
not differentiable, the gradient calculation will fail.

**2. Avoid Discrete Branching:** Standard R conditional statements like
`if (x > 0)` or [`ifelse()`](https://rdrr.io/r/base/ifelse.html) based
on parameter values do not provide derivatives.

- **Problem:** They create "jumps" in the likelihood surface.

- **Solution:** Use smooth approximations. For example, use
  [`fabs`](https://norimune.github.io/BayesRTMB/reference/fabs.md)
  instead of `abs`, or
  [`log_mix`](https://norimune.github.io/BayesRTMB/reference/log_mix.md)
  for mixture logic.

**3. Numerical Stability:** Computations in log-space are preferred to
prevent overflow or underflow. Use the following stable utilities
instead of raw algebraic expressions:

- [`log_sum_exp`](https://norimune.github.io/BayesRTMB/reference/log_sum_exp.md):
  For summing probabilities in log-space.

- [`log1p_exp`](https://norimune.github.io/BayesRTMB/reference/log1p_exp.md):
  For `log(1 + exp(x))`.

- [`inv_logit`](https://norimune.github.io/BayesRTMB/reference/inv_logit.md):
  For mapping real numbers to probabilities.

**4. Vectorization vs. Loops:**

- **Vectorization:** Highly recommended for performance. Standard R
  vectorized arithmetic (`+`, `-`, `*`, `/`, `log`, `exp`) works
  seamlessly with AD.

- **Loops:** Standard `for` loops are safe as long as the operations
  inside are differentiable.

- **Avoid `apply`:** Functions like `apply`, `sapply`, or `lapply` may
  sometimes strip AD attributes and should be replaced with vectorized
  operations or explicit loops.

**5. Matrix Operations:** Use specialized functions for matrix algebra
to maintain efficiency:

- [`quad_form_chol`](https://norimune.github.io/BayesRTMB/reference/quad_form_chol.md):
  For efficient quadratic forms.

- [`log_det_chol`](https://norimune.github.io/BayesRTMB/reference/log_det_chol.md):
  For log-determinants via Cholesky factors.

## See also

[`math_functions`](https://norimune.github.io/BayesRTMB/reference/math_functions.md),
[`distributions`](https://norimune.github.io/BayesRTMB/reference/distributions.md),
[`rtmb_code`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
