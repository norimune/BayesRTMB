# Automatic Differentiation Variational Inference (ADVI)

Automatic Differentiation Variational Inference (ADVI)

## Usage

``` r
ADVI_method(
  model,
  par_list,
  pl_full,
  iter = 3000,
  tol_rel_obj = 0.001,
  window_size = 100,
  num_samples = 1000,
  alpha = 0.01,
  laplace = FALSE,
  print_freq = 500,
  method = c("meanfield", "fullrank", "hybrid"),
  update_progress = NULL,
  update_interval = 100
)
```

## Arguments

- model:

  An RTMB objective function object (\`ad_obj\`).

- par_list:

  A list defining the structure of parameters to be estimated.

- pl_full:

  A list defining the full structure of parameters including random
  effects.

- iter:

  Integer; fixed number of iterations for the optimization. Default is
  3000.

- tol_rel_obj:

  Numeric; relative tolerance for the ELBO to check convergence. Default
  is 0.001.

- window_size:

  Integer; size of the moving window to calculate the median ELBO.
  Default is 100.

- num_samples:

  Integer; number of posterior draws to generate after optimization.
  Default is 1000.

- alpha:

  Numeric; learning rate (step size) for the Adam optimizer. Default is
  0.01.

- laplace:

  Logical; whether Laplace approximation is used. Default is FALSE.

- print_freq:

  Integer; frequency of printing progress to the console. Set to 0 to
  disable. Default is 500.

- method:

  Character; type of variational distribution. One of "meanfield",
  "fullrank", or "hybrid". Default is "meanfield".

- update_progress:

  Optional function to update a progress bar.

- update_interval:

  Integer; interval for updating the progress bar. Default is 100.

## Value

A list containing \`fit\`, \`random_fit\`, \`elbo_history\`,
\`elbo_final\`, \`rel_obj_final\`, and \`converged\`.
