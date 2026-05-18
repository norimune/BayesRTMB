# RTMB-based Multidimensional Unfolding Wrapper

Fits a multidimensional unfolding model for preference/rating data. Rows
are persons or observations and columns are stimuli/items. The model
represents both row scores (\`theta\`) and item locations (\`delta\`) in
a shared D-dimensional space.

## Usage

``` r
rtmb_mdu(
  data,
  ndim = 2,
  distance = c("squared", "euclidean"),
  alpha = c("random", "fix"),
  method = c("rating", "Best", "Best-Worst", "MDS"),
  sets = NULL,
  prior = prior_flat(),
  y_range = NULL,
  init = NULL,
  fixed = NULL,
  view = NULL,
  distance_eps = 1e-04
)
```

## Arguments

- data:

  Numeric matrix or data frame (N rows x M items) for \`method =
  "rating"\`. For choice methods, a list containing \`Best\` and, for
  \`method = "Best-Worst"\`, \`Worst\`. For \`method = "Best-Worst"\`, a
  single matrix/data frame is treated as pre-coded \`Y_dif\` pair
  indices.

- ndim:

  Number of unfolding dimensions.

- distance:

  Character; \`"squared"\` uses squared Euclidean distance (default,
  often easier for optimization), while \`"euclidean"\` uses Euclidean
  distance.

- alpha:

  Character; \`"random"\` estimates item-specific alpha values as random
  effects (default), while \`"fix"\` estimates a single common alpha.

- method:

  Character; \`"rating"\` for continuous ratings, \`"Best"\` for
  best-only choice tasks, \`"Best-Worst"\` for best-worst choice tasks,
  or \`"MDS"\` for fitting a multidimensional scaling model to a
  distance matrix.

- sets:

  Matrix or data frame of presented item sets (P tasks x C items) for
  choice methods.

- prior:

  Prior configuration: \`prior_flat()\`, \`prior_normal()\`, or
  \`prior_weak()\`. \`prior_flat()\` creates a maximum-likelihood model
  suitable for \`classic()\`. The latent coordinates \`delta\` and
  \`theta\` are always treated as random effects with normal scale
  priors, similarly to IRT ability parameters.

- y_range:

  Optional response range for \`method = "rating"\`. If supplied with
  the default flat prior, \`prior_weak()\` is used. For choice methods
  (\`"Best"\` and \`"Best-Worst"\`), \`prior_weak()\` behaves like
  \`prior_normal()\` because there is no observed rating scale.

- init:

  Optional named list of initial values.

- fixed:

  Optional named list of parameter values to fix.

- view:

  Character vector of parameter names to prioritize in summaries.

- distance_eps:

  Small positive constant added to the distance.

## Value

An \`RTMB_Model\` object.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Simulate rating data for Multidimensional Unfolding (MDU)
  set.seed(123)
  N <- 50  # Number of persons
  M <- 10  # Number of items
  D <- 2   # Number of dimensions
  
  # True person and item coordinates in a 2D space
  theta <- matrix(rnorm(N * D), N, D)
  delta <- matrix(rnorm(M * D), M, D)
  
  # Generate distance-like ratings (smaller means more preferred)
  Y <- matrix(NA, N, M)
  for(i in 1:N) {
    for(j in 1:M) {
      Y[i, j] <- sum((theta[i,] - delta[j,])^2) + rnorm(1, 0, 0.5)
    }
  }

  # Fit a 2-dimensional MDU model
  fit_mdu <- rtmb_mdu(Y, ndim = 2, distance = "squared", method = "rating")
  
  # MAP estimation
  map_mdu <- fit_mdu$optimize()
  map_mdu$summary()
  
  # Note: MDU models have many parameters, so MCMC sampling might take time.
  # mcmc_mdu <- fit_mdu$sample(sampling = 500, warmup = 500, chains = 2)
} # }
```
