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
  method = c("rating", "Best", "Best-Worst"),
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
  best-only choice tasks, or \`"Best-Worst"\` for best-worst choice
  tasks.

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
