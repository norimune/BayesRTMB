# Create Initial Values for Multidimensional Unfolding

Create Initial Values for Multidimensional Unfolding

## Usage

``` r
make_init_mdu(
  Y,
  D,
  distance = c("squared", "euclidean"),
  alpha = c("random", "fix"),
  distance_eps = 1e-04,
  min_sigma = 0.1
)
```

## Arguments

- Y:

  Numeric matrix or data frame (N rows x M items).

- D:

  Number of unfolding dimensions.

- distance:

  Character; \`"squared"\` or \`"euclidean"\`.

- alpha:

  Character; \`"random"\` for item-specific alpha initial values or
  \`"fix"\` for a common alpha initial value.

- distance_eps:

  Small positive constant added to the distance.

- min_sigma:

  Minimum initial residual standard deviation.

## Value

A named list of initial values.
