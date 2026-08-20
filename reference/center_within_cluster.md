# Center a Numeric Variable Within Clusters

Subtracts each cluster mean from a numeric vector. Missing cluster
values produce missing centered values. This helper is used by the
regression wrappers in their generated \`setup\` code and can also be
used in hand-written \`rtmb_code()\` models.

## Usage

``` r
center_within_cluster(x, cluster, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector to center.

- cluster:

  Vector identifying the cluster for each element of \`x\`.

- na.rm:

  Logical; whether missing values are removed when calculating cluster
  means.

## Value

A numeric vector with the same length as \`x\`.

## Examples

``` r
center_within_cluster(c(1, 3, 2, 6), c("a", "a", "b", "b"))
#> [1] -1  1 -2  2
```
