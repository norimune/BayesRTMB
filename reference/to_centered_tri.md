# Vector to centered triangular matrix (RTMB compatible)

Vector to centered triangular matrix (RTMB compatible)

## Usage

``` r
to_centered_tri(x, R, C)
```

## Arguments

- x:

  A numeric vector of appropriate length.

- R:

  Number of rows.

- C:

  Number of columns.

## Value

An R x C matrix with column-wise sum-to-zero constraints on lower
elements.
