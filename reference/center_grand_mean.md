# Center a Numeric Variable Around Its Grand Mean

Subtracts the overall mean from a numeric vector. This helper is used by
the regression wrappers in their generated \`setup\` code and can also
be used in hand-written \`rtmb_code()\` models.

## Usage

``` r
center_grand_mean(x, na.rm = TRUE)
```

## Arguments

- x:

  Numeric vector to center.

- na.rm:

  Logical; whether missing values are removed when calculating the mean.

## Value

A numeric vector with the same length as \`x\`.

## Examples

``` r
center_grand_mean(c(1, 2, 3))
#> [1] -1  0  1
```
