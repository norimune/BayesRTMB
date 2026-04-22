# Sort and display factor loadings neatly

Sort and display factor loadings neatly

## Usage

``` r
sort_loadings(loadings, cutoff = 0, round_digits = 3)
```

## Arguments

- loadings:

  Matrix, data frame, or list of factor loadings (if list, the first
  element is used)

- cutoff:

  Absolute loadings below this value will be displayed as blank (default
  is 0.0)

- round_digits:

  Number of decimal places to display (default is 3)

## Value

Sorted loading matrix (returned invisibly, allowing assignment to a
variable)
