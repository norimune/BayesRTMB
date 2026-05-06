# Convert Long Data to Wide Format

A user-friendly wrapper around
[`stats::reshape`](https://rdrr.io/r/stats/reshape.html) to convert data
from long format back to wide format.

## Usage

``` r
to_wide(data, within = "Condition", value = "Value", id = NULL)
```

## Arguments

- data:

  A data frame in long format.

- within:

  The name of the column containing within-subjects factor levels (e.g.,
  "Condition").

- value:

  The name of the column containing measurement values (e.g., "Value").

- id:

  The identifier columns (e.g., "Subject", "Group"). If NULL, inferred
  as all columns except `within` and `value`.

## Value

A data frame in wide format.
