# Convert Wide Data to Long Format

A highly intuitive wrapper around
[`stats::reshape`](https://rdrr.io/r/stats/reshape.html) designed for
psychological research. It converts data from wide format to long format
by identifying within-subjects factors.

## Usage

``` r
to_long(data, within = NULL, label = NULL, value = "Value", id = NULL)
```

## Arguments

- data:

  A data frame in wide format.

- within:

  The columns to gather into long format. Can be:

  - A character vector of column names.

  - A range string like `"time1:time4"`.

  - A prefix string like `"time"` (matches all columns starting with
    "time").

- label:

  The name for the new column that will contain the level names (e.g.,
  "Time"). If NULL, it defaults to the prefix used in `within` or
  "Condition".

- value:

  The name for the new column that will contain the measurement values.
  Default is "Value".

- id:

  The identifier columns that should be repeated for each row. If NULL
  (default), all columns NOT specified in `within` are treated as IDs.

## Value

A data frame in long format.
