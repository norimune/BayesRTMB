# Make Best and Worst Responses from Best-Worst Pair Indices

Converts a matrix of Best-Worst pair indices (\`Y_dif\`) into separate
\`Best\` and \`Worst\` data frames. The returned values are positions
within each row of \`sets\`, not the item labels stored in \`sets\`.

## Usage

``` r
make_bw_from_ydif(Y_dif, sets)

restore_bw_from_ydif(Y_dif, sets)
```

## Arguments

- Y_dif:

  Matrix or data frame of pair indices (N persons x P tasks). Each value
  must be an integer from 1 to \`C \* (C - 1)\`, where \`C\` is the
  number of items per task.

- sets:

  Matrix or data frame of presented item sets (P tasks x C items).

## Value

A list with two data frames: \`Best\` and \`Worst\`.
