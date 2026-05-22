# Make Best-Worst Pair Indices from Best and Worst Responses

Converts separate \`Best\` and \`Worst\` response matrices into
\`Y_dif\` pair indices for Best-Worst MDU models. The returned index is
the position of the ordered \`(best, worst)\` pair among all \`C \* (C -
1)\` position pairs generated from each row of \`sets\`.

## Usage

``` r
make_ydif_from_bw(Best, Worst, sets)
```

## Arguments

- Best:

  Matrix or data frame of best responses (N persons x P tasks). Values
  must be positions within the corresponding row of \`sets\`, from \`1\`
  to \`ncol(sets)\`.

- Worst:

  Matrix or data frame of worst responses (N persons x P tasks). Values
  must be positions within the corresponding row of \`sets\`, from \`1\`
  to \`ncol(sets)\`.

- sets:

  Matrix or data frame of presented item sets (P tasks x C items).

## Value

An integer matrix of pair indices (\`Y_dif\`).
