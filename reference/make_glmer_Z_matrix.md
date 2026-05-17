# Reconstruct an Observation-Level Random-Effect Design Matrix

Converts the transposed block random-effect design matrix used
internally by \`rtmb_glmer()\` into an observation-level matrix. This is
kept as a small utility for users who want to work from the block matrix
returned by \`make_glmer_re_terms()\`.

## Usage

``` r
make_glmer_Z_matrix(Zt, group_idx, N = length(group_idx))
```

## Arguments

- Zt:

  Transposed block random-effect design matrix.

- group_idx:

  Integer group index for each observation.

- N:

  Number of observations. Defaults to \`length(group_idx)\`.

## Value

A matrix with \`N\` rows and one column per random-effect coefficient.
