# Create an AD-compatible vector

Creates a fixed-length vector that is safe to use inside \`rtmb_code()\`
when subsequent assignments may involve RTMB automatic-differentiation
values. This is a safer alternative to \`numeric(length)\` in model
code.

## Usage

``` r
rtmb_vector(value = 0, length, seed = NULL)
```

## Arguments

- value:

  Initial value. Usually a scalar such as \`0\`.

- length:

  Length of the vector.

- seed:

  Optional AD value used to create the container. Usually this can be
  left as \`NULL\`; inside \`rtmb_code()\`, BayesRTMB automatically
  provides an AD seed from the model parameters when available.

## Value

An AD-compatible vector.
