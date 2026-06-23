# Create an AD-compatible array

Creates a fixed-size vector, matrix, or higher-dimensional array that is
safe to use inside \`rtmb_code()\` when subsequent assignments may
involve RTMB automatic-differentiation values. This is a safer
alternative to \`array(value, dim)\` or \`matrix(value, ...)\` in model
code.

## Usage

``` r
rtmb_array(value = 0, dim, seed = NULL)
```

## Arguments

- value:

  Initial value. Usually a scalar such as \`0\`.

- dim:

  Integer dimension vector.

- seed:

  Optional AD value used to create the container. See \`rtmb_vector()\`.

## Value

An AD-compatible array.
