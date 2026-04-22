# Code block for parameter definitions

Defines the list of parameters to be passed to \`rtmb_model\`. By
writing in the format \`variable_name = Dim(...)\` within a block
enclosed in \`\`, evaluation is delayed, enabling strict error checking
during model construction.

## Usage

``` r
parameters_code(expr)
```

## Arguments

- expr:

  A parameter definition expression enclosed in \`\`.

## Value

A lazily evaluated list of parameters.
