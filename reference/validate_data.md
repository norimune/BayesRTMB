# Pre-validation of data and parameters

Before passing data to \`rtmb_model\`, it checks for the presence of
R-specific data types (such as data.frame) and missing values, and
outputs errors or warnings accordingly.

## Usage

``` r
validate_data(dat_list)
```

## Arguments

- dat_list:

  A list of data to be passed to the model (usually containing matrices
  or numeric vectors).

## Value

Returns invisible \`NULL\` on success. Interrupts execution with
\`stop()\` or issues \`warning()\` if issues are found.

## Details

RTMB's automatic differentiation engine requires pure \`matrix\` or
\`numeric\` types. If a user mistakenly passes a \`data.frame\`,
incomprehensible errors occur during the construction of the computation
graph. This function catches such issues in advance.
