# Posterior predictive simulation

Generate replicated outcomes from an estimated BayesRTMB model.
Regression wrappers provide an automatic simulator. For a custom model,
supply an \`rtmb_code()\` object containing a \`generate\` block that
creates and reports one replicated outcome vector. MCMC and variational
fits use posterior draws. A MAP fit uses sampling-based uncertainty when
available and otherwise conditions on its point estimate.

## Usage

``` r
posterior_predict(object, ...)

# S3 method for class 'RTMB_Fit_Base'
posterior_predict(object, ...)
```

## Arguments

- object:

  A BayesRTMB fit object.

- ...:

  Arguments passed to the fit object's \`posterior_predict()\` method.

## Value

A numeric matrix with posterior predictive draws in rows and
observations in columns.
