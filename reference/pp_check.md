# Posterior predictive checks

Compare observed data with replicated outcomes from an estimated
BayesRTMB model. Density overlays are used for continuous outcomes and
binned bar plots for discrete outcomes. When \`type = "auto"\`, the
display is selected from the likelihood's \`\_lpdf\` or \`\_lpmf\`
implementation.

## Usage

``` r
pp_check(object, ...)

# S3 method for class 'RTMB_Fit_Base'
pp_check(object, ...)

# S3 method for class 'rtmb_pp_check'
plot(
  x,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  observed_col = "#1B1B1B",
  predictive_col = "#2C7FB8",
  show_legend = TRUE,
  legend_position = c("auto", "topleft", "topright", "bottomleft", "bottomright"),
  legend_cex = 0.78,
  interval = 0.95,
  ...
)
```

## Arguments

- object:

  A BayesRTMB fit object.

- ...:

  Arguments passed to the fit object's \`pp_check()\` method or to the
  plotting method.

- x:

  An \`rtmb_pp_check\` object.

- main:

  Optional plot title.

- xlab:

  Optional x-axis label.

- ylab:

  Optional y-axis label.

- observed_col:

  Color for observed data.

- predictive_col:

  Color for predictive data and intervals.

- show_legend:

  Logical; display the plot legend.

- legend_position:

  Legend position. \`"auto"\` moves scatter-plot legends away from the
  fitted trend; the other values are passed to \`legend()\`.

- legend_cex:

  Relative text and symbol size for the legend.

- interval:

  Probability covered by posterior predictive intervals.

## Value

An object of class \`rtmb_pp_check\`, returned invisibly after plotting.

## Details

Supplying \`x\` to the fit object's \`pp_check()\` method switches to a
scatter-based check. Observed outcomes are compared with posterior
predictive means and 95 value \`x = ".fitted"\` uses the posterior
predictive mean on the horizontal axis for a model-wide calibration
check.
