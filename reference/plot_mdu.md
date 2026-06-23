# Plot Multidimensional Unfolding Configuration

Draws an item-person map for multidimensional unfolding models. Item
locations are shown as labels, optional blue circles show item
alpha/radius values, and the respondent locations are summarized by a
two- dimensional kernel density contour.

## Usage

``` r
plot_mdu(
  delta,
  theta = NULL,
  item_alpha = NULL,
  phi = NULL,
  dims = c(1, 2),
  radius = NULL,
  signs = c(1, 1),
  item_labels = NULL,
  show_radius = TRUE,
  show_density = TRUE,
  circle_scale = 1,
  alpha = 0.2,
  contour_n = 60,
  distance = c("auto", "squared", "euclidean"),
  point_estimate = c("EAP", "MAP", "mean", "marginal_map", "joint_map"),
  prefer_rotated = TRUE,
  show_phi = NULL,
  main = "MDU Configuration",
  xlab = NULL,
  ylab = NULL,
  ...
)
```

## Arguments

- delta:

  Item coordinate matrix (M items x D dimensions), or a fitted object
  with an \`EAP()\` method.

- theta:

  Person coordinate matrix (N persons x D dimensions). If \`delta\` is a
  fitted object and \`theta\` is \`NULL\`, it is extracted from the fit.

- item_alpha:

  Optional item alpha vector used for item radii. If \`delta\` is a
  fitted object and \`item_alpha\` is \`NULL\`, \`alpha\` is extracted
  when available.

- phi:

  Deprecated alias for \`item_alpha\`.

- dims:

  Integer vector of length 2 specifying dimensions to plot.

- radius:

  Numeric plot radius. If \`NULL\`, a radius is chosen from the
  coordinates.

- signs:

  Numeric vector of length 2. Use \`-1\` to flip an axis.

- item_labels:

  Optional item labels. Defaults to row names or item numbers.

- show_radius:

  Logical; whether to draw radius circles based on item alpha.

- show_density:

  Logical; whether to draw density contours for \`theta\`.

- circle_scale:

  Numeric multiplier for item radii.

- alpha:

  Circle transparency.

- contour_n:

  Grid size passed to \`MASS::kde2d()\`.

- distance:

  Character; \`"auto"\`, \`"squared"\`, or \`"euclidean"\`. Used to
  transform item alpha values into plotted radii. \`"auto"\` uses the
  fit's stored distance when available.

- point_estimate:

  Character; point estimate used when \`delta\` is a fit object. Passed
  to \`estimate()\`.

- prefer_rotated:

  Logical; when \`delta\` is a fitted object, prefer rotated generated
  quantities (\`delta_rot\` and \`theta_rot\`) if available.

- show_phi:

  Deprecated alias for \`show_radius\`.

- main, xlab, ylab:

  Plot title and axis labels.

- ...:

  Additional arguments passed to \`plot()\`.

## Value

Invisibly returns a list with plotted \`delta\`, \`theta\`, and
\`item_alpha\`.
