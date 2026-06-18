#' Plot posterior densities for MCMC samples
#'
#' Draw kernel density plots for each parameter across chains.
#'
#' @param x A 3D array of posterior samples with dimensions
#'   `(iterations, chains, variables)`.
#'   A 2D matrix `(iterations, chains)` is also allowed and is treated
#'   as a single variable.
#' @param mono Logical. If `TRUE`, plot in monochrome. If `FALSE`,
#'   use blue shades.
#'
#' @return No return value. This function is called for its side effect
#'   of plotting.
#'
#' @importFrom graphics plot.new title
#' @export
plot_dens <- function(x, mono = FALSE) {
  if (is.null(dim(x))) {
    stop("`x` must be a matrix or a 3D array.", call. = FALSE)
  }

  dims <- dim(x)

  if (length(dims) == 2) {
    x <- array(x, dim = c(dims[1], dims[2], 1))
    dims <- dim(x)
  }

  if (length(dims) != 3) {
    stop("`x` must have dimension (iterations, chains, variables).", call. = FALSE)
  }

  n_chains <- dims[2]
  n_variables <- dims[3]

  parnames <- dimnames(x)[[3]]
  if (is.null(parnames)) {
    parnames <- paste0("V", seq_len(n_variables))
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  n_cols <- ceiling(sqrt(n_variables))
  n_rows <- ceiling(n_variables / n_cols)

  par(
    mfrow = c(n_rows, n_cols),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  if (isTRUE(mono)) {
    colors <- rep("black", n_chains)
  } else {
    colors <- grDevices::colorRampPalette(c("blue4", "lightblue"))(n_chains)
  }

  for (i in seq_len(n_variables)) {
    dens_list <- vector("list", n_chains)
    max_dens <- 0
    min_x <- Inf
    max_x <- -Inf

    for (n in seq_len(n_chains)) {
      valid_data <- x[, n, i]
      valid_data <- valid_data[is.finite(valid_data)]

      if (length(valid_data) > 1) {
        if (stats::var(valid_data) == 0) {
          valid_data <- jitter(valid_data, amount = 1e-5)
        }
        dens_list[[n]] <- stats::density(valid_data)
        max_dens <- max(max_dens, max(dens_list[[n]]$y))
        min_x <- min(min_x, min(dens_list[[n]]$x))
        max_x <- max(max_x, max(dens_list[[n]]$x))
      }
    }

    valid_idx <- which(!vapply(dens_list, is.null, logical(1)))

    if (length(valid_idx) == 0) {
      plot.new()
      title(main = parnames[i])
      mtext("No finite draws", side = 3, line = -1.5, cex = 0.9)
      next
    }

    plot_idx <- valid_idx[1]

    plot(
      dens_list[[plot_idx]],
      lty = plot_idx,
      col = colors[plot_idx],
      xlab = "",
      ylab = "Density",
      main = parnames[i],
      xlim = c(min_x, max_x),
      ylim = c(0, max_dens * 1.15)
    )

    legend_text <- paste0("chain", plot_idx)
    legend_lty <- plot_idx
    legend_col <- colors[plot_idx]

    if (length(valid_idx) > 1) {
      for (n in valid_idx[-1]) {
        lines(dens_list[[n]], lty = n, col = colors[n])
        legend_text <- c(legend_text, paste0("chain", n))
        legend_lty <- c(legend_lty, n)
        legend_col <- c(legend_col, colors[n])
      }
    }

    legend(
      "topright",
      legend = legend_text,
      lty = legend_lty,
      col = legend_col,
      bg = "transparent",
      bty = "n",
      cex = 0.75
    )
  }

  mtext("Density Plot", outer = TRUE, cex = 1.2, font = 2)
}

#' Plot MCMC trace plots
#'
#' Draw trace plots for each parameter across chains.
#'
#' @param x A 3D array of posterior samples with dimensions
#'   `(iterations, chains, variables)`.
#'   A 2D matrix `(iterations, chains)` is also allowed and is treated
#'   as a single variable.
#' @param mono Logical. If `TRUE`, plot in monochrome. If `FALSE`,
#'   use blue shades.
#'
#' @return No return value. This function is called for its side effect
#'   of plotting.
#'
#' @export
plot_trace <- function(x, mono = FALSE) {
  if (is.null(dim(x))) {
    stop("`x` must be a matrix or a 3D array.", call. = FALSE)
  }

  dims <- dim(x)

  if (length(dims) == 2) {
    x <- array(x, dim = c(dims[1], dims[2], 1))
    dims <- dim(x)
  }

  if (length(dims) != 3) {
    stop("`x` must have dimension (iterations, chains, variables).", call. = FALSE)
  }

  iter <- dims[1]
  n_chains <- dims[2]
  n_variables <- dims[3]

  parnames <- dimnames(x)[[3]]
  if (is.null(parnames)) {
    parnames <- paste0("V", seq_len(n_variables))
  }

  if (isTRUE(mono)) {
    ltys <- seq_len(n_chains)
    colors <- rep("black", n_chains)
  } else {
    ltys <- rep(1, n_chains)
    colors <- grDevices::colorRampPalette(c("blue4", "lightblue"))(n_chains)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Determine grid layout automatically based on the number of variables
  n_cols <- ceiling(sqrt(n_variables))
  n_rows <- ceiling(n_variables / n_cols)

  par(
    mfrow = c(n_rows, n_cols),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (i in seq_len(n_variables)) {
    min_val <- min(x[, , i], na.rm = TRUE)
    max_val <- max(x[, , i], na.rm = TRUE)

    # Add 5% padding to the top and bottom of the Y-axis
    pad <- 0.05 * (max_val - min_val)
    if (is.na(pad) || pad == 0) pad <- 1e-8
    y_limit <- c(min_val - pad, max_val + pad)

    plot(
      seq_len(iter),
      x[, 1, i],
      type = "l",
      ylim = y_limit,
      col = colors[1],
      lty = ltys[1],
      xlab = "Iteration",
      ylab = "Value",
      main = parnames[i]
    )

    legend_text <- "chain1"

    if (n_chains > 1) {
      for (j in 2:n_chains) {
        lines(seq_len(iter), x[, j, i], col = colors[j], lty = ltys[j])
        legend_text <- c(legend_text, paste0("chain", j))
      }
    }

    legend(
      "topright",
      legend = legend_text,
      col = colors,
      lty = ltys,
      bg = "transparent",
      bty = "n",
      cex = 0.75
    )
  }

  mtext("Trace Plot", outer = TRUE, cex = 1.2, font = 2)
}

#' Plot autocorrelation for one variable across chains
#'
#' Draw autocorrelation plots for a selected parameter across all chains.
#'
#' @param x A 3D array of posterior samples with dimensions
#'   `(iterations, chains, variables)`.
#' @param var_idx Integer. Index of the variable to plot.
#'
#' @return No return value. This function is called for its side effect
#'   of plotting.
#'
#' @export
plot_acf <- function(x, var_idx = 1) {
  if (is.null(dim(x)) || length(dim(x)) != 3) {
    stop("`x` must be a 3D array with dimension (iterations, chains, variables).",
         call. = FALSE)
  }

  if (!is.numeric(var_idx) || length(var_idx) != 1) {
    stop("`var_idx` must be a single integer.", call. = FALSE)
  }

  n_variables <- dim(x)[3]
  if (var_idx < 1 || var_idx > n_variables) {
    stop("`var_idx` is out of range.", call. = FALSE)
  }

  parnames <- dimnames(x)[[3]]
  if (is.null(parnames)) {
    parnames <- paste0("V", seq_len(n_variables))
  }

  parname <- parnames[var_idx]
  n_chains <- dim(x)[2]

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(
    mfrow = c(n_chains, 1),
    mar = c(2, 4.1, 1, 2),
    oma = c(1, 0, 2, 0)
  )

  for (i in seq_len(n_chains)) {
    ac_data <- stats::acf(x[, i, var_idx], plot = FALSE)

    plot(
      x = as.vector(ac_data$lag),
      y = as.vector(ac_data$acf),
      type = "h",
      lwd = 2,
      col = "black",
      xlab = "Lag",
      ylab = paste0("chain", i),
      main = ""
    )
    abline(h = 0)
  }

  mtext(paste0("Autocorrelation (", parname, ")"), outer = TRUE, cex = 1.2, font = 2)
}

#' Plot parameter estimates and credible intervals (Forest Plot)
#'
#' @param x A 3D array of posterior samples with dimensions `(iterations, chains, variables)`.
#' @param prob Numeric. Probability mass for the credible interval (default is 0.95).
#' @param point_estimate Character. Point estimate shown as dots. One of
#'   `"median"` (default), `"EAP"`/`"mean"`, `"MAP"`/`"marginal_map"`, or
#'   `"joint_map"`. `"joint_map"` uses the draw with the largest `lp` variable.
#'
#' @return No return value.
#' @export
plot_forest <- function(x, prob = 0.95,
                        point_estimate = c("median", "EAP", "mean", "MAP", "marginal_map", "joint_map")) {
  if (is.null(dim(x)) || length(dim(x)) != 3) {
    stop("`x` must be a 3D array with dimension (iterations, chains, variables).", call. = FALSE)
  }
  point_estimate <- match.arg(point_estimate)

  n_variables <- dim(x)[3]
  parnames <- dimnames(x)[[3]]
  if (is.null(parnames)) {
    parnames <- paste0("V", seq_len(n_variables))
  }

  # Combine all chains and compute summary statistics for each parameter
  lower_prob <- (1 - prob) / 2
  upper_prob <- 1 - lower_prob

  summaries <- t(apply(x, 3, function(v) {
    v_flat <- as.vector(v[!is.na(v)])
    c(
      mean = mean(v_flat),
      median = median(v_flat),
      lower = quantile(v_flat, probs = lower_prob, names = FALSE),
      upper = quantile(v_flat, probs = upper_prob, names = FALSE)
    )
  }))

  point_values <- switch(
    point_estimate,
    median = summaries[, "median"],
    EAP = summaries[, "mean"],
    mean = summaries[, "mean"],
    MAP = apply(x, 3, function(v) .get_marginal_mode(as.vector(v))),
    marginal_map = apply(x, 3, function(v) .get_marginal_mode(as.vector(v))),
    joint_map = {
      lp_idx <- match("lp", parnames)
      if (is.na(lp_idx)) {
        stop("`point_estimate = 'joint_map'` requires a variable named `lp` in `x`.", call. = FALSE)
      }
      lp_vals <- x[, , lp_idx]
      max_idx <- which(lp_vals == max(lp_vals, na.rm = TRUE), arr.ind = TRUE)
      if (nrow(max_idx) == 0L) {
        stop("Failed to find the maximum `lp` value.", call. = FALSE)
      }
      x[max_idx[1, 1], max_idx[1, 2], ]
    }
  )

  point_label <- switch(
    point_estimate,
    median = "median",
    EAP = "EAP",
    mean = "EAP",
    MAP = "MAP",
    marginal_map = "MAP",
    joint_map = "joint MAP"
  )

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Expand left margin to accommodate parameter names
  par(mar = c(5, 6, 3, 2))

  # Reverse Y-axis to order from bottom to top
  y_pos <- seq(n_variables, 1, by = -1)

  plot(
    x = point_values,
    y = y_pos,
    xlim = range(c(summaries[, "lower"], summaries[, "upper"], point_values), na.rm = TRUE),
    ylim = c(0.5, n_variables + 0.5),
    type = "n",
    yaxt = "n",
    ylab = "",
    xlab = "Estimate",
    main = sprintf("Forest Plot (%g%% CI, point = %s)", prob * 100, point_label)
  )

  abline(v = 0, lty = 2, col = "gray50")

  # Draw credible interval lines and selected point estimates
  segments(
    x0 = summaries[, "lower"], y0 = y_pos,
    x1 = summaries[, "upper"], y1 = y_pos,
    lwd = 2, col = "black"
  )
  points(point_values, y_pos, pch = 16, col = "black", cex = 1.2)

  axis(2, at = y_pos, labels = parnames, las = 1, tick = FALSE, line = -0.5)
}

#' Plot Multidimensional Unfolding Configuration
#'
#' @description
#' Draws an item-person map for multidimensional unfolding models. Item
#' locations are shown as labels, optional blue circles show item alpha/radius
#' values, and the respondent locations are summarized by a two-
#' dimensional kernel density contour.
#'
#' @param delta Item coordinate matrix (M items x D dimensions), or a fitted
#'   object with an `EAP()` method.
#' @param theta Person coordinate matrix (N persons x D dimensions). If `delta`
#'   is a fitted object and `theta` is `NULL`, it is extracted from the fit.
#' @param item_alpha Optional item alpha vector used for item radii. If `delta`
#'   is a fitted object and `item_alpha` is `NULL`, `alpha` is extracted when
#'   available.
#' @param phi Deprecated alias for `item_alpha`.
#' @param dims Integer vector of length 2 specifying dimensions to plot.
#' @param radius Numeric plot radius. If `NULL`, a radius is chosen from the
#'   coordinates.
#' @param signs Numeric vector of length 2. Use `-1` to flip an axis.
#' @param item_labels Optional item labels. Defaults to row names or item
#'   numbers.
#' @param show_radius Logical; whether to draw radius circles based on item alpha.
#' @param show_density Logical; whether to draw density contours for `theta`.
#' @param circle_scale Numeric multiplier for item radii.
#' @param alpha Circle transparency.
#' @param contour_n Grid size passed to `MASS::kde2d()`.
#' @param distance Character; `"auto"`, `"squared"`, or `"euclidean"`. Used
#'   to transform item alpha values into plotted radii. `"auto"` uses the fit's
#'   stored distance when available.
#' @param point_estimate Character; point estimate used when `delta` is a fit
#'   object. Passed to `estimate()`.
#' @param prefer_rotated Logical; when `delta` is a fitted object, prefer
#'   rotated generated quantities (`delta_rot` and `theta_rot`) if available.
#' @param show_phi Deprecated alias for `show_radius`.
#' @param main,xlab,ylab Plot title and axis labels.
#' @param ... Additional arguments passed to `plot()`.
#'
#' @return Invisibly returns a list with plotted `delta`, `theta`, and
#'   `item_alpha`.
#' @export
plot_mdu <- function(delta, theta = NULL, item_alpha = NULL, phi = NULL,
                     dims = c(1, 2), radius = NULL, signs = c(1, 1),
                     item_labels = NULL, show_radius = TRUE, show_density = TRUE,
                     circle_scale = 1, alpha = 0.2, contour_n = 60,
                     distance = c("auto", "squared", "euclidean"),
                     point_estimate = c("EAP", "MAP", "mean", "marginal_map", "joint_map"),
                     prefer_rotated = TRUE,
                     show_phi = NULL,
                     main = "MDU Configuration", xlab = NULL, ylab = NULL, ...) {
  point_estimate <- match.arg(point_estimate)
  distance <- match.arg(distance)
  prefer_rotated <- isTRUE(prefer_rotated)
  phi_supplied <- !is.null(phi)

  if (!is.null(show_phi)) {
    warning("'show_phi' is deprecated; use 'show_radius' instead.", call. = FALSE)
    show_radius <- show_phi
  }

  if (!is.matrix(delta) && is.function(delta$estimate)) {
    fit <- delta
    if (distance == "auto") {
      fit_distance <- fit$extra$distance %||% fit$model$extra$distance %||% NULL
      fit_distance <- if (length(fit_distance)) as.character(fit_distance[1]) else NULL
      distance <- if (!is.null(fit_distance) && fit_distance %in% c("squared", "euclidean")) {
        fit_distance
      } else {
        "squared"
      }
    }
    est_type <- switch(
      point_estimate,
      EAP = "EAP",
      mean = "mean",
      MAP = "MAP",
      marginal_map = "marginal_map",
      joint_map = "joint_map"
    )
    est <- fit$estimate(type = est_type, pars = "all", drop = FALSE)
    use_rotated <- prefer_rotated && !is.null(est$delta_rot)
    if (use_rotated) {
      delta <- est$delta_rot
      if (is.null(theta)) {
        if (!is.null(est$theta_rot)) {
          theta <- est$theta_rot
        } else if (!is.null(est$theta)) {
          warning(
            "Using 'delta_rot' but 'theta_rot' was not found. ",
            "Use rotate(..., linked = 'theta') to rotate person coordinates as well.",
            call. = FALSE
          )
          theta <- est$theta
        }
      }
    } else if (is.null(est$delta)) {
      stop("The fitted object does not contain 'delta'.", call. = FALSE)
    } else {
      delta <- est$delta
      if (is.null(theta)) theta <- est$theta
    }
    if (is.null(item_alpha)) item_alpha <- est$alpha
    if (is.null(item_alpha) && is.null(phi)) phi <- est$phi
  }
  if (distance == "auto") distance <- "squared"
  if (is.null(item_alpha) && !is.null(phi)) {
    if (phi_supplied) {
      warning("'phi' is deprecated; use 'item_alpha' instead.", call. = FALSE)
    }
    item_alpha <- phi
  }

  delta <- as.matrix(delta)
  if (is.null(theta)) {
    stop("'theta' must be supplied unless 'delta' is a fitted object containing theta.", call. = FALSE)
  }
  theta <- as.matrix(theta)

  D <- ncol(delta)
  if (length(dims) != 2L || any(dims < 1) || any(dims > D)) {
    stop("'dims' must be an integer vector of length 2 within the coordinate dimensions.", call. = FALSE)
  }
  if (ncol(theta) < max(dims)) {
    stop("'theta' does not have enough columns for the requested 'dims'.", call. = FALSE)
  }
  if (length(signs) != 2L || any(!signs %in% c(-1, 1))) {
    stop("'signs' must be a numeric vector of length 2 containing only -1 or 1.", call. = FALSE)
  }

  x_delta <- signs[1] * delta[, dims[1]]
  y_delta <- signs[2] * delta[, dims[2]]
  x_theta <- signs[1] * theta[, dims[1]]
  y_theta <- signs[2] * theta[, dims[2]]

  if (is.null(item_labels)) {
    item_labels <- rownames(delta)
    if (is.null(item_labels)) item_labels <- seq_len(nrow(delta))
    item_labels <- sub("^[Ii]tem\\s*([0-9]+)$", "\\1", as.character(item_labels))
  }
  if (length(item_labels) != nrow(delta)) {
    stop("'item_labels' must have the same length as nrow(delta).", call. = FALSE)
  }

  finite_coords <- c(x_delta, y_delta, x_theta, y_theta)
  finite_coords <- finite_coords[is.finite(finite_coords)]
  if (is.null(radius)) {
    radius <- max(abs(finite_coords), na.rm = TRUE)
    if (!is.finite(radius) || radius <= 0) radius <- 1
    radius <- radius * 1.1
  }

  if (is.null(xlab)) xlab <- paste0("Dim ", dims[1])
  if (is.null(ylab)) ylab <- paste0("Dim ", dims[2])

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(pty = "s")

  plot(
    c(-radius, radius), c(-radius, radius),
    type = "n", asp = 1, xlab = xlab, ylab = ylab, main = main, ...
  )
  abline(h = 0, v = 0, col = "gray85", lty = 3)

  if (isTRUE(show_density)) {
    ok <- is.finite(x_theta) & is.finite(y_theta)
    if (sum(ok) > 5 && stats::sd(x_theta[ok]) > 0 && stats::sd(y_theta[ok]) > 0) {
      knl <- MASS::kde2d(
        x_theta[ok], y_theta[ok],
        n = contour_n,
        lims = c(-radius, radius, -radius, radius)
      )
      graphics::contour(knl, add = TRUE, drawlabels = TRUE, col = "gray25", lwd = 0.8)
    }
  }

  if (isTRUE(show_radius) && !is.null(item_alpha)) {
    item_alpha <- as.numeric(item_alpha)
    if (length(item_alpha) == 1L) {
      item_radius <- rep(0, nrow(delta))
    } else if (length(item_alpha) == nrow(delta)) {
      item_radius <- pmax(item_alpha - min(item_alpha, na.rm = TRUE), 0)
      if (distance == "squared") item_radius <- sqrt(item_radius)
    } else {
      stop("'item_alpha' must have length 1 or the same length as nrow(delta).", call. = FALSE)
    }
    max_item_radius <- max(item_radius, na.rm = TRUE)
    if (is.finite(max_item_radius) && max_item_radius > 0) {
      item_radius <- item_radius * circle_scale
      circle_col <- grDevices::rgb(0.1, 0.25, 0.9, alpha = alpha)
      angle <- seq(0, 2 * pi, length.out = 160)
      for (m in seq_along(item_radius)) {
        if (is.finite(item_radius[m]) && item_radius[m] > 0) {
          lines(
            x_delta[m] + item_radius[m] * cos(angle),
            y_delta[m] + item_radius[m] * sin(angle),
            col = circle_col
          )
        }
      }
    }
  }

  points(x_theta, y_theta, pch = 16, col = grDevices::rgb(0, 0, 0, 0.12), cex = 0.45)
  text(x_delta, y_delta, labels = item_labels, font = 2, cex = 0.9)

  invisible(list(
    delta = cbind(x_delta, y_delta),
    theta = cbind(x_theta, y_theta),
    item_alpha = item_alpha,
    item_labels = item_labels
  ))
}

#' Plot pairs for posterior samples
#'
#' Draw a scatterplot matrix to examine posterior correlations between parameters.
#'
#' @param x A 3D array of posterior samples with dimensions `(iterations, chains, variables)`.
#' @param pars Character vector or integer vector specifying which parameters to plot.
#'   If NULL (default), up to the first 10 parameters are plotted to prevent overplotting.
#'
#' @return No return value.
#' @export
plot_pairs <- function(x, pars = NULL) {
  if (is.null(dim(x)) || length(dim(x)) != 3) {
    stop("`x` must be a 3D array with dimension (iterations, chains, variables).", call. = FALSE)
  }

  n_variables <- dim(x)[3]
  parnames <- dimnames(x)[[3]]
  if (is.null(parnames)) {
    parnames <- paste0("V", seq_len(n_variables))
  }

  # Filter target parameters for display
  if (!is.null(pars)) {
    if (is.character(pars)) {
      idx <- match(pars, parnames)
      if (any(is.na(idx))) stop("Some specified parameters were not found.", call. = FALSE)
      plot_idx <- idx
    } else if (is.numeric(pars)) {
      plot_idx <- pars
    }
  } else {
    # Default to maximum of 10 parameters (more would be slow and cluttered)
    plot_idx <- seq_len(min(n_variables, 10))
  }

  if (length(plot_idx) < 2) {
    stop("At least 2 parameters are required for a pairs plot.", call. = FALSE)
  }

  # Combine all chains into a 2D matrix
  mat <- apply(x[, , plot_idx, drop = FALSE], 3, as.vector)
  colnames(mat) <- parnames[plot_idx]

  # Custom function to display correlation coefficients in the upper panels
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = "complete.obs")
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    # Adjust text size based on the strength of the correlation
    text(0.5, 0.5, txt, cex = cex.cor * (0.5 + abs(r) / 2))
  }

  # Custom function to draw semi-transparent points to visualize density
  panel.scatter <- function(x, y, ...) {
    points(x, y, pch = 16, col = rgb(0, 0, 0, alpha = 0.2), cex = 0.5)
  }

  pairs(mat,
        lower.panel = panel.scatter,
        upper.panel = panel.cor,
        gap = 0.2,
        main = "Pairs Plot")
}

#' Plot least-squares marginal means
#'
#' Calculate and plot marginal means for a fitted model.
#'
#' @param fit A fitted model object (e.g. `Classic_Fit`).
#' @param specs Character vector specifying the variables for which marginal means are calculated.
#' @param ... Additional arguments passed to \code{\link{lsmeans}} or the plot method.
#'
#' @return Invisibly returns the plotted object.
#' @export
plot_lsmeans <- function(fit, specs, ...) {
  res <- lsmeans(fit, specs, ...)
  plot(res, ...)
}

#' Plot conditional effects
#'
#' Calculate and plot conditional effects for a fitted model.
#'
#' @param fit A fitted model object (e.g., `MCMC_Fit`, `MAP_Fit`, `VB_Fit`).
#' @param effect Name of the variable to visualize (e.g., "X1" or "X1:X2").
#' @param prob Probability for the credible/confidence interval.
#' @param sd_multiplier Numeric multiplier for standard deviation when splitting continuous moderators.
#' @param sd_slice Logical or NULL; controls whether continuous moderators are
#'   evaluated at mean - SD, mean, and mean + SD.
#' @param resolution Grid resolution to calculate for continuous variables.
#' @param ... Additional arguments passed to the plot method.
#'
#' @return Invisibly returns the plotted object of class \code{ce_rtmb}.
#' @export
plot_conditional_effects <- function(fit, effect, prob = 0.95, sd_multiplier = 1,
                                     sd_slice = NULL, resolution = 100, ...) {
  res <- conditional_effects(
    fit, effect,
    prob = prob,
    sd_multiplier = sd_multiplier,
    sd_slice = sd_slice,
    resolution = resolution
  )
  plot(res, ...)
}

#' Plot item/category response curves
#'
#' Calculate and plot item/category response curves for a fitted IRT model.
#'
#' @param fit A fitted IRT model object.
#' @param ... Additional arguments passed to \code{\link{item_curve}} or the plot method.
#'
#' @return Invisibly returns the plotted object of class \code{rtmb_item_curve}.
#' @export
plot_item_curve <- function(fit, ...) {
  res <- item_curve(fit, ...)
  plot(res, ...)
}

#' Plot item information functions
#'
#' Calculate and plot item information functions for a fitted IRT model.
#'
#' @param fit A fitted IRT model object.
#' @param ... Additional arguments passed to \code{\link{item_info}} or the plot method.
#'
#' @return Invisibly returns the plotted object of class \code{rtmb_item_info}.
#' @export
plot_item_info <- function(fit, ...) {
  res <- item_info(fit, ...)
  plot(res, ...)
}

#' Plot test information function
#'
#' Calculate and plot test information function for a fitted IRT model.
#'
#' @param fit A fitted IRT model object.
#' @param ... Additional arguments passed to \code{\link{test_info}} or the plot method.
#'
#' @return Invisibly returns the plotted object of class \code{rtmb_test_info}.
#' @export
plot_test_info <- function(fit, ...) {
  res <- test_info(fit, ...)
  plot(res, ...)
}
