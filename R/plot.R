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
#'
#' @return No return value.
#' @export
plot_forest <- function(x, prob = 0.95) {
  if (is.null(dim(x)) || length(dim(x)) != 3) {
    stop("`x` must be a 3D array with dimension (iterations, chains, variables).", call. = FALSE)
  }

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

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Expand left margin to accommodate parameter names
  par(mar = c(5, 6, 3, 2))

  # Reverse Y-axis to order from bottom to top
  y_pos <- seq(n_variables, 1, by = -1)

  plot(
    x = summaries[, "median"],
    y = y_pos,
    xlim = range(c(summaries[, "lower"], summaries[, "upper"])),
    ylim = c(0.5, n_variables + 0.5),
    type = "n",
    yaxt = "n",
    ylab = "",
    xlab = "Estimate",
    main = sprintf("Forest Plot (%g%% CI)", prob * 100)
  )

  abline(v = 0, lty = 2, col = "gray50")

  # Draw credible interval lines and median points
  segments(
    x0 = summaries[, "lower"], y0 = y_pos,
    x1 = summaries[, "upper"], y1 = y_pos,
    lwd = 2, col = "black"
  )
  points(summaries[, "median"], y_pos, pch = 16, col = "black", cex = 1.2)

  axis(2, at = y_pos, labels = parnames, las = 1, tick = FALSE, line = -0.5)
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
