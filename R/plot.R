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

  par(
    mfrow = c(1, n_variables),
    mar = c(5.1, 4.1, 1, 1),
    oma = c(0, 0, 2, 0)
  )

  if (isTRUE(mono)) {
    colors <- rep("black", n_chains)
  } else {
    colors <- grDevices::colorRampPalette(c("blue4", "lightblue"))(n_chains)
  }

  for (i in seq_len(n_variables)) {
    min_data <- min(x[, , i], na.rm = TRUE)
    max_data <- max(x[, , i], na.rm = TRUE)

    dens_list <- vector("list", n_chains)
    max_dens <- 0

    for (n in seq_len(n_chains)) {
      dens_list[[n]] <- stats::density(x[, n, i], na.rm = TRUE)
      max_dens <- max(max_dens, max(dens_list[[n]]$y, na.rm = TRUE))
    }

    plot(
      dens_list[[1]],
      lty = 1,
      col = colors[1],
      xlab = parnames[i],
      ylab = "",
      main = "",
      xlim = c(min_data, max_data),
      ylim = c(0, max_dens + 0.15 * max_dens)
    )

    legend_text <- "chain1"

    if (n_chains > 1) {
      for (n in 2:n_chains) {
        lines(dens_list[[n]], lty = n, col = colors[n])
        legend_text <- c(legend_text, paste0("chain", n))
      }
    }

    legend(
      "topright",
      legend = legend_text,
      lty = seq_len(n_chains),
      col = colors,
      bg = "transparent",
      bty = "n",
      cex = 0.75
    )
  }

  mtext("Density", outer = TRUE)
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

  par(
    mfrow = c(1, n_variables),
    mar = c(5.1, 4.1, 1, 1),
    oma = c(0, 0, 2, 0)
  )

  for (i in seq_len(n_variables)) {
    min_val <- min(x[, , i], na.rm = TRUE)
    max_val <- max(x[, , i], na.rm = TRUE)

    pad <- 0.1 * (max_val - min_val)
    if (pad == 0) pad <- 1e-8
    y_limit <- c(min_val, max_val + pad)

    plot(
      seq_len(iter),
      x[, 1, i],
      type = "l",
      ylim = y_limit,
      col = colors[1],
      lty = ltys[1],
      xlab = "Iteration",
      ylab = "",
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

  mtext("Trace Plot", outer = TRUE)
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

  mtext(paste0("Autocorrelation (", parname, ")"), outer = TRUE)
}
