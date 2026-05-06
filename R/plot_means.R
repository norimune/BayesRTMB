#' Plot marginal means with error bars
#'
#' plot_means function creates a bar plot of marginal means with error bars.
#'
#' @param object An rtmb_lsmeans object.
#' @param error_bar Type of error bar.
#' @param col Color of the bars.
#' @param main Plot title.
#' @param xlab X axis label.
#' @param ylab Y axis label.
#' @param ... Additional arguments.
#' @export
plot_means <- function(object, error_bar = "se", col = NULL, main = NULL, xlab = NULL, ylab = NULL, ...) {
  if (!inherits(object, "rtmb_lsmeans")) {
    if (!all(c("estimate", "Std. Error") %in% names(object))) {
      stop("object must be an 'rtmb_lsmeans' object or a dataframe with 'estimate' and 'Std. Error'.")
    }
  }
  
  if (length(error_bar) > 1) error_bar <- error_bar[1]
  error_bar <- match.arg(error_bar, c("se", "ci"))
  means <- object$estimate
  names(means) <- rownames(object)
  
  if (error_bar == "se") {
    upper <- means + object$`Std. Error`
    lower <- means - object$`Std. Error`
    if (is.null(ylab)) ylab <- "Mean SE"
  } else {
    upper <- object$`Upper 95%`
    lower <- object$`Lower 95%`
    if (is.null(ylab)) ylab <- "Mean CI"
  }
  
  # Styling defaults
  if (is.null(col)) col <- grDevices::hcl.colors(length(means), "Pastel 1")
  if (is.null(main)) main <- "Marginal Means"
  
  # Handle limits
  ylim <- range(c(0, lower, upper), na.rm = TRUE)
  if (all(means >= 0)) ylim[1] <- 0
  if (all(means <= 0)) ylim[2] <- 0
  
  # Add some padding to ylim
  pad <- diff(ylim) * 0.1
  if (ylim[2] > 0) ylim[2] <- ylim[2] + pad
  if (ylim[1] < 0) ylim[1] <- ylim[1] - pad

  # Draw barplot
  old_par <- par(mar = c(5, 5, 4, 2) + 0.1, mgp = c(3, 0.7, 0), tcl = -0.3)
  on.exit(par(old_par))
  
  b <- barplot(means, ylim = ylim, col = col, main = main, xlab = xlab, ylab = ylab, 
               border = "white", las = 1, cex.axis = 0.9, ...)
  
  # Add a light grid
  abline(h = 0, lwd = 1.5, col = "gray40")
  grid(nx = NA, ny = NULL, col = "gray90", lty = "solid")
  # Redraw bars over grid
  barplot(means, add = TRUE, col = col, border = "white", axes = FALSE, ...)

  # Add error bars
  segments(b, lower, b, upper, lwd = 1.8, col = "gray20")
  # Caps
  eps <- 0.05 * diff(par("usr")[1:2]) / length(means)
  segments(b - eps, lower, b + eps, lower, lwd = 1.8, col = "gray20")
  segments(b - eps, upper, b + eps, upper, lwd = 1.8, col = "gray20")
  
  invisible(b)
}

#' @export
plot.rtmb_lsmeans <- function(x, ...) {
  plot_means(x, ...)
}
