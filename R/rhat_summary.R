#' Summarize MCMC R-hat Values
#'
#' Extract R-hat values from a fitted MCMC object as a numeric vector.
#' The returned object prints a compact convergence summary and can be passed
#' directly to functions such as \code{\link[graphics]{hist}}.
#'
#' @param fit A fitted model object.
#' @param pars Character or numeric vector specifying parameters to include.
#'   If \code{NULL}, all available parameters are summarized.
#' @param chains Numeric vector specifying chains to include.
#' @param best_chains Integer; number of best chains to retain based on mean log-posterior.
#' @param inc_random Logical; whether to include random effects. Default is \code{FALSE}.
#' @param inc_transform Logical; whether to include transformed parameters. Default is \code{TRUE}.
#' @param inc_generate Logical; whether to include generated quantities. Default is \code{TRUE}.
#' @param finite Logical; whether to drop non-finite or missing R-hat values. Default is \code{TRUE}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A numeric vector of R-hat values with class \code{"rhat_summary"}.
#' @export
#'
#' @examples
#' \donttest{
#'   data(debate, package = "BayesRTMB")
#'   mdl <- rtmb_lm(sat ~ talk + perf, data = debate)
#'   fit <- mdl$sample(sampling = 200, warmup = 200, chains = 2)
#'   r <- rhat_summary(fit)
#'   r
#'   hist(r)
#' }
rhat_summary <- function(fit, ...) {
  UseMethod("rhat_summary")
}

#' @export
rhat_summary.default <- function(fit, ...) {
  stop(sprintf("No rhat_summary method for object of class '%s'.", class(fit)[1]), call. = FALSE)
}

#' @export
#' @method rhat_summary mcmc_fit
rhat_summary.mcmc_fit <- function(fit, pars = NULL,
                                  chains = NULL,
                                  best_chains = NULL,
                                  inc_random = FALSE,
                                  inc_transform = TRUE,
                                  inc_generate = TRUE,
                                  finite = TRUE, ...) {
  summ <- fit$summary(
    pars = pars,
    chains = chains,
    best_chains = best_chains,
    max_rows = NULL,
    inc_random = inc_random,
    inc_transform = inc_transform,
    inc_generate = inc_generate
  )
  if (!("rhat" %in% names(summ))) {
    stop("The fit summary does not contain an R-hat column.", call. = FALSE)
  }

  out <- as.numeric(summ$rhat)
  names(out) <- summ$variable
  if (isTRUE(finite)) out <- out[is.finite(out)]
  class(out) <- c("rhat_summary", "numeric")
  attr(out, "thresholds") <- c(warning = 1.01, problem = 1.05)
  out
}

#' @export
print.rhat_summary <- function(x, digits = 3, ...) {
  r <- unclass(x)
  r <- r[is.finite(r)]
  thresholds <- attr(x, "thresholds")
  warning_threshold <- thresholds[["warning"]] %||% 1.01
  problem_threshold <- thresholds[["problem"]] %||% 1.05

  cat("R-hat summary\n")
  if (length(r) == 0L) {
    cat("No finite R-hat values.\n")
    return(invisible(x))
  }

  qs <- stats::quantile(r, probs = c(0, 0.5, 0.9, 0.95, 1), names = FALSE, na.rm = TRUE)
  out <- data.frame(
    n = length(r),
    min = qs[1],
    median = qs[2],
    q90 = qs[3],
    q95 = qs[4],
    max = qs[5],
    n_gt_1.01 = sum(r > warning_threshold),
    n_gt_1.05 = sum(r > problem_threshold),
    check.names = FALSE
  )
  print(format(out, digits = digits), row.names = FALSE, quote = FALSE)
  invisible(x)
}

#' @export
hist.rhat_summary <- function(x, ..., main = "R-hat values", xlab = "R-hat") {
  graphics::hist(unclass(x), ..., main = main, xlab = xlab)
}
