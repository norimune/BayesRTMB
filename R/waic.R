.log_mean_exp <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  m <- max(x)
  m + log(mean(exp(x - m)))
}

.as_log_lik_matrix <- function(x) {
  if (length(dim(x)) == 3L) {
    dn <- dim(x)
    matrix(x, nrow = dn[1L] * dn[2L], ncol = dn[3L])
  } else {
    mat <- as.matrix(x)
    if (ncol(mat) == 0L) stop("No pointwise log_lik values were found.", call. = FALSE)
    mat
  }
}

.compute_waic_from_log_lik <- function(log_lik_draws) {
  mat <- .as_log_lik_matrix(log_lik_draws)
  if (nrow(mat) < 2L) {
    stop(
      "WAIC requires more than one posterior or approximate posterior draw of pointwise log_lik.",
      call. = FALSE
    )
  }

  lppd_i <- apply(mat, 2L, .log_mean_exp)
  p_waic_i <- apply(mat, 2L, stats::var, na.rm = TRUE)
  p_waic_i[!is.finite(p_waic_i)] <- NA_real_
  elpd_i <- lppd_i - p_waic_i

  pointwise <- data.frame(
    lppd = lppd_i,
    p_waic = p_waic_i,
    WAIC = elpd_i,
    `-2 WAIC` = -2 * elpd_i,
    check.names = FALSE
  )

  out <- list(
    WAIC = sum(elpd_i, na.rm = TRUE),
    minus2_WAIC = -2 * sum(elpd_i, na.rm = TRUE),
    lppd = sum(lppd_i, na.rm = TRUE),
    p_waic = sum(p_waic_i, na.rm = TRUE),
    pointwise = pointwise,
    nobs = ncol(mat),
    ndraws = nrow(mat)
  )
  class(out) <- "waic_BayesRTMB"
  out
}

#' @export
print.waic_BayesRTMB <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("WAIC\n")
  cat(sprintf("  WAIC (elpd scale): %.*f\n", digits, x$WAIC))
  cat(sprintf("  -2 WAIC: %.*f\n", digits, x$minus2_WAIC))
  cat(sprintf("  p_waic: %.*f\n", digits, x$p_waic))
  cat(sprintf("  lppd: %.*f\n", digits, x$lppd))
  cat(sprintf("  nobs: %d, draws: %d\n", x$nobs, x$ndraws))
  invisible(x)
}

.waic_from_fit_draws <- function(fit, ...) {
  arr <- tryCatch(
    fit$draws(pars = "log_lik", inc_random = FALSE, inc_transform = FALSE, inc_generate = TRUE, ...),
    error = function(e) NULL
  )
  if (is.null(arr) || length(arr) == 0L || dim(arr)[3L] == 0L) {
    stop(
      "WAIC requires pointwise generated quantity `log_lik`. ",
      "Create the model with WAIC = TRUE or add `log_lik` with generated_quantities().",
      call. = FALSE
    )
  }
  .compute_waic_from_log_lik(arr)
}

.rtmb_block_exprs <- function(ast) {
  if (is.null(ast)) return(list())
  if (is.call(ast) && identical(ast[[1L]], as.name("{"))) as.list(ast)[-1L] else list(ast)
}

.rtmb_merge_generate_ast <- function(existing = NULL, extra = NULL) {
  exprs <- c(.rtmb_block_exprs(existing), .rtmb_block_exprs(extra))
  if (length(exprs) == 0L) return(NULL)
  as.call(c(list(as.name("{")), exprs))
}
