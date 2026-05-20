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

.rtmb_assignment_names <- function(exprs) {
  find_assignments <- function(x) {
    if (!is.call(x)) return(character())
    if ((identical(x[[1L]], as.name("<-")) || identical(x[[1L]], as.name("="))) &&
        is.name(x[[2L]])) {
      return(as.character(x[[2L]]))
    }
    unique(unlist(lapply(as.list(x)[-1L], find_assignments), use.names = FALSE))
  }
  unique(unlist(lapply(exprs, find_assignments), use.names = FALSE))
}

.rtmb_list_call_add_log_lik <- function(ret_call, log_lik_name = "log_lik") {
  if (is.call(ret_call) && identical(ret_call[[1L]], as.name("list"))) {
    args <- as.list(ret_call)[-1L]
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- rep("", length(args))
    if (!log_lik_name %in% arg_names) {
      args[[length(args) + 1L]] <- as.name(log_lik_name)
      arg_names <- c(arg_names, log_lik_name)
      names(args) <- arg_names
    }
    return(as.call(c(list(as.name("list")), args)))
  }
  NULL
}

.rtmb_waic_generate_ast <- function(existing = NULL, waic = NULL, log_lik_name = "log_lik") {
  waic_exprs <- .rtmb_block_exprs(waic)
  existing_exprs <- .rtmb_block_exprs(existing)
  if (length(waic_exprs) == 0L) return(existing)

  if (length(existing_exprs) == 0L) {
    ret_call <- as.call(c(list(as.name("list")), stats::setNames(list(as.name(log_lik_name)), log_lik_name)))
    return(as.call(c(list(as.name("{")), waic_exprs, list(as.call(list(as.name("return"), ret_call))))))
  }

  last_idx <- length(existing_exprs)
  last_expr <- existing_exprs[[last_idx]]
  if (is.call(last_expr) && identical(last_expr[[1L]], as.name("return"))) {
    ret_expr <- last_expr[[2L]]
    ret_call <- .rtmb_list_call_add_log_lik(ret_expr, log_lik_name)
    if (!is.null(ret_call)) {
      existing_exprs[[last_idx]] <- as.call(list(as.name("return"), ret_call))
      return(as.call(c(list(as.name("{")), existing_exprs[-last_idx], waic_exprs, existing_exprs[last_idx])))
    }
    if (is.name(ret_expr)) {
      add_expr <- as.call(list(as.name("<-"), as.call(list(as.name("$"), ret_expr, as.name(log_lik_name))), as.name(log_lik_name)))
      return(as.call(c(list(as.name("{")), existing_exprs[-last_idx], waic_exprs, list(add_expr), existing_exprs[last_idx])))
    }
  }

  if (is.call(last_expr) && identical(last_expr[[1L]], as.name("list"))) {
    ret_call <- .rtmb_list_call_add_log_lik(last_expr, log_lik_name)
    existing_exprs[[last_idx]] <- as.call(list(as.name("return"), ret_call))
    return(as.call(c(list(as.name("{")), existing_exprs[-last_idx], waic_exprs, existing_exprs[last_idx])))
  }

  if (is.name(last_expr)) {
    add_expr <- as.call(list(as.name("<-"), as.call(list(as.name("$"), last_expr, as.name(log_lik_name))), as.name(log_lik_name)))
    return(as.call(c(list(as.name("{")), existing_exprs[-last_idx], waic_exprs, list(add_expr), list(as.call(list(as.name("return"), last_expr))))))
  }

  returned_names <- .rtmb_assignment_names(existing_exprs)
  ret_args <- lapply(returned_names, as.name)
  names(ret_args) <- returned_names
  ret_args[[length(ret_args) + 1L]] <- as.name(log_lik_name)
  names(ret_args)[length(ret_args)] <- log_lik_name
  ret_call <- as.call(c(list(as.name("list")), ret_args))
  as.call(c(list(as.name("{")), existing_exprs, waic_exprs, list(as.call(list(as.name("return"), ret_call)))))
}
