#' Summarize MCMC Draws Stored in an Array
#'
#' Summarizes posterior draws stored in the same array layout used by
#' \code{\link{MCMC_Fit}}: iterations by chains by variables. The returned
#' object uses the same columns and print method as
#' \code{MCMC_Fit$summary()}.
#'
#' @param draws A numeric three-dimensional array with dimensions
#'   iterations by chains by variables.
#' @param pars Optional numeric or character vector selecting variables.
#'   Character values may be full variable names (for example, `"beta[1]"`) or
#'   base names (for example, `"beta"`). Prefix character names with `"-"` to
#'   exclude them.
#' @param chains Optional numeric vector selecting chains. Positive and negative
#'   integer indexing are supported.
#' @param max_rows Maximum number of variables to include. Use `NULL` to include
#'   all variables.
#' @param digits Number of decimal places used when printing the result.
#'
#' @return A data frame with class `"summary_BayesRTMB"` containing posterior
#'   means, standard deviations, marginal MAP estimates, 95 percent intervals,
#'   bulk and tail effective sample sizes, and split R-hat values.
#'
#' @examples
#' set.seed(123)
#' draws <- array(
#'   rnorm(200 * 2 * 2),
#'   dim = c(200, 2, 2),
#'   dimnames = list(NULL, c("chain1", "chain2"), c("alpha", "beta"))
#' )
#' summary_mcmc(draws)
#'
#' @export
summary_mcmc <- function(draws, pars = NULL, chains = NULL,
                         max_rows = 10, digits = 2) {
  if (!is.array(draws) || length(dim(draws)) != 3L) {
    stop(
      "'draws' must be a three-dimensional array with dimensions ",
      "iterations x chains x variables.",
      call. = FALSE
    )
  }
  if (!is.numeric(draws)) {
    stop("'draws' must be numeric.", call. = FALSE)
  }
  if (any(dim(draws) == 0L)) {
    stop("'draws' must have at least one iteration, chain, and variable.",
         call. = FALSE)
  }

  if (!is.null(max_rows)) {
    if (length(max_rows) != 1L || !is.numeric(max_rows) ||
        !is.finite(max_rows) || max_rows < 1 || max_rows != as.integer(max_rows)) {
      stop("'max_rows' must be NULL or a positive integer.", call. = FALSE)
    }
    max_rows <- as.integer(max_rows)
  }
  if (length(digits) != 1L || !is.numeric(digits) ||
      !is.finite(digits) || digits < 0 || digits != as.integer(digits)) {
    stop("'digits' must be a non-negative integer.", call. = FALSE)
  }
  digits <- as.integer(digits)

  n_chains <- dim(draws)[2]
  chain_idx <- seq_len(n_chains)
  if (!is.null(chains)) {
    if (!is.numeric(chains) || anyNA(chains) || any(!is.finite(chains)) ||
        any(chains != as.integer(chains))) {
      stop("'chains' must be a numeric vector of integer indices.",
           call. = FALSE)
    }
    chain_idx <- tryCatch(
      seq_len(n_chains)[as.integer(chains)],
      error = function(e) NULL
    )
    if (is.null(chain_idx) || length(chain_idx) == 0L || anyNA(chain_idx)) {
      stop("The specified chains were not found.", call. = FALSE)
    }
  }
  draws <- draws[, chain_idx, , drop = FALSE]

  n_vars <- dim(draws)[3]
  param_names <- dimnames(draws)[[3]]
  if (is.null(param_names)) {
    param_names <- paste0("V", seq_len(n_vars))
  } else {
    missing_names <- is.na(param_names) | !nzchar(param_names)
    param_names[missing_names] <- paste0("V", which(missing_names))
  }

  target_idx <- seq_len(n_vars)
  if (!is.null(pars)) {
    if (is.numeric(pars)) {
      if (anyNA(pars) || any(!is.finite(pars)) ||
          any(pars != as.integer(pars))) {
        stop("'pars' must contain integer indices.", call. = FALSE)
      }
      target_idx <- tryCatch(
        seq_len(n_vars)[as.integer(pars)],
        error = function(e) NULL
      )
      if (is.null(target_idx) || length(target_idx) == 0L ||
          anyNA(target_idx)) {
        stop("The index specified in 'pars' was not found or is invalid.",
             call. = FALSE)
      }
    } else if (is.character(pars)) {
      base_names <- sub("\\[.*\\]$", "", param_names)
      if (any(startsWith(pars, "-"))) {
        exclude_names <- sub("^-", "", pars)
        target_idx <- which(
          !(param_names %in% exclude_names | base_names %in% exclude_names)
        )
      } else {
        target_idx <- unique(unlist(lapply(pars, function(par) {
          which(param_names == par | base_names == par)
        }), use.names = FALSE))
      }
      if (length(target_idx) == 0L) {
        stop("The variable name specified in 'pars' was not found.",
             call. = FALSE)
      }
    } else {
      stop("'pars' must be either numeric or character.", call. = FALSE)
    }
  } else {
    lp_idx <- which(param_names == "lp")
    target_idx <- c(lp_idx, setdiff(target_idx, lp_idx))
  }

  if (!is.null(max_rows)) {
    target_idx <- head(target_idx, max_rows)
  }

  result <- vector("list", length(target_idx))
  n_iter <- dim(draws)[1]
  n_selected_chains <- dim(draws)[2]

  for (i in seq_along(target_idx)) {
    p <- target_idx[i]
    mat_p <- matrix(
      draws[, , p],
      nrow = n_iter,
      ncol = n_selected_chains
    )
    valid_vec <- as.vector(mat_p)
    valid_vec <- valid_vec[is.finite(valid_vec)]

    if (length(valid_vec) == 0L) {
      result[[i]] <- data.frame(
        variable = param_names[p],
        mean = NA_real_,
        sd = NA_real_,
        map = NA_real_,
        q2.5 = NA_real_,
        q97.5 = NA_real_,
        ess_bulk = NA_real_,
        ess_tail = NA_real_,
        rhat = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    sd_val <- stats::sd(valid_vec)
    if (is.na(sd_val) || sd_val < 1e-10) {
      map_val <- valid_vec[1]
      q95 <- rep(valid_vec[1], 2)
      rhat_val <- NA_real_
      ebulk_val <- NA_real_
      etail_val <- NA_real_
    } else {
      map_val <- map_est(valid_vec)
      q95 <- quantile95(valid_vec)
      rhat_val <- r_hat(mat_p)
      ebulk_val <- ess_bulk(mat_p)
      etail_val <- ess_tail95(mat_p)
    }

    result[[i]] <- data.frame(
      variable = param_names[p],
      mean = mean(valid_vec),
      sd = sd_val,
      map = map_val,
      q2.5 = unname(q95[1]),
      q97.5 = unname(q95[2]),
      ess_bulk = ebulk_val,
      ess_tail = etail_val,
      rhat = rhat_val,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, result)
  numeric_cols <- vapply(result, is.numeric, logical(1))
  result[numeric_cols] <- lapply(result[numeric_cols], function(x) {
    x[abs(x) < 1e-12 & !is.na(x)] <- 0
    x
  })

  class(result) <- c("summary_BayesRTMB", "data.frame")
  attr(result, "digits") <- digits
  result
}
