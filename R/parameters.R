#' Parameter Types and Constraints in RTMB Models
#'
#' @description
#' When declaring parameters in the \code{parameters} block using \code{Dim()},
#' you can specify a \code{type} to automatically apply structural constraints.
#' These constraints ensure that parameters remain within their valid mathematical
#' space during estimation.
#'
#' @details
#' \strong{Standard Types:}
#' \itemize{
#'   \item \code{"vector"}, \code{"matrix"}, \code{"array"}: Standard unconstrained containers.
#'   Use \code{lower} and \code{upper} for element-wise bounds.
#' }
#'
#' \strong{Ordering Constraints:}
#' \itemize{
#'   \item \code{"ordered"}: A vector where $x_1 < x_2 < \dots < x_K$.
#'   \item \code{"positive_ordered"}: A vector where $0 < x_1 < x_2 < \dots < x_K$.
#' }
#'
#' \strong{Constrained Vectors:}
#' \itemize{
#'   \item \code{"simplex"}: A vector where all elements are >= 0 and sum(x) = 1.
#'   \item \code{"centered"}: A vector of length K where sum(x) = 0. (Estimated using K-1 degrees of freedom).
#' }
#'
#' \strong{Correlation and Covariance Matrices:}
#' \itemize{
#'   \item \code{"corr_matrix"}: A symmetric, positive-definite correlation matrix (diagonal elements are 1).
#'   \item \code{"cov_matrix"}: A symmetric, positive-definite covariance matrix.
#'   \item \code{"CF_corr"}: The Cholesky factor of a correlation matrix.
#'   \item \code{"CF_cov"}: The Cholesky factor of a covariance matrix (diagonal elements are positive).
#' }
#'
#' \strong{Special Structural Matrices:}
#' \itemize{
#'   \item \code{"lower_tri"}: A lower-triangular matrix.
#'   \item \code{"positive_lower_tri"}: A lower-triangular matrix with positive diagonal elements.
#'   \item \code{"centered_matrix"}: A matrix where each column sums to zero.
#'   \item \code{"centered_tri"}: A triangular matrix with column-wise sum-to-zero constraints (often used for identification in factor analysis).
#' }
#'
#' @name parameter_types
#' @seealso \code{\link{rtmb_code}}, \code{\link{math_functions}}, \code{\link{distributions}}
NULL

#' Define parameter dimensions and types
#'
#' @param dim Dimensions of the parameter.
#' @param type Type of the parameter.
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @param random Logical; whether it is a random effect.
#'
#' @export
Dim <- function(dim = 1, type = NULL, lower = NULL, upper = NULL, random = FALSE) {
  # 1. Basic Type Handling
  if (is.null(type)) {
    if (length(dim) == 1) type <- "vector"
    else if (length(dim) == 2) type <- "matrix"
    else type <- "array"
  }

  bounds_type <- "none"
  if (type %in% c("ordered", "positive_ordered", "simplex", "corr_matrix",
                  "cov_matrix", "CF_corr", "CF_cov", "centered",
                  "centered_matrix", "lower_tri", "lower_tri_stz",
                  "centered_tri", "positive_centered_tri",
                  "positive_lower_tri")) {
    bounds_type <- type
  } else if (!is.null(lower) && is.null(upper)) bounds_type <- "lower"
  else if (is.null(lower) && !is.null(upper)) bounds_type <- "upper"
  else if (!is.null(lower) && !is.null(upper)) bounds_type <- "interval"

  # 2. Multi-dimensional recursion (K x P x P convention)
  if (length(dim) > 2 && bounds_type %in% c("corr_matrix", "CF_corr", "cov_matrix", "CF_cov",
                                          "lower_tri", "lower_tri_stz", "positive_lower_tri",
                                          "centered_matrix", "centered_tri", "positive_centered_tri")) {
    K <- dim[1]
    sub_dim <- dim[-1]
    sub_p <- Dim(dim = sub_dim, type = type, lower = lower, upper = upper, random = random)

    return(list(
      dim        = dim,
      length     = prod(dim),
      unc_length = sub_p$unc_length * K,
      type       = type,
      bounds     = bounds_type,
      lower      = lower,
      upper      = upper,
      random     = random
    ))
  }

  # 3. Standard 1D/2D Calculation
  calc_length <- prod(dim)
  calc_unc_length <- calc_length

  if (bounds_type %in% c("corr_matrix", "CF_corr")) {
    R <- dim[1]
    C <- if (length(dim) > 1) dim[2] else dim[1]
    if (bounds_type == "corr_matrix") {
      if (R != C) stop(sprintf("[Parameter definition error] Variable specified as type = 'corr_matrix', but it does not have square matrix dimensions (Specified dimension: %s).", paste(dim, collapse = " x ")))
      calc_unc_length <- R * (R - 1) / 2
    } else {
      calc_unc_length <- 0
      if (R >= 2) {
        for (i in 2:R) {
          calc_unc_length <- calc_unc_length + min(i - 1, C - 1)
        }
      }
    }
  } else if (bounds_type %in% c("cov_matrix", "CF_cov")) {
    R <- dim[1]
    C <- if (length(dim) > 1) dim[2] else dim[1]
    if (bounds_type == "cov_matrix") {
      if (R != C) stop(sprintf("[Parameter definition error] Variable specified as type = 'cov_matrix', but it does not have square matrix dimensions (Specified dimension: %s).", paste(dim, collapse = " x ")))
      calc_unc_length <- R * (R + 1) / 2
    } else {
      calc_unc_length <- 0
      for (i in 1:R) {
        calc_unc_length <- calc_unc_length + min(i, C)
      }
    }
  } else if (bounds_type %in% c("simplex", "centered")) {
    calc_unc_length <- calc_length - 1
  } else if (bounds_type == "centered_matrix") {
    R <- dim[1]
    C <- dim[2]
    calc_unc_length <- (R - 1) * C
  } else if (bounds_type %in% c("centered_tri", "positive_centered_tri")) {
    R <- dim[1]
    C <- dim[2]
    calc_unc_length <- 0
    for (d in 1:C) {
      if (d <= R - 1) calc_unc_length <- calc_unc_length + (R - d)
    }
  } else if (bounds_type %in% c("lower_tri", "lower_tri_stz", "positive_lower_tri")) {
    R <- dim[1]
    C <- if (length(dim) > 1) dim[2] else dim[1]
    if (R >= C) {
      calc_unc_length <- C * (C + 1) / 2 + (R - C) * C
    } else {
      calc_unc_length <- R * (R + 1) / 2
    }
    if (bounds_type == "lower_tri_stz") calc_unc_length <- calc_unc_length - 1
  }

  list(
    dim        = dim,
    length     = calc_length,
    unc_length = calc_unc_length,
    type       = type,
    bounds     = bounds_type,
    lower      = lower,
    upper      = upper,
    random     = random
  )
}


generate_flat_names <- function(base_name, dims, names_def = NULL) {
  len <- prod(dims)
  if (is.null(dims)) dims <- len

  if (len == 1) {
    if (!is.null(names_def) && length(names_def) == 1) {
      return(paste0(base_name, "[", names_def, "]"))
    }
    return(base_name)
  }

  if (length(dims) > 1) {
    if (!is.null(names_def)) {
      if (is.list(names_def) && length(names_def) == length(dims)) {
        # When names are specified as a list per dimension
        grid <- do.call(expand.grid, c(names_def, list(stringsAsFactors = FALSE)))
        return(paste0(base_name, "[", apply(grid, 1, paste, collapse = ","), "]"))
      } else if (is.atomic(names_def) && length(names_def) == dims[1]) {
        # When names are specified as a vector
        if (length(dims) == 2 && dims[1] == dims[2]) {
          # Square matrix (e.g., CF_corr): Apply to both rows and columns
          grid <- expand.grid(names_def, names_def, stringsAsFactors = FALSE)
        } else {
          # Non-square matrix: Apply only to rows
          dim_list <- lapply(dims, seq_len)
          dim_list[[1]] <- names_def
          grid <- do.call(expand.grid, dim_list)
        }
        return(paste0(base_name, "[", apply(grid, 1, paste, collapse = ","), "]"))
      }
    }
    # When names are not specified, or sizes do not match
    grid <- do.call(expand.grid, lapply(dims, seq_len))
    return(paste0(base_name, "[", apply(grid, 1, paste, collapse = ","), "]"))
  } else {
    # For 1D vectors
    if (!is.null(names_def) && length(names_def) == len) {
      return(paste0(base_name, "[", names_def, "]"))
    }
    return(paste0(base_name, "[", seq_len(len), "]"))
  }
}



# Function to restore a flat vector in constrained space to a list
constrained_vector_to_list <- function(vec, par_list) {
  res <- list()
  idx <- 1
  for (name in names(par_list)) {
    p <- par_list[[name]]
    len <- p$length
    if (len > 0) {
      val <- vec[idx:(idx + len - 1)]
      if (length(p$dim) > 1 && prod(p$dim) == len) {
        dim(val) <- p$dim
      }
      res[[name]] <- val
      idx <- idx + len
    } else {
      res[[name]] <- numeric(0)
    }
  }
  return(res)
}

# calc_log_jacobian
calc_log_jacobian <- function(para_unc_list, par_list, only_random = FALSE) {
  lj <- 0
  for (name in names(par_list)) {
    p <- par_list[[name]]
    val_unc <- para_unc_list[[name]]
    b_type <- p$bounds

    if (only_random && !isTRUE(p$random)) next

    if (b_type == "lower" || b_type == "upper") {
      lj <- lj + sum(val_unc)
    } else if (b_type == "interval") {
      lj <- lj + sum(log(p$upper - p$lower) - val_unc - 2 * log(1 + exp(-val_unc)))
    } else if (b_type == "ordered") {
      if (p$length > 1) {
        if (is.matrix(val_unc)) {
          lj <- lj + sum(val_unc[, -1])
        } else {
          lj <- lj + sum(val_unc[2:p$length])
        }
      }
    } else if (b_type == "positive_ordered") {
      lj <- lj + sum(val_unc)
    } else if (b_type == "simplex") {
      K <- p$length
      z <- 1 / (1 + exp(-(val_unc - log(1 / (K - seq_len(K - 1))))))
      lj <- lj + sum(log(z) + log(1 - z))
      if (K > 2) {
        for (k in 1:(K - 2)) {
          lj <- lj + (K - k - 1) * log(1 - z[k])
        }
      }
    } else if (b_type %in% c("centered", "centered_matrix", "centered_tri")) {
      diff_dim <- p$length - p$unc_length
      lj <- lj + (diff_dim / 2) * log(2 * pi)
    } else if (b_type == "positive_centered_tri") {
      diff_dim <- p$length - p$unc_length
      lj <- lj + (diff_dim / 2) * log(2 * pi)
      R <- p$dim[1]
      C <- p$dim[2]
      idx <- 1
      for (d in 1:min(C, R - 1)) {
        n_d <- R - d + 1
        u1 <- val_unc[idx]
        lj <- lj + u1 + 0.5 * log(n_d / (n_d - 1))
        idx <- idx + n_d - 1
      }
    } else if (b_type %in% c("corr_matrix", "CF_corr")) {
      K_rep <- if (length(p$dim) > 2) p$dim[1] else 1
      P_mat <- if (length(p$dim) > 2) p$dim[2] else p$dim[1]

      unc_per_slice <- length(val_unc) / K_rep
      for (k in 1:K_rep) {
        v_u <- val_unc[((k-1)*unc_per_slice + 1):(k*unc_per_slice)]
        lj_term <- 0
        idx_u <- 1
        for (i in 2:P_mat) {
          for (j in 1:(i - 1)) {
            # CF_corr と corr_matrix で正しいヤコビアンを使い分ける修正
            if (b_type == "corr_matrix") {
              lj_term <- lj_term - 2 * (P_mat - i + 1) * log(cosh(v_u[idx_u]))
            } else { # CF_corr
              lj_term <- lj_term - (i - j + 1) * log(cosh(v_u[idx_u]))
            }
            idx_u <- idx_u + 1
          }
        }
        lj <- lj + lj_term
      }
    } else if (b_type %in% c("cov_matrix", "CF_cov")) {
      K_rep <- if (length(p$dim) > 2) p$dim[1] else 1
      P_mat <- if (length(p$dim) > 2) p$dim[2] else p$dim[1]

      unc_per_slice <- length(val_unc) / K_rep
      for (k in 1:K_rep) {
        v_u <- val_unc[((k-1)*unc_per_slice + 1):(k*unc_per_slice)]
        idx_u <- 1
        for (i in 1:P_mat) {
          for (j in 1:i) {
            if (i == j) {
              if (b_type == "cov_matrix") lj <- lj + (P_mat - i + 2) * v_u[idx_u]
              else lj <- lj + v_u[idx_u]
            }
            idx_u <- idx_u + 1
          }
        }
      }
    } else if (b_type %in% c("lower_tri", "positive_lower_tri")) {
      K_rep <- if (length(p$dim) > 2) p$dim[1] else 1
      R_mat <- if (length(p$dim) > 2) p$dim[2] else p$dim[1]
      C_mat <- if (length(p$dim) > 2) p$dim[3] else (if (length(p$dim) > 1) p$dim[2] else p$dim[1])

      unc_per_slice <- length(val_unc) / K_rep
      for (k in 1:K_rep) {
        v_u <- val_unc[((k-1)*unc_per_slice + 1):(k*unc_per_slice)]
        idx_u <- 1
        for (i in 1:R_mat) {
          max_j <- min(i, C_mat)
          for (j in 1:max_j) {
            if (i == j && b_type == "positive_lower_tri") {
              lj <- lj + v_u[idx_u]
            }
            idx_u <- idx_u + 1
          }
        }
      }
    }
  }
  return(lj)
}
#' Generate Random Initial Values
#'
#' @description Generates random initial values for model parameters by drawing from a uniform distribution
#' on the unconstrained scale and then transforming them to the constrained scale.
#'
#' @param pl_full A list containing the full parameter structure, including the `init` field which serves as a template.
#' @param par_list A list defining the structure and constraints of the parameters, including `unc_length`.
#' @param range Numeric; the range for the uniform distribution `[-range, range]` used for generating unconstrained values. Default is 2.
#' @return A numeric vector of initial values where `NA` elements are replaced with randomly generated constrained values.
generate_random_init <- function(pl_full, par_list, range = 2) {
  unc_list <- list()
  for(name in names(par_list)) {
    unc_list[[name]] <- runif(par_list[[name]]$unc_length, min = -range, max = range)
  }

  con_list <- to_constrained(unc_list, par_list)

  # Group-wise NA checking instead of element-wise NA checking to preserve structures
  init_list <- constrained_vector_to_list(pl_full$init, par_list)
  for (name in names(par_list)) {
    # If ANY element in the parameter block is NA, we replace the ENTIRE block with random values
    if (any(is.na(init_list[[name]]))) {
      init_list[[name]] <- con_list[[name]]
    }
  }

  return(unlist(init_list, use.names = FALSE))
}

parse_parameters <- function(par_list, par_names = NULL) {
  P <- if (length(par_list) > 0) sum(sapply(par_list, function(x) x$length)) else 0

  init_vec    <- numeric(P)
  lower_vec   <- rep(NA, P)
  upper_vec   <- rep(NA, P)
  bounds_type <- character(P)
  flat_names  <- character(P)

  idx <- 1
  for (name in names(par_list)) {
    p <- par_list[[name]]
    len <- p$length
    end_idx <- idx + len - 1
    names_def <- if (!is.null(par_names)) par_names[[name]] else NULL
    flat_names[idx:end_idx] <- generate_flat_names(name, p$dim, names_def)

    if (!is.null(p$lower)) lower_vec[idx:end_idx] <- p$lower
    if (!is.null(p$upper)) upper_vec[idx:end_idx] <- p$upper

    b_type <- "none"
    if (!is.null(p$lower) && is.null(p$upper)) b_type <- "lower"
    if (is.null(p$lower) && !is.null(p$upper)) b_type <- "upper"
    if (!is.null(p$lower) && !is.null(p$upper)) b_type <- "interval"
    bounds_type[idx:end_idx] <- b_type

    if (!is.null(p$init)) {
      init_vec[idx:end_idx] <- p$init
    } else {
      init_vec[idx:end_idx] <- NA
    }

    idx <- idx + len
  }

  return(list(
    names  = flat_names,
    init   = init_vec,
    lower  = lower_vec,
    upper  = upper_vec,
    bounds = bounds_type,
    idx    = list(
      lower = which(bounds_type == "lower"),
      upper = which(bounds_type == "upper"),
      intv  = which(bounds_type == "interval")
    )
  ))
}

unpack_parameters <- function(para, par_list) {
  res <- list()
  idx <- 1
  for (name in names(par_list)) {
    p <- par_list[[name]]
    len <- p$length
    val <- para[idx:(idx + len - 1)]
    if (length(p$dim) > 1) {
      dim(val) <- p$dim
    }
    res[[name]] <- val
    idx <- idx + len
  }
  return(res)
}

#' Select parameters from a list by name or index
#' @param all_list A named list of parameters.
#' @param pars A character or numeric vector.
#' @return A subset of the list.
#' @keywords internal
select_parameters <- function(all_list, pars) {
  if (is.null(pars)) return(all_list)
  if (identical(pars, "parameters") || identical(pars, "all")) return(all_list)
  
  # 1. Numeric indexing (positive or negative)
  if (is.numeric(pars)) {
    # We use a trick to handle both positive and negative indices correctly in one go
    # However, base R already handles this: list[pars]
    # But we want to ensure we don't return out-of-bounds NAs if possible
    res <- tryCatch(all_list[pars], error = function(e) NULL)
    if (is.null(res)) return(all_list[0]) # Return empty list if error
    return(res[names(res) != ""]) # Filter out any NA names that might appear from out-of-bounds positive indices
  }
  
  # 2. Character indexing
  if (is.character(pars)) {
    is_negative <- grepl("^-", pars)
    if (any(is_negative)) {
      if (!all(is_negative)) {
        stop("Cannot mix positive and negative parameter names in selection.")
      }
      exclude_names <- gsub("^-", "", pars)
      return(all_list[!(names(all_list) %in% exclude_names)])
    } else {
      # Standard positive selection
      return(all_list[names(all_list) %in% pars])
    }
  }
  
  return(all_list)
}
