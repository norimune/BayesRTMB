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
#'   \item \code{"simplex"}: A vector where all elements are $\ge 0$ and $\sum x = 1$.
#'   \item \code{"sum_to_zero"}: A vector of length $K$ where $\sum x = 0$. (Estimated using $K-1$ degrees of freedom).
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
  if (is.null(type)) {
    if (length(dim) == 1) type <- "vector"
    else if (length(dim) == 2) type <- "matrix"
    else type <- "array"
  }

  bounds_type <- "none"
  if (type %in% c("ordered", "positive_ordered", "simplex", "corr_matrix",
                  "cov_matrix", "CF_corr", "CF_cov", "sum_to_zero",
                  "centered_matrix", "lower_tri", "lower_tri_stz",
                  "centered_tri", "positive_centered_tri",
                  "positive_lower_tri")) {
    bounds_type <- type
  } else if (!is.null(lower) && is.null(upper)) bounds_type <- "lower"
  else if (is.null(lower) && !is.null(upper)) bounds_type <- "upper"
  else if (!is.null(lower) && !is.null(upper)) bounds_type <- "interval"

  calc_length <- prod(dim)
  calc_unc_length <- calc_length

  if (bounds_type %in% c("corr_matrix", "CF_corr")) {
    R <- dim[1]
    C <- if (length(dim) > 1) dim[2] else dim[1]
    if (bounds_type == "corr_matrix") {
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
      calc_unc_length <- R * (R + 1) / 2
    } else {
      calc_unc_length <- 0
      for (i in 1:R) {
        calc_unc_length <- calc_unc_length + min(i, C)
      }
    }
  } else if (bounds_type %in% c("simplex", "sum_to_zero")) {
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
    lower      = lower,
    upper      = upper,
    random     = random,
    bounds     = bounds_type
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

stz_basis <- function(K) {
  Q <- matrix(0, nrow = K, ncol = K - 1)
  for (j in 1:(K - 1)) {
    Q[1:j, j] <- 1 / sqrt(j * (j + 1))
    Q[j + 1, j] <- -j / sqrt(j * (j + 1))
  }
  return(Q)
}

# Function to restore a flat vector in constrained space to a list
constrained_vector_to_list <- function(vec, par_list) {
  res <- list()
  idx <- 1
  for (name in names(par_list)) {
    p <- par_list[[name]]
    len <- p$length
    val <- vec[idx:(idx + len - 1)]

    if (length(p$dim) > 1 && !(p$bounds %in% c("corr_matrix", "CF_corr", "cov_matrix", "CF_cov", "lower_tri", "lower_tri_stz"))) {
      dim(val) <- p$dim
    } else if (p$bounds %in% c("corr_matrix", "cov_matrix")) {
      dim(val) <- c(p$dim[1], p$dim[1])
    } else if (p$bounds %in% c("CF_corr", "CF_cov", "lower_tri", "lower_tri_stz", "positive_lower_tri")) {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      dim(val) <- c(R, C)
    }

    res[[name]] <- val
    idx <- idx + len
  }
  return(res)
}

# Function to restore a flat vector in unconstrained space to a list (*Added due to omission)
unconstrained_vector_to_list <- function(vec, par_list) {
  res <- list()
  idx <- 1
  for (name in names(par_list)) {
    p <- par_list[[name]]
    len <- p$unc_length
    if (len > 0) {
      val <- vec[idx:(idx + len - 1)]
      res[[name]] <- val
      idx <- idx + len
    } else {
      res[[name]] <- numeric(0)
    }
  }
  return(res)
}

# Convert to unconstrained space
to_unconstrained <- function(para_orig_list, par_list) {
  para_unc <- list()
  for (name in names(par_list)) {
    p <- par_list[[name]]
    val_orig <- para_orig_list[[name]]
    b_type <- p$bounds

    if (b_type == "lower") {
      para_unc[[name]] <- log(val_orig - p$lower)
    } else if (b_type == "upper") {
      para_unc[[name]] <- log(p$upper - val_orig)
    } else if (b_type == "interval") {
      prob <- (val_orig - p$lower) / (p$upper - p$lower)
      para_unc[[name]] <- log(prob / (1 - prob))
    } else if (b_type == "ordered") {
      y <- numeric(p$length)
      y[1] <- val_orig[1]
      if (p$length > 1) {
        for (i in 2:p$length) {
          y[i] <- log(val_orig[i] - val_orig[i-1])
        }
      }
      para_unc[[name]] <- y
    } else if (b_type == "positive_ordered") {
      y <- numeric(p$length)
      y[1] <- log(val_orig[1])
      if (p$length > 1) {
        for (i in 2:p$length) {
          y[i] <- log(val_orig[i] - val_orig[i-1])
        }
      }
      para_unc[[name]] <- y
    } else if (b_type == "simplex") {
      K <- p$length
      y <- numeric(K - 1)
      stick_len <- 1
      for (k in 1:(K - 1)) {
        z_k <- if (stick_len > 1e-10) val_orig[k] / stick_len else 0
        z_k <- min(max(z_k, 1e-10), 1 - 1e-10)
        y[k] <- log(z_k / (1 - z_k)) - log(1 / (K - k))
        stick_len <- stick_len - val_orig[k]
      }
      para_unc[[name]] <- y
    } else if (b_type == "sum_to_zero") {
      K <- p$length
      Q <- matrix(0, nrow = K, ncol = K - 1)
      for (j in 1:(K - 1)) {
        Q[1:j, j] <- 1 / sqrt(j * (j + 1))
        Q[j + 1, j] <- -j / sqrt(j * (j + 1))
      }
      para_unc[[name]] <- as.numeric(crossprod(Q, val_orig))

    } else if (b_type == "centered_matrix") {
      R <- p$dim[1]
      C <- p$dim[2]
      Q <- matrix(0, nrow = R, ncol = R - 1)
      for (j in 1:(R - 1)) {
        Q[1:j, j] <- 1 / sqrt(j * (j + 1))
        Q[j + 1, j] <- -j / sqrt(j * (j + 1))
      }
      Y <- crossprod(Q, val_orig)
      para_unc[[name]] <- as.numeric(Y)

    } else if (b_type == "centered_tri") {
      R <- p$dim[1]
      C <- p$dim[2]
      y <- numeric(p$unc_length)
      idx <- 1
      for (d in 1:min(C, R - 1)) {
        x_d <- val_orig[d:R, d]
        n_d <- length(x_d)
        Q_d <- stz_basis(n_d)
        z_d <- as.numeric(crossprod(Q_d, x_d))
        y[idx:(idx + n_d - 2)] <- z_d
        idx <- idx + n_d - 1
      }
      para_unc[[name]] <- y

    } else if (b_type == "positive_centered_tri") {
      R <- p$dim[1]
      C <- p$dim[2]
      y <- numeric(p$unc_length)
      idx <- 1
      for (d in 1:min(C, R - 1)) {
        x_d <- val_orig[d:R, d]
        n_d <- length(x_d)

        u1 <- log(x_d[1])
        y[idx] <- u1

        if (n_d > 2) {
          v <- x_d[2:n_d] + x_d[1] / (n_d - 1)
          Q_d <- stz_basis(n_d - 1)
          z_d <- as.numeric(crossprod(Q_d, v))
          y[(idx + 1):(idx + n_d - 2)] <- z_d
        }
        idx <- idx + n_d - 1
      }
      para_unc[[name]] <- y

    } else if (b_type %in% c("corr_matrix", "CF_corr")) {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      L <- if (b_type == "corr_matrix") t(chol(val_orig + diag(1e-8, R))) else val_orig
      z <- numeric(p$unc_length)
      idx <- 1
      if (R >= 2) {
        for (i in 2:R) {
          prod_term <- 1
          max_j <- min(i - 1, C - 1)
          if (max_j > 0) {
            for (j in 1:max_j) {
              z_val <- L[i, j] / sqrt(prod_term)
              z_val <- min(max(z_val, -0.999999), 0.999999)
              z[idx] <- z_val
              prod_term <- prod_term * (1 - z_val^2)
              idx <- idx + 1
            }
          }
        }
      }
      para_unc[[name]] <- atanh(z)
    } else if (b_type %in% c("cov_matrix", "CF_cov")) {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      L <- if (b_type == "cov_matrix") t(chol(val_orig + diag(1e-8, R))) else val_orig
      y <- numeric(p$unc_length)
      idx <- 1
      for (i in 1:R) {
        max_j <- min(i, C)
        for (j in 1:max_j) {
          if (i == j) {
            y[idx] <- log(L[i, j])
          } else {
            y[idx] <- L[i, j]
          }
          idx <- idx + 1
        }
      }
      para_unc[[name]] <- y
    } else if (b_type %in% c("lower_tri", "positive_lower_tri")) {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      y <- numeric(p$unc_length)
      idx <- 1
      for (i in 1:R) {
        max_j <- min(i, C)
        for (j in 1:max_j) {
          y[idx] <- val_orig[i, j]
          idx <- idx + 1
        }
      }
      para_unc[[name]] <- y
    } else if (b_type == "lower_tri_stz") {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      y <- numeric(p$unc_length)
      idx <- 1
      last_i <- R
      last_j <- min(R, C)
      for (i in 1:R) {
        max_j <- min(i, C)
        for (j in 1:max_j) {
          if (i == last_i && j == last_j) break
          y[idx] <- val_orig[i, j]
          idx <- idx + 1
        }
      }
      para_unc[[name]] <- y
    } else {
      para_unc[[name]] <- as.numeric(val_orig)
    }
  }
  return(para_unc)
}

to_constrained <- function(para_unc_list, par_list) {
  para <- para_unc_list
  for (name in names(par_list)) {
    p <- par_list[[name]]
    val_unc <- para_unc_list[[name]]
    b_type <- p$bounds

    # Zero (scalar) to preserve RTMB's advector attribute
    ad_zero <- val_unc[1] * 0

    if (b_type == "lower") {
      para[[name]] <- p$lower + exp(val_unc)
    } else if (b_type == "upper") {
      para[[name]] <- p$upper - exp(val_unc)
    } else if (b_type == "interval") {
      prob <- plogis(val_unc)
      para[[name]] <- p$lower + (p$upper - p$lower) * prob
    } else if (b_type == "ordered") {
      if (p$length > 1) {
        para[[name]] <- val_unc[1] + cumsum(c(ad_zero, exp(val_unc[-1])))
      } else {
        para[[name]] <- val_unc
      }
    } else if (b_type == "positive_ordered") {
      para[[name]] <- cumsum(exp(val_unc))
    } else if (b_type == "simplex") {
      K <- p$length
      z <- 1 / (1 + exp(-(val_unc - log(1 / (K - seq_len(K - 1))))))
      ad_one <- ad_zero + 1
      rem_prev <- c(ad_one, cumprod(1 - z)[- (K - 1)])

      x_first <- rem_prev * z
      x_last <- cumprod(1 - z)[K - 1]
      para[[name]] <- c(x_first, x_last)

    } else if (b_type == "sum_to_zero") {
      K <- p$length
      Q <- matrix(0, nrow = K, ncol = K - 1)
      for (j in 1:(K - 1)) {
        Q[1:j, j] <- 1 / sqrt(j * (j + 1))
        Q[j + 1, j] <- -j / sqrt(j * (j + 1))
      }
      para[[name]] <- as.numeric(Q %*% val_unc)

    } else if (b_type == "centered_matrix") {
      R <- p$dim[1]
      C <- p$dim[2]
      Q <- matrix(0, nrow = R, ncol = R - 1)
      for (j in 1:(R - 1)) {
        Q[1:j, j] <- 1 / sqrt(j * (j + 1))
        Q[j + 1, j] <- -j / sqrt(j * (j + 1))
      }
      mat_unc <- matrix(val_unc, nrow = R - 1, ncol = C)
      para[[name]] <- Q %*% mat_unc

    } else if (b_type == "centered_tri") {
      R <- p$dim[1]
      C <- p$dim[2]
      L <- matrix(ad_zero, nrow = R, ncol = C)
      idx <- 1
      for (d in 1:min(C, R - 1)) {
        n_d <- R - d + 1
        z_d <- val_unc[idx:(idx + n_d - 2)]
        Q_d <- stz_basis(n_d)
        x_d <- as.numeric(Q_d %*% z_d)
        L[d:R, d] <- x_d
        idx <- idx + n_d - 1
      }
      para[[name]] <- L

    } else if (b_type == "positive_centered_tri") {
      R <- p$dim[1]
      C <- p$dim[2]
      L <- matrix(ad_zero, nrow = R, ncol = C)
      idx <- 1
      for (d in 1:min(C, R - 1)) {
        n_d <- R - d + 1
        x_d <- rep(ad_zero, n_d)

        u1 <- val_unc[idx]
        x1 <- exp(u1)
        x_d[1] <- x1

        if (n_d > 2) {
          z_d <- val_unc[(idx + 1):(idx + n_d - 2)]
          Q_d <- stz_basis(n_d - 1)
          v <- as.numeric(Q_d %*% z_d)
          x_d[2:n_d] <- v - x1 / (n_d - 1)
        } else if (n_d == 2) {
          x_d[2] <- -x1
        }
        L[d:R, d] <- x_d
        idx <- idx + n_d - 1
      }
      para[[name]] <- L

    } else if (b_type %in% c("corr_matrix", "CF_corr")) {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      z <- tanh(val_unc)
      L <- matrix(ad_zero, nrow = R, ncol = C)
      L[1, 1] <- 1
      idx <- 1
      if (R >= 2) {
        for (i in 2:R) {
          prod_term <- 1
          max_j <- min(i - 1, C - 1)
          if (max_j > 0) {
            for (j in 1:max_j) {
              L[i, j] <- z[idx] * sqrt(prod_term)
              prod_term <- prod_term * (1 - z[idx]^2)
              idx <- idx + 1
            }
          }
          last_j <- min(i, C)
          if (last_j > max_j) {
            L[i, last_j] <- sqrt(abs(prod_term))
          }
        }
      }
      if (b_type == "corr_matrix") para[[name]] <- L %*% t(L)
      else para[[name]] <- L

    } else if (b_type %in% c("cov_matrix", "CF_cov")) {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      L <- matrix(ad_zero, nrow = R, ncol = C)
      idx <- 1
      for (i in 1:R) {
        max_j <- min(i, C)
        for (j in 1:max_j) {
          if (i == j) L[i, j] <- exp(val_unc[idx])
          else L[i, j] <- val_unc[idx]
          idx <- idx + 1
        }
      }
      if (b_type == "cov_matrix") para[[name]] <- L %*% t(L)
      else para[[name]] <- L

    } else if (b_type %in% c("lower_tri", "positive_lower_tri")) {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      L <- matrix(ad_zero, nrow = R, ncol = C)
      idx <- 1
      for (i in 1:R) {
        max_j <- min(i, C)
        for (j in 1:max_j) {
          L[i, j] <- val_unc[idx]
          idx <- idx + 1
        }
      }
      para[[name]] <- L

    } else if (b_type == "lower_tri_stz") {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      L <- matrix(ad_zero, nrow = R, ncol = C)
      idx <- 1
      sum_val <- 0
      last_i <- R
      last_j <- min(R, C)
      for (i in 1:R) {
        max_j <- min(i, C)
        for (j in 1:max_j) {
          if (i == last_i && j == last_j) break
          L[i, j] <- val_unc[idx]
          sum_val <- sum_val + val_unc[idx]
          idx <- idx + 1
        }
      }
      L[last_i, last_j] <- -sum_val
      para[[name]] <- L
    }

    if (!(b_type %in% c("corr_matrix", "CF_corr", "cov_matrix", "CF_cov", "lower_tri", "lower_tri_stz", "positive_lower_tri")) && length(p$dim) > 1) {
      dim(para[[name]]) <- p$dim
    }
  }
  return(para)
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
        lj <- lj + sum(val_unc[2:p$length])
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
    } else if (b_type %in% c("sum_to_zero", "centered_matrix", "centered_tri")) {
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
        # Jacobian log(|J|) = u1 + 0.5 * log(n_d / (n_d - 1))
        lj <- lj + u1 + 0.5 * log(n_d / (n_d - 1))
        idx <- idx + n_d - 1
      }
    } else if (b_type %in% c("corr_matrix", "CF_corr")) {
      K <- p$dim[1]
      lj_term <- 0
      idx <- 1
      for (i in 2:K) {
        for (j in 1:(i - 1)) {
          lj_term <- lj_term - 2 * (K - i + 1) * log(cosh(val_unc[idx]))
          idx <- idx + 1
        }
      }
      lj <- lj + lj_term
    } else if (b_type == "cov_matrix") {
      K <- p$dim[1]
      idx <- 1
      for (i in 1:K) {
        for (j in 1:i) {
          if (i == j) lj <- lj + (K - i + 2) * val_unc[idx]
          idx <- idx + 1
        }
      }
    } else if (b_type == "CF_cov") {
      K <- p$dim[1]
      idx <- 1
      for (i in 1:K) {
        for (j in 1:i) {
          if (i == j) lj <- lj + val_unc[idx]
          idx <- idx + 1
        }
      }
    } else if (b_type == "positive_lower_tri") {
      R <- p$dim[1]
      C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]
      idx <- 1
      for (i in 1:R) {
        max_j <- min(i, C)
        for (j in 1:max_j) {
          if (i == j) lj <- lj + val_unc[idx] # Add only diagonal elements
          idx <- idx + 1
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
  random_pool <- unlist(con_list, use.names = FALSE)

  init_final <- pl_full$init
  idx_na <- is.na(init_final)
  if (any(idx_na)) {
    init_final[idx_na] <- random_pool[idx_na]
  }
  return(init_final)
}

parse_parameters <- function(par_list, par_names = NULL) {
  P <- sum(sapply(par_list, function(x) x$length))

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
