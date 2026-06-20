#' Mathematical and Matrix Utility Functions for RTMB Models
#'
#' @description
#' This page summarizes the utility functions provided to assist in writing
#' efficient and numerically stable model code within `rtmb_code`. These functions
#' are specifically designed to be compatible with RTMB's Automatic Differentiation (AD).
#'
#' @details
#' \strong{1. Link and Inverse Functions:}
#' \itemize{
#'   \item \code{logit(x)}: Computes the logit transformation \code{log(x/(1-x))}.
#'   \item \code{inv_logit(x)}: Computes the inverse logit (logistic) transformation \code{1/(1+exp(-x))}.
#' }
#'
#' \strong{2. Numerical Stability Utilities:}
#' These functions use specialized algorithms to prevent overflow/underflow in AD calculations.
#' \itemize{
#'   \item \code{log_sum_exp(x)}: Safely computes \code{log(sum(exp(x)))} using the log-sum-exp trick.
#'   \item \code{log1p_exp(x)}: Computes \code{log(1 + exp(x))} stably.
#'   \item \code{log1m_exp(x)}: Computes \code{log(1 - exp(x))} stably for \code{x < 0}.
#'   \item \code{log_softmax(x)}: Computes the log of the softmax function.
#'   \item \code{fabs(x)}: A smooth version of the absolute value function \code{abs(x)} to ensure differentiability at zero.
#' }
#'
#' \strong{3. Matrix and Vector Transformations:}
#' Used to convert unconstrained vectors into structured matrices (e.g., for factor analysis or identification constraints).
#' \itemize{
#'   \item \code{sum_to_zero(x)}: Transforms a vector of length K-1 into a vector of length K that sums to zero.
#'   \item \code{to_lower_tri(x, M, D)}: Fills a matrix of size M x D with elements from vector \code{x} in a lower-triangular fashion.
#'   \item \code{to_centered_matrix(x, R, C)}: Creates an R x C matrix where each column sums to zero.
#'   \item \code{to_centered_tri(x, R, C)}: Creates an R x C matrix with column-wise sum-to-zero constraints on the lower elements (useful for identification in factor analysis).
#'   \item \code{rtmb_vector(value, length)}: Creates an AD-compatible vector initialized with \code{value}.
#'   \item \code{rtmb_array(value, dim)}: Creates an AD-compatible array initialized with \code{value}.
#' }
#'
#' \strong{4. Linear Algebra for AD:}
#' \itemize{
#'   \item \code{log_det_chol(L)}: Calculates the log-determinant of a covariance matrix from its Cholesky factor \code{L}.
#'   \item \code{quad_form_chol(x, L)}: Computes the quadratic form x^T Sigma^(-1) x using the Cholesky factor \code{L}.
#'   \item \code{distance(x, y)}: Computes the Euclidean distance between two vectors with a small epsilon for stability.
#' }
#' @name math_functions
#' @family utilities
#' @import RTMB
NULL

#' Create an AD-compatible vector
#'
#' @description
#' Creates a fixed-length vector that is safe to use inside `rtmb_code()` when
#' subsequent assignments may involve RTMB automatic-differentiation values.
#' This is a safer alternative to `numeric(length)` in model code.
#'
#' @param value Initial value. Usually a scalar such as `0`.
#' @param length Length of the vector.
#' @param seed Optional AD value used to create the container. Usually this can
#'   be left as `NULL`; inside `rtmb_code()`, BayesRTMB automatically provides an
#'   AD seed from the model parameters when available.
#' @return An AD-compatible vector.
#' @export
rtmb_vector <- function(value = 0, length, seed = NULL) {
  if (missing(length)) {
    stop("'length' must be supplied.", call. = FALSE)
  }
  if (!is.numeric(length) || base::length(length) != 1L || is.na(length) || length < 0) {
    stop("'length' must be a non-negative scalar.", call. = FALSE)
  }
  length <- as.integer(length)
  if (length == 0L) {
    return(RTMB::advector(numeric(0)))
  }

  if (is.null(seed)) {
    seed <- .rtmb_auto_ad_seed(parent.frame())
  }

  if (!is.null(seed)) {
    if (base::length(seed) < 1L) {
      stop("'seed' must contain at least one value.", call. = FALSE)
    }
    if (base::length(value) != 1L && base::length(value) != length) {
      stop("'value' must be length 1 or have the requested vector length.", call. = FALSE)
    }
    return(seed[1] * rep(0, length) + value)
  }

  if (inherits(value, "advector")) {
    if (base::length(value) == 1L) return(rep(value, length))
    return(value)
  }

  if (base::length(value) == 1L) {
    return(RTMB::advector(rep(value, length)))
  }
  if (base::length(value) != length) {
    stop("'value' must be length 1 or have the requested vector length.", call. = FALSE)
  }
  RTMB::advector(value)
}

.rtmb_auto_ad_seed <- function(envir) {
  if (exists(".rtmb_ad_seed", envir = envir, inherits = FALSE)) {
    seed <- get(".rtmb_ad_seed", envir = envir, inherits = FALSE)
    if (!is.null(seed)) return(seed)
  }
  NULL
}

#' Create an AD-compatible array
#'
#' @description
#' Creates a fixed-size vector, matrix, or higher-dimensional array that is safe
#' to use inside `rtmb_code()` when subsequent assignments may involve RTMB
#' automatic-differentiation values. This is a safer alternative to
#' `array(value, dim)` or `matrix(value, ...)` in model code.
#'
#' @param value Initial value. Usually a scalar such as `0`.
#' @param dim Integer dimension vector.
#' @param seed Optional AD value used to create the container. See
#'   `rtmb_vector()`.
#' @return An AD-compatible array.
#' @export
rtmb_array <- function(value = 0, dim, seed = NULL) {
  if (missing(dim)) {
    stop("'dim' must be supplied.", call. = FALSE)
  }
  if (!is.numeric(dim) || any(is.na(dim)) || any(dim < 0)) {
    stop("'dim' must be a non-negative numeric vector.", call. = FALSE)
  }
  dim <- as.integer(dim)
  len <- prod(dim)
  out <- rtmb_vector(value, len, seed = seed)
  dim(out) <- dim
  out
}

#' Inverse logit function
#'
#' @param x A numeric vector.
#' @return The inverse logit of x.
#' @export
inv_logit <- function(x) plogis(x)

#' Logit function
#'
#' @param x A numeric vector of probabilities.
#' @return The logit of x.
#' @export
logit <- function(x) qlogis(x)

#' Euclidean distance
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param eps A small value added for numerical stability (default is 1e-8).
#' @return The Euclidean distance between x and y.
#' @export
distance <- function(x, y, eps = 1e-8) {
  sqrt(sum((x - y)^2) + eps)
}

#' Squared Euclidean distance
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param eps A small value added for numerical stability (default is 1e-8).
#' @return The squared Euclidean distance between x and y.
#' @export
squared_distance <- function(x, y, eps = 1e-8) {
  sum((x - y)^2) + eps
}

#' Log-sum-exp function
#'
#' @param x A numeric vector or scalar.
#' @param y An optional numeric scalar.
#' @return The log of the sum of the exponentials.
#' @export
log_sum_exp <- function(x, y = NULL) {
  if (is.null(y)) {
    # Check if x is a matrix
    if (is.matrix(x)) {
      return(log_sum_exp_matrix(x))
    }
    max_val <- max(x)
    return(max_val + log(sum(exp(x - max_val))))
  } else {
    max_val <- (x + y + abs(x - y)) / 2
    return(max_val + log(exp(x - max_val) + exp(y - max_val)))
  }
}

#' Log-sum-exp function for matrices (row-wise)
#'
#' @param M A numeric matrix or advector matrix.
#' @return A numeric vector of row-wise log-sum-exp.
#' @export
log_sum_exp_matrix <- function(M) {
  # This function hides 'apply' from the BayesRTMB parser
  max_val <- apply(M, 1, max)
  return(max_val + log(rowSums(exp(M - max_val))))
}

#' Log of one minus x
#'
#' @param x A numeric vector.
#' @return The natural logarithm of 1 - x.
#' @export
log1m <- function(x) {
  log1p(-x)
}

#' Log of one plus exponential of x
#'
#' @param x A numeric vector.
#' @return The natural logarithm of 1 + exp(x).
#' @export
log1p_exp <- function(x) {
  max_val <- (x + sqrt(x^2 + 1e-10)) / 2
  return(max_val + log(exp(x - max_val) + exp(-max_val)))
}

#' Log of one minus exponential of x
#'
#' @param x A numeric vector.
#' @return The natural logarithm of 1 - exp(x).
#' @export
log1m_exp <- function(x) {
  # Stable implementation of log(1 - exp(x)) for x < 0
  ifelse(x > -0.6931472, log(-expm1(x)), log1p(-exp(x)))
}

#' Log mixture of two probabilities
#'
#' @param theta Mixing proportion (between 0 and 1).
#' @param lp1 Log-probability of the first component.
#' @param lp2 Log-probability of the second component.
#' @return The log-probability of the mixture.
#' @export
log_mix <- function(theta, lp1, lp2) {
  l1 <- log(theta) + lp1
  l2 <- log1p(-theta) + lp2
  return(log_sum_exp(l1, l2))
}

#' Softmax function
#'
#' @param x A numeric vector.
#' @return A numeric vector of softmax probabilities.
#' @export
softmax <- function(x) {
  max_x <- max(x)
  exp_x <- exp(x - max_x)
  return(exp_x / sum(exp_x))
}

#' Log-softmax function
#'
#' @param x A numeric vector.
#' @return A numeric vector of log-softmax probabilities.
#' @export
log_softmax <- function(x) {
  return(x - log_sum_exp(x))
}

#' Smooth absolute value function
#'
#' @param x A numeric vector.
#' @param eps A small value added for numerical stability (default is 1e-8).
#' @return The smooth absolute value of x.
#' @export
fabs <- function(x, eps = 1e-8) {
  sqrt(x^2 + eps)
}

#' Log determinant of a Cholesky factor
#'
#' @param L A lower triangular Cholesky factor matrix.
#' @return The log determinant of the corresponding covariance matrix.
#' @export
log_det_chol <- function(L) {
  return(2 * sum(log(diag(L))))
}

#' Quadratic form using a Cholesky factor
#'
#' @param x A numeric vector.
#' @param L A lower triangular Cholesky factor matrix.
#' @return The value of the quadratic form.
#' @export
quad_form_chol <- function(x, L) {
  z <- solve(L, x)
  return(sum(z^2))
}

#' Quadratic form with a diagonal matrix
#'
#' @param m A numeric matrix.
#' @param v A numeric vector representing the diagonal elements.
#' @return The value of the quadratic form.
#' @export
quad_form_diag <- function(m, v) {
  v_col <- matrix(as.vector(v), ncol = 1)
  return(m * (v_col %*% t(v_col)))
}
#' stz basis function
#'
#' @param K A numeric value
#' @return A transformation matrix
#' @export
stz_basis <- function(K) {
  Q <- matrix(0, nrow = K, ncol = K - 1)
  for (j in 1:(K - 1)) {
    Q[1:j, j] <- 1 / sqrt(j * (j + 1))
    Q[j + 1, j] <- -j / sqrt(j * (j + 1))
  }
  return(Q)
}

#' Sum-to-zero transformation
#'
#' Transforms a vector of length K-1 to a vector of length K that sums to zero
#' using an orthogonal basis matrix.
#'
#' @param x A numeric vector of length K-1.
#' @return A numeric vector of length K whose elements sum to zero.
#' @export
sum_to_zero <- function(x) {
  K <- length(x) + 1
  Q <- stz_basis(K)
  return(as.vector(Q %*% x))
}

#' Vector to lower triangular matrix (RTMB compatible)
#'
#' @param x A numeric or advector.
#' @param M Number of rows.
#' @param D Number of columns.
#' @return An M x D lower triangular matrix.
to_lower_tri <- function(x, M, D) {
  y <- matrix(x[1] * 0, nrow = M, ncol = D)
  y[lower.tri(y, diag = TRUE)] <- x
  return(y)
}
#' Vector to centered matrix (RTMB compatible)
#'
#' @param x A numeric vector of length (R-1) * C.
#' @param R Number of rows.
#' @param C Number of columns.
#' @return An R x C matrix where each column sums to zero.
#' @export
to_centered_matrix <- function(x, R, C) {
  mat_unc <- matrix(x[1] * 0, nrow = R - 1, ncol = C)
  mat_unc[] <- x
  Q <- stz_basis(R)
  return(Q %*% mat_unc)
}
#' Vector to centered triangular matrix (RTMB compatible)
#'
#' @param x A numeric vector of appropriate length.
#' @param R Number of rows.
#' @param C Number of columns.
#' @return An R x C matrix with column-wise sum-to-zero constraints on lower elements.
#' @export
to_centered_tri <- function(x, R, C) {
  y <- matrix(x[1] * 0, nrow = R, ncol = C)
  max_d <- min(C, R - 1)
  pos <- 1

  for (d in 1:max_d) {
    n_d <- R - d + 1
    n_params <- n_d - 1

    if (pos + n_params - 1 <= length(x)) {
      z_d <- x[pos:(pos + n_params - 1)]
      Q_d <- stz_basis(n_d)
      y[d:R, d] <- Q_d %*% z_d
      pos <- pos + n_params
    }
  }

  return(y)
}
