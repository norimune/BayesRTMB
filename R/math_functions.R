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
    max_val <- max(as.numeric(x))
    return(max_val + log(sum(exp(x - max_val))))
  } else {
    max_val <- (x + y + abs(x - y)) / 2
    return(max_val + log(exp(x - max_val) + exp(y - max_val)))
  }
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
  return(log(1 - exp(x) + 1e-10))
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
  l2 <- log(1 - theta) + lp2
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
#' Sum-to-zero transformation
#'
#' Transforms a vector of length K-1 to a vector of length K that sums to zero
#' using an orthogonal basis matrix.
#'
#' @param x A numeric vector of length K-1.
#' @return A numeric vector of length K whose elements sum to zero.
#' @export
sum_to_zero <- function(x) {
  # パラメータ数 K を特定 (入力は K-1 次元)
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
  pos <- 1
  for (i in 1:M) {
    max_j <- min(i, D)
    for (j in 1:max_j) {
      y[i, j] <- x[pos]
      pos <- pos + 1
    }
  }
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
  # 1. 無制約ベクトルを (R-1) x C の行列に整形
  # RTMB::advector(0) で初期化して型を合わせる
  mat_unc <- matrix(RTMB::advector(0), nrow = R - 1, ncol = C)
  mat_unc[] <- x
  Q <- stz_basis(R)

  # 結果も advector 属性を持つ行列として返される
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
  y <- matrix(RTMB::advector(0), nrow = R, ncol = C)
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
