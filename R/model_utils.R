#' @title Parameter Conversion and Utility Functions
#'
#' @description
#' These functions provide tools for converting parameters between constrained (natural)
#' and unconstrained (mathematical) spaces, building RTMB objective functions,
#' and performing numerical differentiation.
#'
#' @name model_utils
NULL

#' @rdname model_utils
#' @description Convert a list of constrained parameters to unconstrained space.
#' @param para_orig_list A named list of parameters in their natural space.
#' @param par_list The parameter definition list (from \code{RTMB_Model$par_list}).
#' @return A named list of parameters in unconstrained space.
#' @export
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
      if (is.matrix(val_orig)) {
        y <- val_orig
        y[, 1] <- val_orig[, 1]
        if (ncol(val_orig) > 1) {
          for (k in 2:ncol(val_orig)) {
            y[, k] <- log(val_orig[, k] - val_orig[, k - 1])
          }
        }
        para_unc[[name]] <- y
      } else {
        y <- numeric(p$length)
        y[1] <- val_orig[1]
        if (p$length > 1) {
          for (i in 2:p$length) {
            y[i] <- log(val_orig[i] - val_orig[i-1])
          }
        }
        para_unc[[name]] <- y
      }
    } else if (b_type == "positive_ordered") {
      if (is.matrix(val_orig)) {
        y <- val_orig
        y[, 1] <- log(val_orig[, 1])
        if (ncol(val_orig) > 1) {
          for (k in 2:ncol(val_orig)) {
            y[, k] <- log(val_orig[, k] - val_orig[, k - 1])
          }
        }
        para_unc[[name]] <- y
      } else {
        y <- numeric(p$length)
        y[1] <- log(val_orig[1])
        if (p$length > 1) {
          for (i in 2:p$length) {
            y[i] <- log(val_orig[i] - val_orig[i-1])
          }
        }
        para_unc[[name]] <- y
      }
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
    } else if (b_type == "centered") {
      K <- p$length
      Q <- stz_basis(K)
      para_unc[[name]] <- as.numeric(crossprod(Q, val_orig))

    } else if (b_type == "centered_matrix") {
      R <- p$dim[1]
      C <- p$dim[2]
      Q <- stz_basis(R)
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
      if (length(p$dim) == 3) {
        K <- p$dim[1]; R <- p$dim[2]; C <- p$dim[3]
      } else {
        R <- p$dim[1]; C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]; K <- 1
      }
      z_all <- numeric(p$unc_length)
      idx <- 1
      for (k in 1:K) {
        val_slice <- if (K > 1) {
          if (length(p$dim) == 3) matrix(val_orig[k, , ], R, C) else val_orig
        } else {
          val_orig
        }
        L <- if (b_type == "corr_matrix") t(chol(val_slice + diag(1e-8, R))) else val_slice
        for (i in 2:R) {
          prod_term <- 1
          max_j <- min(i - 1, C - 1)
          if (max_j > 0) {
            for (j in 1:max_j) {
              z_all[idx] <- L[i, j] / sqrt(prod_term)
              prod_term <- prod_term * (1 - z_all[idx]^2)
              idx <- idx + 1
            }
          }
        }
      }
      para_unc[[name]] <- atanh(pmin(pmax(z_all, -0.9999), 0.9999))

    } else if (b_type %in% c("cov_matrix", "CF_cov")) {
      if (length(p$dim) == 3) {
        K <- p$dim[1]; R <- p$dim[2]; C <- p$dim[3]
      } else {
        R <- p$dim[1]; C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]; K <- 1
      }
      y <- numeric(p$unc_length)
      idx <- 1
      for (k in 1:K) {
        val_slice <- if (K > 1) {
          if (length(p$dim) == 3) matrix(val_orig[k, , ], R, C) else val_orig
        } else {
          val_orig
        }
        L <- if (b_type == "cov_matrix") t(chol(val_slice + diag(1e-8, R))) else val_slice
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
      }
      para_unc[[name]] <- y

    } else if (b_type %in% c("lower_tri", "positive_lower_tri")) {
      if (length(p$dim) == 3) {
        K <- p$dim[1]; R <- p$dim[2]; C <- p$dim[3]
      } else {
        R <- p$dim[1]; C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]; K <- 1
      }
      y <- numeric(p$unc_length)
      idx <- 1
      for (k in 1:K) {
        val_slice <- if (K > 1) {
          if (length(p$dim) == 3) matrix(val_orig[k, , ], R, C) else val_orig
        } else {
          val_orig
        }
        for (i in 1:R) {
          max_j <- min(i, C)
          for (j in 1:max_j) {
            y[idx] <- val_slice[i, j]
            idx <- idx + 1
          }
        }
      }
      para_unc[[name]] <- y

    } else if (b_type == "lower_tri_stz") {
      if (length(p$dim) == 3) {
        K <- p$dim[1]; R <- p$dim[2]; C <- p$dim[3]
      } else {
        R <- p$dim[1]; C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]; K <- 1
      }
      y <- numeric(p$unc_length)
      idx <- 1
      for (k in 1:K) {
        val_slice <- if (K > 1) {
          if (length(p$dim) == 3) matrix(val_orig[k, , ], R, C) else val_orig
        } else {
          val_orig
        }
        last_i <- R
        last_j <- min(R, C)
        for (i in 1:R) {
          max_j <- min(i, C)
          for (j in 1:max_j) {
            if (i == last_i && j == last_j) break
            y[idx] <- val_slice[i, j]
            idx <- idx + 1
          }
        }
      }
      para_unc[[name]] <- y
    } else {
      para_unc[[name]] <- as.numeric(val_orig)
    }
  }
  return(para_unc)
}

#' @rdname model_utils
#' @description Convert a list of unconstrained parameters to constrained space.
#' @param para_unc_list A named list of parameters in unconstrained space.
#' @param par_list The parameter definition list.
#' @return A named list of parameters in constrained space.
#' @export
to_constrained <- function(para_unc_list, par_list) {
  para <- para_unc_list
  for (name in names(par_list)) {
    p <- par_list[[name]]
    val_unc <- para_unc_list[[name]]
    b_type <- p$bounds

    ad_zero <- if (length(val_unc) > 0) val_unc[1] * 0 else 0

    if (b_type == "lower") {
      para[[name]] <- p$lower + exp(val_unc)
    } else if (b_type == "upper") {
      para[[name]] <- p$upper - exp(val_unc)
    } else if (b_type == "interval") {
      prob <- plogis(val_unc)
      para[[name]] <- p$lower + (p$upper - p$lower) * prob
    } else if (b_type == "ordered") {
      if (p$length > 1) {
        if (is.matrix(val_unc)) {
          res <- val_unc
          res[, 1] <- val_unc[, 1]
          if (ncol(val_unc) > 1) {
            for (k in 2:ncol(val_unc)) {
              res[, k] <- res[, k - 1] + exp(val_unc[, k])
            }
          }
          para[[name]] <- res
        } else {
          para[[name]] <- val_unc[1] + cumsum(c(ad_zero, exp(val_unc[-1])))
        }
      } else {
        para[[name]] <- val_unc
      }
    } else if (b_type == "positive_ordered") {
      if (is.matrix(val_unc)) {
        res <- val_unc
        res[, 1] <- exp(val_unc[, 1])
        if (ncol(val_unc) > 1) {
          for (k in 2:ncol(val_unc)) {
            res[, k] <- res[, k - 1] + exp(val_unc[, k])
          }
        }
        para[[name]] <- res
      } else {
        para[[name]] <- cumsum(exp(val_unc))
      }
    } else if (b_type == "simplex") {
      K <- p$length
      z <- 1 / (1 + exp(-(val_unc - log(1 / (K - seq_len(K - 1))))))
      ad_one <- ad_zero + 1
      rem_prev <- c(ad_one, cumprod(1 - z)[- (K - 1)])
      x_first <- rem_prev * z
      x_last <- cumprod(1 - z)[K - 1]
      para[[name]] <- c(x_first, x_last)

    } else if (b_type == "centered") {
      K <- p$length
      Q <- stz_basis(K)
      para[[name]] <- as.numeric(Q %*% val_unc)

    } else if (b_type == "centered_matrix") {
      R <- p$dim[1]; C <- p$dim[2]
      Q <- stz_basis(R)
      mat_unc <- matrix(val_unc, nrow = R - 1, ncol = C)
      para[[name]] <- Q %*% mat_unc

    } else if (b_type == "centered_tri") {
      R <- p$dim[1]; C <- p$dim[2]
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
      R <- p$dim[1]; C <- p$dim[2]
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
      if (length(p$dim) == 3) {
        K <- p$dim[1]; R <- p$dim[2]; C <- p$dim[3]
      } else {
        R <- p$dim[1]; C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]; K <- 1
      }
      z_all <- tanh(val_unc)

      if (K > 1) {
        res_list <- vector("list", K)
        unc_per_slice <- length(val_unc) / K
        for (k in 1:K) {
          z <- z_all[((k-1)*unc_per_slice + 1):(k*unc_per_slice)]
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
          if (b_type == "corr_matrix") {
            res_list[[k]] <- as.vector(L %*% t(L))
          } else {
            res_list[[k]] <- as.vector(L)
          }
        }
        res_array_flat <- do.call(c, res_list)
        idx_base <- array(1:(K * R * C), dim = c(R, C, K))
        idx_perm <- as.vector(aperm(idx_base, c(3, 1, 2)))

        res_array <- res_array_flat[idx_perm]
        dim(res_array) <- p$dim
        para[[name]] <- res_array
      } else {
        z <- z_all
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
      }

    } else if (b_type %in% c("cov_matrix", "CF_cov")) {
      if (length(p$dim) == 3) {
        K <- p$dim[1]; R <- p$dim[2]; C <- p$dim[3]
      } else {
        R <- p$dim[1]; C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]; K <- 1
      }

      if (K > 1) {
        res_list <- vector("list", K)
        unc_per_slice <- length(val_unc) / K
        for (k in 1:K) {
          v_u <- val_unc[((k-1)*unc_per_slice + 1):(k*unc_per_slice)]
          L <- matrix(ad_zero, nrow = R, ncol = C)
          idx <- 1
          for (i in 1:R) {
            max_j <- min(i, C)
            for (j in 1:max_j) {
              if (i == j) L[i, j] <- exp(v_u[idx])
              else L[i, j] <- v_u[idx]
              idx <- idx + 1
            }
          }
          if (b_type == "cov_matrix") {
            res_list[[k]] <- as.vector(L %*% t(L))
          } else {
            res_list[[k]] <- as.vector(L)
          }
        }
        res_array_flat <- do.call(c, res_list)
        idx_base <- array(1:(K * R * C), dim = c(R, C, K))
        idx_perm <- as.vector(aperm(idx_base, c(3, 1, 2)))

        res_array <- res_array_flat[idx_perm]
        dim(res_array) <- p$dim
        para[[name]] <- res_array
      } else {
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
      }

    } else if (b_type %in% c("lower_tri", "positive_lower_tri")) {
      if (length(p$dim) == 3) {
        K <- p$dim[1]; R <- p$dim[2]; C <- p$dim[3]
      } else {
        R <- p$dim[1]; C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]; K <- 1
      }

      if (K > 1) {
        res_list <- vector("list", K)
        unc_per_slice <- length(val_unc) / K
        for (k in 1:K) {
          v_u <- val_unc[((k-1)*unc_per_slice + 1):(k*unc_per_slice)]
          L <- matrix(ad_zero, nrow = R, ncol = C)
          idx <- 1
          for (i in 1:R) {
            max_j <- min(i, C)
            for (j in 1:max_j) {
              L[i, j] <- v_u[idx]
              idx <- idx + 1
            }
          }
          res_list[[k]] <- as.vector(L)
        }
        res_array_flat <- do.call(c, res_list)
        idx_base <- array(1:(K * R * C), dim = c(R, C, K))
        idx_perm <- as.vector(aperm(idx_base, c(3, 1, 2)))

        res_array <- res_array_flat[idx_perm]
        dim(res_array) <- p$dim
        para[[name]] <- res_array
      } else {
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
      }

    } else if (b_type == "lower_tri_stz") {
      if (length(p$dim) == 3) {
        K <- p$dim[1]; R <- p$dim[2]; C <- p$dim[3]
      } else {
        R <- p$dim[1]; C <- if (length(p$dim) > 1) p$dim[2] else p$dim[1]; K <- 1
      }

      if (K > 1) {
        res_list <- vector("list", K)
        unc_per_slice <- length(val_unc) / K
        for (k in 1:K) {
          v_u <- val_unc[((k-1)*unc_per_slice + 1):(k*unc_per_slice)]
          L <- matrix(ad_zero, nrow = R, ncol = C)
          idx <- 1
          sum_val <- 0
          last_i <- R
          last_j <- min(R, C)
          for (i in 1:R) {
            max_j <- min(i, C)
            for (j in 1:max_j) {
              if (i == last_i && j == last_j) break
              L[i, j] <- v_u[idx]
              sum_val <- sum_val + v_u[idx]
              idx <- idx + 1
            }
          }
          L[last_i, last_j] <- -sum_val
          res_list[[k]] <- as.vector(L)
        }
        res_array_flat <- do.call(c, res_list)
        idx_base <- array(1:(K * R * C), dim = c(R, C, K))
        idx_perm <- as.vector(aperm(idx_base, c(3, 1, 2)))

        res_array <- res_array_flat[idx_perm]
        dim(res_array) <- p$dim
        para[[name]] <- res_array
      } else {
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
    } else {
      para[[name]] <- val_unc
    }

    if (!(b_type %in% c("corr_matrix", "CF_corr", "cov_matrix", "CF_cov", "lower_tri", "lower_tri_stz", "positive_lower_tri")) && length(p$dim) > 1) {
      dim(para[[name]]) <- p$dim
    }
  }
  return(para)
}

#' @rdname model_utils
#' @description Convert a flat unconstrained vector to a named list.
#' @param vec A numeric vector in unconstrained space.
#' @param par_list The parameter definition list.
#' @return A named list of unconstrained parameters.
#' @export
unconstrained_vector_to_list <- function(vec, par_list) {
  res <- list()
  idx <- 1
  for (name in names(par_list)) {
    p <- par_list[[name]]
    len <- p$unc_length
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

#' @rdname model_utils
#' @description Convert a flat constrained vector to a named list.
#' @param vec A numeric vector in constrained space.
#' @param par_list The parameter definition list.
#' @return A named list of constrained parameters.
#' @export
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

#' @rdname model_utils
#' @description Alias for to_unconstrained (for backward compatibility).
#' @export
to_unconstrained_list <- to_unconstrained

#' @rdname model_utils
#' @description Alias for to_constrained (for backward compatibility).
#' @export
to_constrained_list <- to_constrained

#' Simple Numerical Jacobian
#'
#' @description Calculates the numerical Jacobian matrix of a function using central differences.
#' @param f A function that takes a numeric vector and returns a numeric vector.
#' @param x The point at which to evaluate the Jacobian.
#' @param h The relative step size for finite differences. Default is 1e-4.
#' @return A matrix of partial derivatives.
#' @keywords internal
simple_jacobian <- function(f, x, h = 1e-4) {
  p <- length(x)
  f0 <- f(x)
  res <- matrix(0, nrow = length(f0), ncol = p)
  for (i in 1:p) {
    h_i <- max(abs(x[i]), 1e-4) * h
    x_plus <- x_minus <- x
    x_plus[i] <- x[i] + h_i
    x_minus[i] <- x[i] - h_i
    res[, i] <- (f(x_plus) - f(x_minus)) / (2 * h_i)
  }
  return(res)
}

#' Factory function to create the internal AD objective function
#'
#' @description Creates a closure that can be passed to parallel workers
#' without carrying the entire R6 object.
#' @param data Model data.
#' @param par_list Parameter definitions.
#' @param log_prob User log-probability function.
#' @param transform User transform function.
#' @param jacobian_target Character; "all", "random", or "none".
#' @return A function of one argument (unconstrained parameter list).
#' @keywords internal
.make_f_ad <- function(data, par_list, log_prob, transform = NULL, jacobian_target = "all") {
  data <- data
  par_list <- par_list
  log_prob <- log_prob
  transform <- transform
  jacobian_target <- jacobian_target

  function(y_unc_list) {
    para <- to_constrained(y_unc_list, par_list)

    if (!is.null(transform)) {
      tran_res <- transform(data, para)
      para <- c(para, tran_res)
      
      if (length(tran_res) > 0) {
        for (n in names(tran_res)) {
          assign(n, tran_res[[n]])
          eval(substitute(RTMB::ADREPORT(V), list(V = as.name(n))))
        }
      }
    }

    lp <- log_prob(data, para)

    if (jacobian_target == "all") {
      lj <- calc_log_jacobian(y_unc_list, par_list, only_random = FALSE)
      return(-(lp + lj))
    } else if (jacobian_target == "random") {
      lj <- calc_log_jacobian(y_unc_list, par_list, only_random = TRUE)
      return(-(lp + lj))
    } else {
      return(-lp)
    }
  }
}
