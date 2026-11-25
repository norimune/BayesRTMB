test_that("matrix parameter Jacobians agree with numerical derivatives", {
  specs <- list(
    list(type = "corr_matrix", dim = c(3, 3)),
    list(type = "corr_matrix", dim = c(2, 4, 4)),
    list(type = "cov_matrix", dim = c(2, 2)),
    list(type = "cov_matrix", dim = c(2, 3, 3)),
    list(type = "CF_corr", dim = c(4, 3)),
    list(type = "CF_cov", dim = c(4, 3)),
    list(type = "positive_lower_tri", dim = c(3, 2))
  )
  for (spec in specs) {
    p <- list(x = Dim(spec$dim, type = spec$type))
    u <- seq(-0.7, 0.8, length.out = p$x$unc_length)
    forward <- function(u) {
      x <- to_constrained(list(x = u), p)$x
      diagonal <- !spec$type %in% c("corr_matrix", "CF_corr")
      elements <- function(A) {
        keep <- lower.tri(A, diag = diagonal)
        if (spec$type == "CF_corr") keep <- keep & col(A) < ncol(A)
        A[keep]
      }
      if (length(spec$dim) == 2L) return(elements(x))
      unlist(lapply(seq_len(spec$dim[1]), function(k) elements(x[k, , ])))
    }
    step <- 1e-6
    J <- sapply(seq_along(u), function(i) {
      d <- numeric(length(u))
      d[i] <- step
      (forward(u + d) - forward(u - d)) / (2 * step)
    })
    expected <- as.numeric(determinant(J, logarithm = TRUE)$modulus)
    expect_equal(calc_log_jacobian(list(x = u), p), expected, tolerance = 1e-8)
  }
})

test_that("positive lower triangular transforms preserve shape and round-trip", {
  for (dims in list(c(2, 2), c(1, 3, 2), c(2, 3, 2))) {
    p <- list(x = Dim(dims, type = "positive_lower_tri"))
    u <- list(x = seq(-1, -0.2, length.out = p$x$unc_length))
    x <- to_constrained(u, p)
    expect_equal(dim(x$x), dims)
    if (length(dims) == 2L) {
      expect_true(all(diag(x$x) > 0))
    } else {
      expect_true(all(x$x[, 1, 1] > 0 & x$x[, 2, 2] > 0))
    }
    expect_equal(to_unconstrained(x, p), u, tolerance = 1e-12)

    objective <- function(par) sum(to_constrained(par, p)$x)
    ad <- RTMB::MakeADFun(objective, u, silent = TRUE)
    expect_equal(ad$fn(ad$par), objective(u), tolerance = 1e-12)
    expect_true(all(is.finite(ad$gr(ad$par))))
  }
})
