test_that("normalized centered priors do not receive extra Jacobian constants", {
  code <- rtmb_code(
    parameters = { x <- Dim(4, type = "centered") },
    model = { x ~ centered_multi_normal(1) }
  )
  model <- rtmb_model(list(), code, init = list(x = rep(0, 4)), silent = TRUE)
  ad <- model$get_ad_obj(jacobian_target = "all")
  u <- c(0.3, -0.5, 0.7)
  expect_equal(-ad$fn(u), sum(dnorm(u, log = TRUE)), tolerance = 1e-12)

  for (type in c("centered_matrix", "centered_tri", "positive_centered_tri")) {
    p <- list(x = Dim(c(4, 2), type = type))
    u <- seq(-0.4, 0.4, length.out = p$x$unc_length)
    transform <- function(v) as.numeric(to_constrained(list(x = v), p)$x)
    J <- sapply(seq_along(u), function(i) {
      d <- numeric(length(u))
      d[i] <- 1e-6
      (transform(u + d) - transform(u - d)) / 2e-6
    })
    expected <- as.numeric(determinant(crossprod(J))$modulus) / 2
    expect_equal(calc_log_jacobian(list(x = u), p), expected, tolerance = 1e-8)
  }
})

test_that("a tied map counts the Jacobian once per shared bounded coordinate", {
  code <- rtmb_code(
    parameters = { x <- Dim(2, lower = 0) },
    model = { x ~ exponential(1) }
  )
  model <- rtmb_model(list(), code, init = list(x = c(1, 1)), silent = TRUE)
  map <- list(x = factor(c(1, 1)))
  ad <- model$get_ad_obj(jacobian_target = "all", map = map)
  u <- log(2)
  expect_equal(-ad$fn(u), -2 * exp(u) + u, tolerance = 1e-12)

  values <- list(x = rep(u, 2))
  expect_equal(calc_log_jacobian(values, list(x = Dim(2, upper = 0)), map = map), u)
  expect_equal(
    calc_log_jacobian(values, list(x = Dim(2, lower = 0, upper = 1)), map = map),
    log(plogis(u) * plogis(-u))
  )
  expect_error(
    calc_log_jacobian(values, list(x = Dim(2, lower = c(0, 1))), map = map),
    "identical bounds"
  )
  expect_error(
    calc_log_jacobian(values, list(x = Dim(2, type = "ordered")), map = map),
    "tied map entries is not supported"
  )
})

test_that("structured fixed values are validated and cannot be fixed element-wise", {
  code <- rtmb_code(
    parameters = {
      cuts <- Dim(3, type = "ordered")
      L <- Dim(c(3, 3), type = "positive_lower_tri")
      theta <- Dim(3, type = "simplex")
    },
    model = {
      cuts ~ normal(0, 2)
      L ~ normal(0, 1)
      theta ~ dirichlet(c(2, 2, 2))
    }
  )
  model <- rtmb_model(
    list(), code, init = list(cuts = c(0, 1, 2), L = diag(3), theta = c(0.2, 0.3, 0.5)),
    silent = TRUE
  )
  expect_error(model$fixed_model(list("cuts[2]" = 1.5)), "Fix the entire parameter")
  expect_error(model$fixed_model(list("L[5]" = 1.5)), "Fix the entire parameter")
  expect_error(model$fixed_model(list(theta = c(0.2, 0.3, 0.9))), "simplex.*constraint")
  expect_error(model$fixed_model(list(theta = c(0.2, 0.8))), "Expected 3 fixed values")
  fixed <- model$fixed_model(list(cuts = c(0, 1.5, 2), theta = c(0.1, 0.3, 0.6)), silent = TRUE)
  ad <- fixed$get_ad_obj()
  x <- fixed$to_constrained(ad$env$parList())
  expect_equal(unname(x$cuts), c(0, 1.5, 2))
  expect_equal(unname(x$theta), c(0.1, 0.3, 0.6))
})

test_that("boundary transformations preserve representable interior values", {
  p <- list(x = Dim(3, type = "simplex"))
  x <- list(x = c(1 - 1e-12, 2e-13, 8e-13))
  expect_equal(log(to_constrained(to_unconstrained(x, p), p)$x), log(x$x), tolerance = 1e-12)

  for (type in c("cov_matrix", "CF_corr")) {
    p <- list(x = Dim(c(2, 2), type = type))
    u <- list(x = if (type == "cov_matrix") c(log(1e-6), 0, 0) else 20)
    expect_equal(to_unconstrained(to_constrained(u, p), p), u, tolerance = 1e-12)
  }

  p <- list(x = Dim(lower = 0, upper = 1))
  u <- list(x = -710)
  expect_equal(to_unconstrained(to_constrained(u, p), p), u, tolerance = 1e-12)
  expect_equal(calc_log_jacobian(u, p), -710)
  ad <- RTMB::MakeADFun(function(v) calc_log_jacobian(v, p), list(x = 0), silent = TRUE)
  expect_equal(ad$fn(-710), -710)
  expect_equal(as.numeric(ad$gr(-710)), 1)
})

test_that("constant and singleton dimensions are handled consistently", {
  expect_equal(to_constrained(list(x = numeric(0)), list(x = Dim(1, type = "simplex")))$x, 1)
  expect_equal(to_constrained(list(x = numeric(0)), list(x = Dim(1, type = "centered")))$x, 0)
  expect_equal(dim(stz_basis(1)), c(1L, 0L))
  for (type in c("centered_matrix", "centered_tri", "positive_centered_tri")) {
    p <- list(x = Dim(c(1, 2), type = type))
    u <- list(x = numeric(0))
    expect_equal(to_constrained(u, p)$x, matrix(0, 1, 2))
    expect_equal(to_unconstrained(list(x = matrix(0, 1, 2)), p), u)
    expect_error(Dim(c(2, 3, 3), type = type), "requires two dimensions")
  }
  for (type in c("corr_matrix", "CF_corr", "cov_matrix", "CF_cov", "lower_tri_stz")) {
    p <- list(x = Dim(c(1, 3, 3), type = type))
    u <- list(x = seq(-0.4, 0.6, length.out = p$x$unc_length))
    x <- to_constrained(u, p)
    expect_equal(dim(x$x), c(1, 3, 3))
    expect_equal(to_unconstrained(x, p), u, tolerance = 1e-12)
  }
})
