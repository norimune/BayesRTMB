test_that("rtmb_vector and rtmb_array create AD-compatible containers", {
  v <- rtmb_vector(0, 3)
  a <- rtmb_array(0, dim = c(2, 3))
  v_seed <- rtmb_vector(0, 3, seed = RTMB::advector(0))
  a_seed <- rtmb_array(0, dim = c(2, 3), seed = RTMB::advector(0))

  expect_true(methods::is(v, "advector"))
  expect_true(methods::is(a, "advector"))
  expect_true(methods::is(v_seed, "advector"))
  expect_true(methods::is(a_seed, "advector"))
  expect_equal(length(v), 3)
  expect_equal(dim(a), c(2L, 3L))
  expect_equal(length(v_seed), 3)
  expect_equal(dim(a_seed), c(2L, 3L))
})

test_that("rtmb containers can be assigned to inside rtmb_code", {
  dat <- list(
    N = 4L,
    M = 3L,
    y = c(0.1, -0.2, 0.3, 0.4)
  )

  code <- rtmb_code(
    parameters = {
      theta <- Dim()
    },
    model = {
      x <- rtmb_vector(0, N)
      A <- rtmb_array(0, dim = c(N, M))
      for (i in 1:N) {
        x[i] <- theta + i * 0.01
        for (j in 1:M) {
          A[i, j] <- theta + i * 0.01 + j * 0.02
        }
      }
      y ~ normal(x, 1)
      theta ~ normal(0, 1)
    }
  )

  expect_s3_class(rtmb_model(dat, code, silent = TRUE), "RTMB_Model")
})

test_that("softmax helpers support AD containers", {
  x_num <- c(0.2, -0.1, 0.4)
  x <- RTMB::advector(x_num)

  expect_equal(softmax(x_num), exp(x_num) / sum(exp(x_num)), tolerance = 1e-8)
  expect_equal(log_softmax(x_num), log(exp(x_num) / sum(exp(x_num))), tolerance = 1e-8)
  expect_true(methods::is(log_sum_exp(x), "advector"))
  expect_true(methods::is(softmax(x), "advector"))
  expect_true(methods::is(log_softmax(x), "advector"))
  expect_true(methods::is(.rtmb_c(0, x), "advector"))
  expect_true(methods::is(log_sum_exp(.rtmb_c(0, x)), "advector"))
  expect_true(methods::is(softmax(.rtmb_c(0, x)), "advector"))

  A <- rtmb_array(0, dim = c(2, 3), seed = RTMB::advector(0))
  A[1, ] <- RTMB::advector(c(0.1, 0.2, 0.3))
  A[2, ] <- RTMB::advector(c(-0.1, 0.4, 0.6))
  expect_true(methods::is(log_sum_exp(A), "advector"))
  expect_length(log_sum_exp(A), 2)
})

test_that("softmax helpers support baseline categories inside rtmb_code", {
  dat <- list(
    N = 3L,
    K = 3L,
    K_prob = 2L,
    X = matrix(c(0, 1, 0, 1, 0, 1), nrow = 3, ncol = 2),
    y = c(1L, 2L, 3L)
  )

  code <- rtmb_code(
    setup = {
      N <- N
      K <- K
      K_prob <- K_prob
      X <- X
      y <- y
    },
    parameters = {
      z <- Dim(c(K_prob, K - 1))
    },
    model = {
      eta <- X %*% z
      lp <- z[1] * 0
      for (i in 1:N) {
        lp <- lp + log_softmax(c(0, eta[i, ]))[y[i]]
      }
      z ~ normal(0, 1)
    }
  )

  expect_s3_class(rtmb_model(dat, code, silent = TRUE), "RTMB_Model")
})

test_that("report selects multiple transformed values", {
  fn <- transform_code({
    a <- 1
    b <- 2
    c <- 3
    report(a, b)
  })

  out <- fn(list(.dummy = 1), list(.par = 2))
  expect_named(out, c("a", "b"))
  expect_equal(out$a, 1)
  expect_equal(out$b, 2)
  expect_false("c" %in% names(out))
})

test_that("namespaced report selects transformed values", {
  fn <- transform_code({
    a <- 1
    b <- 2
    BayesRTMB::report(b)
  })

  out <- fn(list(.dummy = 1), list(.par = 2))
  expect_named(out, "b")
  expect_equal(out$b, 2)
})
