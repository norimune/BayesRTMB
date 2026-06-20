test_that("rtmb_vector and rtmb_array create AD-compatible containers", {
  v <- rtmb_vector(0, 3)
  a <- rtmb_array(0, dim = c(2, 3))

  expect_true(methods::is(v, "advector"))
  expect_true(methods::is(a, "advector"))
  expect_equal(length(v), 3)
  expect_equal(dim(a), c(2L, 3L))
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
