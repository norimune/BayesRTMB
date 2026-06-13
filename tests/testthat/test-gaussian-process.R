test_that("gaussian_process_lpdf vectorizes matrix rows", {
  set.seed(1)
  x <- seq_len(5)
  y <- matrix(rnorm(15), nrow = 3, ncol = 5)

  row_lp <- vapply(seq_len(nrow(y)), function(i) {
    gaussian_process_lpdf(
      y[i, ],
      x = x,
      mean = 0,
      magnitude = 1.2,
      smoothing = 2.1,
      noise = 0.3
    )
  }, numeric(1))

  expect_equal(
    gaussian_process_lpdf(
      y,
      x = x,
      mean = 0,
      magnitude = 1.2,
      smoothing = 2.1,
      noise = 0.3
    ),
    sum(row_lp),
    tolerance = 1e-10
  )

  expect_equal(
    gaussian_process_lpdf(
      y,
      x = x,
      mean = 0,
      magnitude = 1.2,
      smoothing = 2.1,
      noise = 0.3,
      sum = FALSE
    ),
    row_lp,
    tolerance = 1e-10
  )
})

test_that("gaussian_process_lpdf checks coordinate length", {
  y <- matrix(rnorm(6), nrow = 2, ncol = 3)
  expect_error(
    gaussian_process_lpdf(y, x = 1:2),
    "coordinate rows"
  )
})
