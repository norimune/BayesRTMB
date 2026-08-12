test_that("centering helpers retain vector shape and missing values", {
  x <- c(1, 3, NA, 2, 6)

  grand <- center_grand_mean(x)
  expect_equal(grand, x - mean(x, na.rm = TRUE))
  expect_null(dim(grand))
  expect_true(is.na(grand[3]))

  cluster <- c("a", "a", "a", "b", "b")
  within <- center_within_cluster(x, cluster)
  expect_equal(within[c(1, 2)], c(-1, 1))
  expect_true(is.na(within[3]))
  expect_equal(within[c(4, 5)], c(-2, 2))
  expect_null(dim(within))
})

test_that("centering helpers validate their inputs", {
  expect_error(center_grand_mean(factor(c("a", "b"))), "numeric vector")
  expect_error(center_grand_mean(matrix(1:4, 2)), "numeric vector")
  expect_error(
    center_within_cluster(1:3, c("a", "b")),
    "same length"
  )

  centered <- center_within_cluster(c(1, 2, 3), c("a", NA, "a"))
  expect_equal(centered[c(1, 3)], c(-1, 1))
  expect_true(is.na(centered[2]))
})

test_that("wrapper setup displays explicit GMC and CWC operations", {
  dat <- data.frame(
    y = c(1, 2, 3, 4, 5, 6),
    x = c(1, 2, 3, 10, 11, 12),
    z = c(2, 4, 6, 20, 22, 24),
    id = factor(rep(c("a", "b"), each = 3))
  )

  mdl <- rtmb_glmer(
    y ~ x + z + (1 | id),
    data = dat,
    gmc = "x",
    cwc = list(id, "z")
  )
  out <- capture.output(mdl$print_code())

  expect_true(any(grepl("x <- center_grand_mean\\(x\\)", out)))
  expect_true(any(grepl("z <- center_within_cluster\\(z, id\\)", out)))
  expect_equal(mean(mdl$raw_data$x), 0, tolerance = 1e-12)
  expect_equal(
    as.numeric(tapply(mdl$raw_data$z, mdl$raw_data$id, mean)),
    c(0, 0),
    tolerance = 1e-12
  )
})
