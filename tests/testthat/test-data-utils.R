test_that("to_long preserves range column order before sorting", {
  dat <- data.frame(
    id = 1:2,
    setNames(as.data.frame(matrix(seq_len(24), nrow = 2)), paste0("X", 1:12))
  )

  res <- to_long(dat, within = "X1:X12", id = "id")

  expect_identical(levels(res$Condition), paste0("X", 1:12))
  expect_identical(as.character(res$Condition[res$id == 1]), paste0("X", 1:12))
})

test_that("to_long preserves input row order by default", {
  dat <- data.frame(
    id = c(2L, 1L),
    X1 = c(20L, 10L),
    X2 = c(21L, 11L)
  )

  res <- to_long(dat, within = "X1:X2", id = "id")

  expect_identical(res$id, c(2L, 2L, 1L, 1L))
  expect_identical(as.character(res$Condition), c("X1", "X2", "X1", "X2"))
  expect_identical(res$Value, c(20L, 21L, 10L, 11L))
})

test_that("to_long can sort by id and label when requested", {
  dat <- data.frame(
    id = c(2L, 1L),
    X1 = c(20L, 10L),
    X2 = c(21L, 11L)
  )

  res <- to_long(dat, within = "X1:X2", id = "id", sort = TRUE)

  expect_identical(res$id, c(1L, 1L, 2L, 2L))
  expect_identical(as.character(res$Condition), c("X1", "X2", "X1", "X2"))
  expect_identical(res$Value, c(10L, 11L, 20L, 21L))
})

test_that("to_long handles multiple within groups", {
  dat <- data.frame(
    id = 1:2,
    X1 = 1:2,
    X2 = 3:4,
    X3 = 5:6,
    Y1 = 11:12,
    Y2 = 13:14,
    Y3 = 15:16
  )

  res <- to_long(
    dat,
    within = c("X1:X3", "Y1:Y3"),
    label = "time",
    value = c("x", "y"),
    id = "id"
  )

  expect_identical(names(res), c("id", "time", "x", "y"))
  expect_identical(levels(res$time), c("X1", "X2", "X3"))
  expect_identical(as.character(res$time[res$id == 1]), c("X1", "X2", "X3"))
  expect_identical(res$x[res$id == 1], c(1L, 3L, 5L))
  expect_identical(res$y[res$id == 1], c(11L, 13L, 15L))
})

test_that("to_long concatenates list of column name vectors for one value", {
  dat <- data.frame(
    id = 1:2,
    v1 = 1:2,
    v2 = 3:4,
    w1 = 11:12,
    w2 = 13:14
  )
  A_name <- c("v1", "v2")
  B_name <- c("w1", "w2")

  res <- to_long(
    dat,
    within = list(A_name, B_name),
    label = "item",
    value = "score",
    id = "id"
  )

  expect_identical(names(res), c("id", "item", "score"))
  expect_identical(levels(res$item), c("v1", "v2", "w1", "w2"))
  expect_identical(as.character(res$item[res$id == 1]), c("v1", "v2", "w1", "w2"))
  expect_identical(res$score[res$id == 1], c(1L, 3L, 11L, 13L))
})

test_that("to_long warns when multiple within groups have different lengths", {
  dat <- data.frame(
    id = 1:2,
    X1 = 1:2,
    X2 = 3:4,
    Y1 = 11:12
  )

  expect_warning(
    expect_error(
      to_long(dat, within = c("X1:X2", "Y1"), value = c("x", "y"), id = "id"),
      "same number of columns"
    ),
    "same number of columns"
  )
})

test_that("to_wide preserves single-value column names", {
  dat <- data.frame(
    id = c(1L, 1L, 2L, 2L),
    time = c("X1", "X2", "X1", "X2"),
    score = c(10L, 20L, 30L, 40L)
  )

  res <- to_wide(dat, within = "time", value = "score", id = "id")

  expect_identical(names(res), c("id", "X1", "X2"))
  expect_identical(res$X1, c(10L, 30L))
  expect_identical(res$X2, c(20L, 40L))
})

test_that("to_wide handles multiple value columns", {
  dat <- data.frame(
    id = c(1L, 1L, 2L, 2L),
    time = c("X1", "X2", "X1", "X2"),
    x = c(10L, 20L, 30L, 40L),
    y = c(11L, 21L, 31L, 41L)
  )

  res <- to_wide(dat, within = "time", value = c("x", "y"), id = "id")

  expect_identical(names(res), c("id", "x.X1", "y.X1", "x.X2", "y.X2"))
  expect_identical(res$x.X1, c(10L, 30L))
  expect_identical(res$x.X2, c(20L, 40L))
  expect_identical(res$y.X1, c(11L, 31L))
  expect_identical(res$y.X2, c(21L, 41L))
})
