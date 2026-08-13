test_that(".data exposes the original data object only during setup", {
  dat <- data.frame(
    a = c(1, 2, 3, 4),
    b = c(2, 3, 4, 5)
  )

  code <- rtmb_code(
    setup = {
      Y <- as.matrix(.data)
      column_sum <- a + b
      N <- nrow(Y)
      P <- ncol(Y)
    },
    parameters = {
      mu <- Dim(P)
    },
    model = {
      for (j in 1:P) {
        Y[, j] ~ normal(mu[j], 1)
      }
    }
  )

  model <- rtmb_model(data = dat, code = code, silent = TRUE)

  expect_equal(model$data$Y, as.matrix(dat))
  expect_equal(model$data$column_sum, dat$a + dat$b)
  expect_false(".data" %in% names(model$data))
})

test_that(".data allows a matrix to be passed directly to rtmb_model", {
  dat <- matrix(seq_len(8), nrow = 4, ncol = 2)

  code <- rtmb_code(
    setup = {
      Y <- as.matrix(.data)
      N <- nrow(Y)
      P <- ncol(Y)
    },
    parameters = {
      mu <- Dim(P)
    },
    model = {
      for (j in 1:P) {
        Y[, j] ~ normal(mu[j], 1)
      }
    }
  )

  model <- rtmb_model(data = dat, code = code, silent = TRUE)

  expect_equal(model$data$Y, dat)
  expect_false(".data" %in% names(model$data))
})

test_that(".data is a reserved read-only setup binding", {
  code <- rtmb_code(
    setup = {
      N <- 1
    },
    parameters = {
      mu <- Dim()
    },
    model = {
      mu ~ normal(0, 1)
    }
  )

  expect_error(
    rtmb_model(data = list(.data = 1), code = code, silent = TRUE),
    "reserved"
  )

  assign_code <- rtmb_code(
    setup = {
      .data <- NULL
      N <- 1
    },
    parameters = {
      mu <- Dim()
    },
    model = {
      mu ~ normal(0, 1)
    }
  )

  expect_error(
    rtmb_model(data = list(x = 1), code = assign_code, silent = TRUE),
    "\\.data"
  )
})
