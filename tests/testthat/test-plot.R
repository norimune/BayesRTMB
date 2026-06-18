test_that("plot_mdu compacts default Item labels", {
  delta <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  theta <- matrix(c(0.2, 0.1, 0.8, 0.9, 0.4, 0.6), ncol = 2, byrow = TRUE)
  rownames(delta) <- c("Item1", "Item2")

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  res <- plot_mdu(
    delta, theta,
    show_density = FALSE,
    show_radius = FALSE,
    main = ""
  )

  expect_equal(res$item_labels, c("1", "2"))
})

test_that("plot_mdu prefers rotated MDU coordinates from fitted objects", {
  delta <- matrix(c(0, 0, 1, 0), ncol = 2, byrow = TRUE)
  theta <- matrix(c(0.2, 0, 0.8, 0), ncol = 2, byrow = TRUE)
  delta_rot <- matrix(c(0, 1, 0, 2), ncol = 2, byrow = TRUE)
  theta_rot <- matrix(c(0, 1.2, 0, 1.8), ncol = 2, byrow = TRUE)

  fit <- list(
    extra = list(distance = "squared"),
    estimate = function(...) {
      list(
        delta = delta,
        theta = theta,
        delta_rot = delta_rot,
        theta_rot = theta_rot
      )
    }
  )

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  res_rot <- plot_mdu(fit, show_density = FALSE, show_radius = FALSE, main = "")
  res_raw <- plot_mdu(
    fit,
    prefer_rotated = FALSE,
    show_density = FALSE,
    show_radius = FALSE,
    main = ""
  )

  expect_equal(unname(res_rot$delta), delta_rot)
  expect_equal(unname(res_rot$theta), theta_rot)
  expect_equal(unname(res_raw$delta), delta)
  expect_equal(unname(res_raw$theta), theta)
})

test_that("plot_mdu reads MDU distance from fitted model metadata", {
  delta <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  theta <- matrix(c(0.2, 0.1, 0.8, 0.9, 0.4, 0.6), ncol = 2, byrow = TRUE)

  fit <- list(
    model = list(extra = list(distance = "euclidean")),
    estimate = function(...) {
      list(delta = delta, theta = theta, alpha = c(1, 4))
    }
  )

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot_mdu(fit, show_density = FALSE, main = ""))
})

test_that("plot_mdu supports deprecated show_phi alias", {
  delta <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  theta <- matrix(c(0.2, 0.1, 0.8, 0.9, 0.4, 0.6), ncol = 2, byrow = TRUE)

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_warning(
    plot_mdu(delta, theta, show_density = FALSE, show_phi = FALSE, main = ""),
    "deprecated"
  )
})
