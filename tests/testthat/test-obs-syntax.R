test_that("obs() lhs with one value is equivalent to shorthand sampling syntax", {
  dat <- list(
    y = c(-0.2, 0.1, 0.4),
    mu = c(0, 0.2, 0.3)
  )
  par <- list(sigma = 1.2)

  fn_obs <- model_code({
    obs(y) ~ normal(mu, sigma)
  })
  fn_plain <- model_code({
    y ~ normal(mu, sigma)
  })

  expect_equal(fn_obs(dat, par), fn_plain(dat, par), tolerance = 1e-12)
})

test_that("obs() lhs passes multiple observed values to user-defined lpdf", {
  paired_lpdf <- function(x, z, theta, sum = TRUE) {
    res <- -0.5 * (x + z - theta)^2
    if (sum) sum(res) else res
  }

  dat <- list(
    y = c(0.2, -0.3, 0.5),
    z = c(1.0, 0.5, -0.2)
  )
  par <- list(theta = 0.4)

  fn_obs <- model_code({
    obs(y, z) ~ paired(theta)
  })
  fn_old <- model_code({
    y ~ paired(z, theta)
  })

  expect_equal(fn_obs(dat, par), paired_lpdf(dat$y, dat$z, par$theta), tolerance = 1e-12)
  expect_equal(fn_obs(dat, par), fn_old(dat, par), tolerance = 1e-12)
})

test_that("obs() lhs works inside rtmb_code and get_n_obs", {
  dat <- list(
    y = c(0.2, -0.3, 0.5),
    z = c(1.0, 0.5, -0.2)
  )

  code <- rtmb_code(
    setup = {
      paired_lpdf <- function(x, z, theta, sum = TRUE) {
        res <- -0.5 * (x + z - theta)^2
        if (sum) sum(res) else res
      }
    },
    parameters = {
      theta <- Dim()
    },
    model = {
      obs(y, z) ~ paired(theta)
      theta ~ normal(0, 1)
    }
  )

  mdl <- rtmb_model(dat, code, init = list(theta = 0), silent = TRUE)
  expect_s3_class(mdl, "RTMB_Model")
  expect_equal(mdl$get_n_obs(), length(dat$y))
})

test_that("empty obs() lhs gives a clear error", {
  expect_error(
    model_code({
      obs() ~ normal(0, 1)
    }),
    "obs\\(\\) on the left side"
  )
})

test_that("rtmb_code captures external setup helpers used by sampling syntax", {
  paired_external_lpdf <- function(x, z, theta, sum = TRUE) {
    res <- -0.5 * (x + z - theta)^2
    if (sum) sum(res) else res
  }

  dat <- list(
    y = c(0.2, -0.3, 0.5),
    z = c(1.0, 0.5, -0.2)
  )

  code <- rtmb_code(
    setup = {
      paired_lpdf <- paired_external_lpdf
    },
    parameters = {
      theta <- Dim()
    },
    model = {
      obs(y, z) ~ paired(theta)
      theta ~ normal(0, 1)
    }
  )

  mdl <- rtmb_model(dat, code, init = list(theta = 0), silent = TRUE)
  expect_s3_class(mdl, "RTMB_Model")
})

test_that("rtmb_code captures recursive dependencies of external setup functions", {
  helper_scale <- 0.5
  helper_shift <- function(x, z) helper_scale * (x + z)
  paired_external_lpdf <- function(x, z, theta, sum = TRUE) {
    res <- -0.5 * (helper_shift(x, z) - theta)^2
    if (sum) sum(res) else res
  }

  dat <- list(
    y = c(0.2, -0.3, 0.5),
    z = c(1.0, 0.5, -0.2)
  )

  code <- rtmb_code(
    setup = {
      paired_lpdf <- paired_external_lpdf
    },
    parameters = {
      theta <- Dim()
    },
    model = {
      obs(y, z) ~ paired(theta)
      theta ~ normal(0, 1)
    }
  )

  rm(helper_scale, helper_shift, paired_external_lpdf)

  mdl <- rtmb_model(dat, code, init = list(theta = 0), silent = TRUE)
  expect_s3_class(mdl, "RTMB_Model")
})

test_that("setup-defined functions retain external helpers after setup cleanup", {
  setup_helper_scale <- 0.5
  setup_helper <- function(x, z) setup_helper_scale * (x + z)

  dat <- list(
    y = c(0.2, -0.3, 0.5),
    z = c(1.0, 0.5, -0.2)
  )

  code <- rtmb_code(
    setup = {
      paired_lpdf <- function(x, z, theta, sum = TRUE) {
        res <- -0.5 * (setup_helper(x, z) - theta)^2
        if (sum) sum(res) else res
      }
    },
    parameters = {
      theta <- Dim()
    },
    model = {
      obs(y, z) ~ paired(theta)
      theta ~ normal(0, 1)
    }
  )

  rm(setup_helper_scale, setup_helper)

  mdl <- rtmb_model(dat, code, init = list(theta = 0), silent = TRUE)
  expect_s3_class(mdl, "RTMB_Model")
})

test_that("data-supplied lpdf functions can use package lpmf helpers", {
  ez_lpdf <- function(correct, n_acc, p, sum = TRUE) {
    binomial_lpmf(correct, n_acc, p, sum = sum)
  }

  dat <- list(
    correct = 3,
    n_acc = 5,
    ez_lpdf = ez_lpdf
  )

  code <- rtmb_code(
    parameters = {
      eta <- Dim()
    },
    transform = {
      p <- inv_logit(eta)
    },
    model = {
      obs(correct, n_acc) ~ ez(p)
      eta ~ normal(0, 1)
    }
  )

  mdl <- rtmb_model(dat, code, init = list(eta = 0), silent = TRUE)
  expect_s3_class(mdl, "RTMB_Model")
})
