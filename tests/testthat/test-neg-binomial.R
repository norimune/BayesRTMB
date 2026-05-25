test_that("negative binomial 2 log pmf matches stats parameterization", {
  x <- c(0, 1, 3, 8)
  mu <- c(0.5, 1.2, 4.0, 9.5)
  size <- 2.75

  expect_equal(
    neg_binomial_2_lpmf(x, mu, size, sum = FALSE),
    stats::dnbinom(x, size = size, mu = mu, log = TRUE),
    tolerance = 1e-12
  )
})

test_that("negative binomial GLM is AD compatible", {
  set.seed(1)
  dat <- data.frame(
    HR = rnbinom(20, mu = 10, size = 2),
    Steal = rnorm(20),
    Walk = rnorm(20),
    K = rnorm(20)
  )

  mdl <- rtmb_glm(
    HR ~ Steal + Walk + K,
    data = dat,
    family = "neg_binomial",
    prior = prior_normal()
  )
  obj <- mdl$get_ad_obj()

  expect_true(all(is.finite(obj$gr(obj$par))))
})
