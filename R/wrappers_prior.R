# Internal helper to merge priors without removing NULL values
.merge_prior <- function(default, user) {
  if (is.null(user)) return(default)
  res <- user
  for (name in names(default)) {
    if (!(name %in% names(user))) {
      res[[name]] <- default[[name]]
    }
  }
  return(res)
}


#' Specify a flat prior
#'
#' @description
#' `prior_flat()` specifies that no additional prior density is added by the
#' wrapper. This is the only prior type allowed for `classic()`.
#'
#' @return A list with class `"rtmb_prior"` and `type = "flat"`.
#' @export
prior_flat <- function() {
  res <- list(type = "flat")
  class(res) <- "rtmb_prior"
  res
}


#' Specify normal/exponential priors for MAP and Bayesian inference
#'
#' @description
#' `prior_normal()` specifies normal priors for location parameters and
#' exponential priors for scale parameters. It is intended for MAP and Bayesian
#' inference, not for `classic()`.
#'
#' @param Intercept_sd Standard deviation for the intercept prior. If `NULL`, no intercept prior is added.
#' @param b_sd Standard deviation for coefficient priors. If `NULL`, no coefficient prior is added.
#' @param mu_sd Standard deviation for mean/intercept priors. If `NULL`, no mean prior is added.
#' @param sigma_rate Rate for residual standard deviation priors. If `NULL`, no sigma prior is added.
#' @param tau_rate Rate for random-effect standard deviation priors. If `NULL`, no tau prior is added.
#' @param ... Optional wrapper-specific hyperparameters.
#'
#' @return A list with class `"rtmb_prior"` and `type = "normal"`.
#' @export
prior_normal <- function(
  Intercept_sd = 10,
  mu_sd = 10,
  b_sd = 10,
  sigma_rate = 5,
  tau_rate = 5,
  ...
) {
  res <- list(
    type = "normal",
    Intercept_sd = Intercept_sd,
    mu_sd = mu_sd,
    b_sd = b_sd,
    sigma_rate = sigma_rate,
    tau_rate = tau_rate,
    ...
  )
  class(res) <- "rtmb_prior"
  res
}


#' Specify a flat prior
#'
#' @description
#' Deprecated alias of `prior_flat()`.
#'
#' @return A list with class `"rtmb_prior"` and `type = "flat"`.
#' @export
prior_uniform <- function() {
  prior_flat()
}

#' Specify a weakly informative prior
#'
#' @param sd_ratio Ratio of the prior standard deviation to the half-range of the response variable. Default is 0.5.
#' @param max_beta Maximum expected effect size. Default is 1.0.
#' @param ... Optional hyperparameters for other distributions.
#' @return A list with class "rtmb_prior"
#' @export
prior_weak <- function(sd_ratio = 0.5, max_beta = 1.0, ...) {
  res <- list(type = "weak", sd_ratio = sd_ratio, max_beta = max_beta, ...)
  class(res) <- "rtmb_prior"
  return(res)
}

#' Specify a JZS (Jeffrey-Zellner-Siow) prior for t-tests
#'
#' @param r Numeric; scale parameter for the Cauchy prior on effect size delta. Default is 0.707.
#' @return A list with class "rtmb_prior"
#' @export
prior_jzs <- function(r = 0.707) {
  res <- list(type = "jzs", r = r)
  class(res) <- "rtmb_prior"
  return(res)
}

#' Specify a Spike-and-Slab prior for variable selection
#'
#' @param ssp_ratio Prior probability of inclusion for each variable. Default is 0.25.
#' @param max_beta Maximum expected effect size. Default is 1.0.
#' @param ... Optional hyperparameters for other distributions.
#' @return A list with class "rtmb_prior"
#' @export
prior_ssp <- function(ssp_ratio = 0.25, max_beta = 1.0, ...) {
  res <- list(type = "ssp", ssp_ratio = ssp_ratio, max_beta = max_beta, ...)
  class(res) <- "rtmb_prior"
  return(res)
}

#' Specify a Regularized Horseshoe prior for continuous shrinkage
#'
#' @param expected_vars Expected number of non-zero variables. Default is 3.
#' @param slab_scale Scale parameter for the slab distribution. Default is 2.0.
#' @param slab_df Degrees of freedom for the slab distribution. Default is 4.0.
#' @param ... Optional hyperparameters for other distributions.
#' @return A list with class "rtmb_prior"
#' @export
prior_rhs <- function(expected_vars = 3, slab_scale = 2.0, slab_df = 4.0, ...) {
  res <- list(type = "rhs", expected_vars = expected_vars, slab_scale = slab_scale, slab_df = slab_df, ...)
  class(res) <- "rtmb_prior"
  return(res)
}
