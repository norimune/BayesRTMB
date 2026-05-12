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


#' Specify a uniform or manual prior
#'
#' @param Intercept_sd Standard deviation for the intercept prior (Normal). Default is NULL (flat).
#' @param b_sd Standard deviation for the coefficients prior (Normal). Default is NULL (flat).
#' @param sigma_rate Rate for the residual standard deviation prior (Exponential). Default is NULL (flat).
#' @param tau_rate Rate for the random effects standard deviation prior (Exponential). Default is NULL (flat).
#' @param ... Optional hyperparameters
#' @return A list with class "rtmb_prior"
#' @export
prior_uniform <- function(Intercept_sd = NULL, b_sd = NULL, mu_sd = NULL, sigma_rate = NULL, tau_rate = NULL, ...) {
  res <- list(type = "uniform", Intercept_sd = Intercept_sd, b_sd = b_sd, mu_sd = mu_sd, sigma_rate = sigma_rate, tau_rate = tau_rate, ...)
  class(res) <- "rtmb_prior"
  return(res)
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
