#' RTMB-based Linear Mixed Model (LMM) wrapper function
#'
#' @param formula Formula
#' @param data Data frame
#' @param laplace Logical; whether to marginalize random effects using Laplace approximation
#' @param prior An object of class `"rtmb_prior"`.
#' Use `prior_flat()` for no prior, `prior_normal()` for default normal/exponential priors,
#' or `prior_weak()`, `prior_rhs()`, `prior_ssp()` for weakly informative or regularized Bayesian inference.
#' Default is `prior_flat()`.
#' @param y_range Theoretical minimum and maximum values of the response variable
#' @param init Initial values

#' @param gmc Character vector of variable names for GMC
#' @param centering Alias for `gmc`.
#' @param cwc List for CWC
#' @param view Character vector of parameter names to prioritize in summary.
#' @param factors Character vector of variable names to be treated as factors.
#' @param contrasts Character string specifying the contrast type ("treatment" or "sum").
#' @param sigma_by Character vector specifying variables to group residual variance by (heteroscedasticity).
#' @param resid_corr Residual correlation structure (e.g., "ar1", "cs", "un", "toep").
#' @param resid_time Variable name for time points in residual correlation.
#' @param resid_group Variable name for grouping in residual correlation.
#' @param within Optional list for wide-to-long conversion. For repeated measures data in wide format,
#' specify the factor names and their levels, e.g., \code{list(Time = 4)} or \code{list(A = 2, B = 3)}.
#' The total number of levels must match the number of columns in \code{cbind()} on the LHS.
#' If omitted and the LHS is \code{cbind()}, the within-factor name is inferred from RHS variables not present in the data.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param missing Missing value handling strategy: "listwise".
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Additional arguments passed to \code{rtmb_model()}.
#'
#' @return RTMB_Model object
#' @export
#' @example inst/examples/ex_lm.R
rtmb_lmer <- function(formula, data, laplace = TRUE,
                       prior = prior_flat(),
                       y_range = NULL,
                       init = NULL,
                       fixed = NULL,
                       gmc = NULL,
                       centering = NULL,
                       cwc = NULL,
                       view = NULL,
                       sigma_by = NULL,
                       factors = NULL,
                       contrasts = "treatment",
                       resid_corr = NULL,
                       resid_time = NULL,
                       resid_group = NULL,
                       within = NULL,
                       missing = c("listwise", "fiml"),
                       WAIC = FALSE,
                       ...) {
  missing <- match.arg(missing)
  rtmb_glmer(formula = formula, data = data, family = "gaussian",
             laplace = laplace,
             prior = prior,
             y_range = y_range,
             init = init,
             fixed = fixed,
             gmc = gmc,
             centering = centering,
             cwc = cwc,
             view = view,
             sigma_by = sigma_by,
             factors = factors,
             contrasts = contrasts,
             resid_corr = resid_corr,
             resid_time = resid_time,
             resid_group = resid_group,
             within = within,
             missing = missing,
             WAIC = WAIC,
             .force_sum = TRUE)
}

#' RTMB-based GLM wrapper function (no random effects)
#'
#' @param formula Formula
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param prior An object of class `"rtmb_prior"`.
#' Use `prior_flat()` for no prior, `prior_normal()` for default normal/exponential priors,
#' or `prior_weak()`, `prior_rhs()`, `prior_ssp()` for weakly informative or regularized Bayesian inference.
#' Default is `prior_flat()`.
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param init List of initial values

#' @param gmc Character vector of variable names for GMC
#' @param centering Alias for `gmc`.
#' @param factors Character vector of variable names to be treated as factors.
#' @param contrasts Character string specifying the contrast type ("treatment" or "sum").
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param view Optional character vector of parameter names to show first in summary.
#' @param missing Missing value handling strategy: "listwise".
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Additional arguments passed to \code{rtmb_glmer()}.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glm <- function(formula, data, family = "gaussian",
                       prior = prior_flat(),
                       y_range = NULL,
                       init = NULL, fixed = NULL,
                       gmc = NULL,
                       centering = NULL,
                       view = NULL,
                       factors = NULL,
                       contrasts = "treatment",
                       missing = c("listwise", "fiml"), WAIC = FALSE, ...) {
  missing <- match.arg(missing)
  rtmb_glmer(formula = formula, data = data, family = family,
             laplace = FALSE,
             prior = prior,
             y_range = y_range,
             init = init,
             fixed = fixed,
             gmc = gmc,
             centering = centering,
             view = view,
             factors = factors,
             contrasts = contrasts, missing = missing, WAIC = WAIC, ...)
}

#' RTMB-based Linear Regression wrapper function
#'
#' @param formula Formula (e.g., Y ~ X1 + X2)
#' @param data Data frame
#' @param prior An object of class `"rtmb_prior"`.
#' Use `prior_flat()` for no prior, `prior_normal()` for default normal/exponential priors,
#' or `prior_weak()`, `prior_rhs()`, `prior_ssp()` for weakly informative or regularized Bayesian inference.
#' Default is `prior_flat()`.
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param init List of initial values.

#' @param gmc Character vector of variable names for GMC
#' @param centering Alias for `gmc`.
#' @param factors Character vector of variable names to be treated as factors.
#' @param contrasts Character string specifying the contrast type ("treatment" or "sum").
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param view Optional character vector of parameter names to show first in summary.
#' @param missing Missing value handling strategy: "listwise".
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Additional arguments passed to \code{rtmb_glmer()}.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_lm <- function(formula, data,
                    prior = prior_flat(),
                    y_range = NULL,
                    init = NULL, fixed = NULL,
                    gmc = NULL,
                    centering = NULL,
                    view = NULL,
                    factors = NULL,
                    contrasts = "treatment",
                    missing = c("listwise", "fiml"), WAIC = FALSE, ...) {
  missing <- match.arg(missing)
  rtmb_glmer(formula = formula, data = data, family = "gaussian",
             laplace = FALSE,
             prior = prior,
             y_range = y_range,
             init = init,
             fixed = fixed,
             gmc = gmc,
             centering = centering,
             view = view,
             factors = factors,
             contrasts = contrasts, missing = missing, WAIC = WAIC, ...)
}
