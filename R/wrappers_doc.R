#' Common Features and Arguments of RTMB Wrapper Functions
#'
#' @description
#' The RTMB wrapper functions (`rtmb_lm`, `rtmb_glm`, `rtmb_glmer`, `rtmb_fa`, etc.)
#' share a unified interface designed to make Bayesian and Frequentist inference
#' accessible through familiar R formulas and standard model specifications.
#'
#' @details
#' All wrapper functions in this package are built upon the same core engine.
#' This ensures that regardless of the model type, the workflow for estimation,
#' summary, and expansion remains consistent.
#'
#' \strong{1. Unified Inference Methods:}
#' Every model object returned by a wrapper function provides the following methods:
#' \itemize{
#'   \item \code{$optimize()}: Performs MAP estimation (comparable to MLE).
#'   \item \code{$sample()}: Performs MCMC sampling (NUTS) for full Bayesian inference.
#'   \item \code{$variational()}: Performs Variational Inference (ADVI) for fast posterior approximation.
#' }
#'
#' \strong{2. Prior API and Regularization:}
#' You can specify the prior distribution using the `prior` argument:
#' \itemize{
#'   \item \code{prior_flat()}: No additional prior (default). Required for \code{classic()} inference.
#'   \item \code{prior_normal()}: Manual normal/exponential priors for parameters.
#'   \item \code{prior_weak()}: Weakly informative priors based on \code{y_range}.
#'   \item \code{prior_rhs()}: Regularized Horseshoe prior for continuous shrinkage.
#'   \item \code{prior_ssp()}: Spike-and-Slab prior for sparse variable selection.
#' }
#' \emph{Note: When using `prior_weak()`, `prior_rhs()`, or `prior_ssp()`, you must specify \code{y_range = c(min, max)}
#' to let the model set appropriate global scales for the priors.}
#'
#' \strong{3. Weakly Informative Priors (\code{y_range}):}
#' By providing the theoretical range of your response variable via \code{y_range}
#' while keeping the default \code{prior = prior_flat()},
#' the wrappers automatically switch to \code{prior_weak()}. These priors
#' are designed to be broad enough to cover any reasonable value but narrow enough
#' to stabilize the estimation and prevent the sampler from wandering into
#' non-sensical parameter space.
#'
#' \strong{4. Fixed vs. Random Effects:}
#' For mixed-effect models (e.g., \code{rtmb_glmer}), random effects are marginalized
#' using the Laplace approximation during \code{$optimize(laplace = TRUE)}.
#' When using \code{$sample()}, random effects are treated as unknown parameters
#' and sampled alongside fixed effects.
#'
#' \strong{5. Null Model Creation (\code{null}):}
#' You can specify a \code{null} argument (e.g., \code{null = "x1 ~ normal(0, 0.1)"})
#' in the wrappers to simultaneously create a restricted version of your model.
#' This is particularly useful for computing Bayes Factors or performing
#' model comparisons.
#'
#' @name rtmb_wrappers
#' @family wrappers
NULL
