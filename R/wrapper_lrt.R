#' Fit a Latent Rank Theory (LRT) Model
#'
#' @description
#' Fits a Latent Rank Theory model, which is a mixture model with ordered ranks
#' and Gaussian Process smoothing on the mean profiles.
#'
#' @param formula A formula specifying the response variable(s).
#' @param k Number of ranks (mixture components).
#' @param rank_coords Optional numeric vector of coordinates for each rank. Default is 1:k.
#' @param covariance Covariance structure: "diagonal", "diagonal_equal", "full", "full_equal", or "full_equal_corr".
#' @param data A data frame containing the variables.
#' @param magnitude Signal standard deviation for the GP prior. If NULL, it is estimated.
#' @param smoothing Length-scale for the GP prior. If NULL, it is estimated.
#' @param noise Measurement noise for the GP prior (default is 0.01).
#' @param prob_smoothing Logical; whether to apply smoothing to the class membership probabilities.
#' @param link Link function for class probabilities: "ordered" or "sequential".
#' @param prior Prior configuration: `prior_flat()`, `prior_normal()`,
#'   `prior_weak()`, `prior_rhs()`, or `prior_ssp()`. Default is
#'   `prior_flat()`. If `y_range` is supplied with the default flat prior,
#'   the wrapper automatically switches to `prior_weak()`.
#' @param y_range Optional numeric vector or matrix defining the theoretical range (min, max) of response variables.
#'   Specifying this automatically enables weakly informative priors if `prior` is `prior_flat()`.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param two_stage Logical; if TRUE, estimate the latent-rank measurement model
#'   first and then estimate the rank regression with delta-method uncertainty
#'   propagation. Currently supported for `$optimize()` only.
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Additional arguments passed to `rtmb_model`.
#' @return A \code{RTMB_Model} object.
#' @export
rtmb_lrt <- function(formula, k = 3, data = NULL,
                     rank_coords = NULL,
                     covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
                     magnitude = NULL, smoothing = NULL, noise = 0.01,
                     prob_smoothing = FALSE,
                     link = c("ordered", "sequential"),
                     prior = prior_flat(), y_range = NULL, fixed = NULL,
                     two_stage = FALSE, WAIC = FALSE, ...) {

  dots <- list(...)
  if ("2stage" %in% names(dots)) {
    two_stage <- isTRUE(dots[["2stage"]])
  }
  wrapper_optimize_defaults <- list()
  if (!is.null(dots$num_estimate)) {
    wrapper_optimize_defaults$num_estimate <- dots$num_estimate
  }

  if (is.null(rank_coords)) rank_coords <- 1:k

  if (is.null(prior)) {
    prior <- prior_flat()
  }

  # Automatically switch to prior_weak() if y_range is provided and prior is default flat
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && identical(prior$type, "flat")) {
    prior <- prior_weak()
  }

  prior <- .validate_prior_type(
    prior,
    allowed = c("flat", "normal", "weak", "rhs", "ssp"),
    context = "rtmb_lrt()"
  )

  # NSE for formula: handle case where formula is just a variable name in data
  formula_expr <- substitute(formula)
  formula_val <- try(formula, silent = TRUE)

  if (inherits(formula_val, "try-error") ||
      (!inherits(formula_val, "formula") && !is.matrix(formula_val) && !is.vector(formula_val))) {
    if (!is.null(data) && all(all.vars(formula_expr) %in% names(data))) {
      formula <- as.formula(call("~", formula_expr, 1))
    } else if (!inherits(formula_val, "try-error")) {
      formula <- formula_val
    } else {
      stop(paste0("Object '", deparse(formula_expr), "' not found."))
    }
  } else {
    formula <- formula_val
  }

  covariance <- match.arg(covariance)
  link <- match.arg(link)
  K_mix <- k
  prior_type <- prior$type

  if (prior_type == "flat") {
    # no merge
  } else if (prior_type == "normal") {
    normal_defaults <- list(
      Intercept_sd = NULL,
      mu_sd = NULL,
      b_sd = NULL,
      sigma_rate = NULL,
      lkj_eta = NULL,
      mag_rate = NULL,
      smooth_rate = NULL,
      cutpoint_sd = NULL
    )
    prior <- .merge_prior(normal_defaults, prior)
  } else {
    default_prior <- list(Intercept_sd = 10, mu_sd = 10, b_sd = 10, sigma_rate = 1, lkj_eta = 1.0, mag_rate = 1.0, smooth_rate = 1.0, cutpoint_sd = 2.5)
    prior <- .merge_prior(default_prior, prior)
  }

  # Sync mu_sd and Intercept_sd
  if (prior_type == "normal") {
    if (is.null(prior$mu_sd) && !is.null(prior$Intercept_sd)) {
      prior$mu_sd <- prior$Intercept_sd
    }
  }

  regularization <- if (prior_type %in% c("rhs", "ssp")) prior_type else "none"
  use_weak_info <- prior_type %in% c("weak", "rhs", "ssp")

  if (use_weak_info) {
    if (is.null(prior$max_beta)) prior$max_beta <- 1.0
    if (is.null(prior$expected_vars)) prior$expected_vars <- 3
    if (is.null(prior$slab_scale)) prior$slab_scale <- 2.0
    if (is.null(prior$slab_df)) prior$slab_df <- 4.0
    if (is.null(prior$ssp_ratio)) prior$ssp_ratio <- 0.25
  }

  setup_from_formula <- inherits(formula, "formula") && !is.null(data)

  # Data parsing
  if (is.matrix(formula) || is.vector(formula)) {
    Y_mat <- as.matrix(formula)
    X_prob <- NULL
  } else if (!inherits(formula, "formula")) {
    var_name <- as.character(formula)
    if (!is.null(data) && var_name %in% names(data)) {
      Y_mat <- as.matrix(data[[var_name]])
    } else {
      Y_mat <- as.matrix(formula)
    }
    X_prob <- NULL
  } else {
    mf_y <- model.frame(formula, data = data)
    Y_mat <- model.response(mf_y)
    if (is.null(Y_mat)) Y_mat <- as.matrix(mf_y)

    # Check if formula has RHS covariates for prob
    rhs <- formula[[3]]
    is_empty_rhs <- (is.numeric(rhs) && rhs == 1) || (is.symbol(rhs) && rhs == "1")

    if (!is_empty_rhs) {
      X_prob <- model.matrix(formula, data = data)
      if ("(Intercept)" %in% colnames(X_prob)) {
        X_prob <- X_prob[, colnames(X_prob) != "(Intercept)", drop = FALSE]
      }
      if (ncol(X_prob) == 0) X_prob <- NULL
    } else {
      X_prob <- NULL
    }
  }

  if (setup_from_formula) {
    setup_vars <- all.vars(formula)
    missing_vars <- setdiff(setup_vars, names(data))
    if (length(missing_vars) > 0) {
      stop(
        "The following variables in formula are not found in data: ",
        paste(missing_vars, collapse = ", "),
        call. = FALSE
      )
    }
    setup_df <- as.data.frame(data)[, setup_vars, drop = FALSE]
    class(setup_df) <- c("rtmb_setup_df", class(setup_df))
  }

  Y_mat <- as.matrix(Y_mat)
  N_obs <- nrow(Y_mat)
  P_dim <- ncol(Y_mat)
  if (is.null(colnames(Y_mat))) colnames(Y_mat) <- paste0("Y", 1:P_dim)
  has_cov_prob <- !is.null(X_prob)
  K_prob <- if (has_cov_prob) ncol(X_prob) else 0

  if (isTRUE(two_stage) && !has_cov_prob) {
    stop("two_stage = TRUE requires covariates on the right-hand side of the formula.", call. = FALSE)
  }

  multivariate <- P_dim > 1
  is_diag <- covariance %in% c("diagonal", "diagonal_equal")
  is_sigma_equal <- covariance %in% c("diagonal_equal", "full_equal", "full_equal_corr")

  # --- 1. Setup ---
  if (setup_from_formula) {
    setup_exprs <- list(
      as.name("{"),
      quote(mf <- model.frame(formula, df)),
      quote(Y <- model.response(mf)),
      quote(if (is.null(Y)) Y <- as.matrix(mf)),
      quote(Y <- as.matrix(Y)),
      quote(N <- nrow(Y)),
      quote(P <- ncol(Y)),
      quote(K <- K),
      quote(rank_coords <- rank_coords)
    )
  } else {
    setup_exprs <- list(
      as.name("{"),
      quote(N <- N),
      quote(K <- K),
      quote(P <- P),
      quote(rank_coords <- rank_coords)
    )
  }
  if (has_cov_prob) {
    if (setup_from_formula) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_prob <- model.matrix(formula, mf))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(if ("(Intercept)" %in% colnames(X_prob)) X_prob <- X_prob[, colnames(X_prob) != "(Intercept)", drop = FALSE])
      setup_exprs[[length(setup_exprs) + 1]] <- quote(K_prob <- ncol(X_prob))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_means <- matrix(colMeans(X_prob), 1, K_prob))
    } else {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(K_prob <- K_prob)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_means <- matrix(X_means, 1, K_prob))
    }
    if (regularization == "rhs") {
      # RHS setup (cap expected_vars to prevent zero-division)
      p0_rhs <- min(prior$expected_vars, K_prob - 1)
      if (p0_rhs < 1) p0_rhs <- 1
      setup_exprs[[length(setup_exprs) + 1]] <- if (setup_from_formula) {
        bquote(tau0 <- .(p0_rhs) / (K_prob - .(p0_rhs)) / sqrt(N))
      } else {
        quote(tau0 <- tau0)
      }
      setup_exprs[[length(setup_exprs) + 1]] <- if (setup_from_formula) bquote(half_slab_df <- .(prior$slab_df) / 2) else quote(half_slab_df <- half_slab_df)
      setup_exprs[[length(setup_exprs) + 1]] <- if (setup_from_formula) bquote(half_slab_scale2 <- .(prior$slab_scale)^2 / 2) else quote(half_slab_scale2 <- half_slab_scale2)
    } else if (regularization == "ssp") {
      # SSP setup
      setup_exprs[[length(setup_exprs) + 1]] <- if (setup_from_formula) bquote(tau_scale <- .(prior$max_beta) / 1.96) else quote(tau_scale <- tau_scale)
    }
  }
  setup_ast <- as.call(setup_exprs)

  # --- 2. Parameters ---
  param_exprs <- list(as.name("{"))
  param_exprs[[length(param_exprs) + 1]] <- bquote(mu_p <- Dim(c(.(P_dim), .(K_mix)), type = "ordered"))

  # Sigma parameterization based on covariance
  if (is_sigma_equal) {
    param_exprs[[length(param_exprs) + 1]] <- bquote(sigma <- Dim(.(P_dim), lower = 0))
  } else {
    param_exprs[[length(param_exprs) + 1]] <- bquote(sigma <- Dim(c(.(K_mix), .(P_dim)), lower = 0))
  }

  if (multivariate && !is_diag) {
    if (covariance == "full") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr <- Dim(c(.(K_mix), .(P_dim), .(P_dim)), type = "CF_corr"))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr <- Dim(c(.(P_dim), .(P_dim)), type = "CF_corr"))
    }
  }

  if (has_cov_prob) {
    if (link == "ordered") {
      if (regularization == "rhs") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(z_b <- Dim(.(K_prob)))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_lambda_b <- Dim(.(K_prob), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_tau_b <- Dim(lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(c2_b <- Dim(lower = 0))
      } else if (regularization == "ssp") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b_raw <- Dim(.(K_prob)))
        param_exprs[[length(param_exprs) + 1]] <- bquote(r_b <- Dim(.(K_prob), lower = 0, upper = 1))
        param_exprs[[length(param_exprs) + 1]] <- bquote(tau_b <- Dim(lower = 0))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b <- Dim(.(K_prob)))
      }
      param_exprs[[length(param_exprs) + 1]] <- bquote(cutpoints <- Dim(.(K_mix - 1), type = "ordered"))
    } else {
      # Sequential: Separate b for each step
      if (regularization == "rhs") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(z_b <- Dim(c(.(K_prob), .(K_mix - 1))))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_lambda_b <- Dim(c(.(K_prob), .(K_mix - 1)), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_tau_b <- Dim(.(K_mix - 1), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(c2_b <- Dim(.(K_mix - 1), lower = 0))
      } else if (regularization == "ssp") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b_raw <- Dim(c(.(K_prob), .(K_mix - 1))))
        param_exprs[[length(param_exprs) + 1]] <- bquote(r_b <- Dim(c(.(K_prob), .(K_mix - 1)), lower = 0, upper = 1))
        param_exprs[[length(param_exprs) + 1]] <- bquote(tau_b <- Dim(.(K_mix - 1), lower = 0))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b <- Dim(c(.(K_prob), .(K_mix - 1))))
      }
      param_exprs[[length(param_exprs) + 1]] <- bquote(alpha <- Dim(.(K_mix - 1)))
    }
  } else {
    if (prob_smoothing) {
      param_exprs[[length(param_exprs) + 1]] <- bquote(logit_theta <- Dim(.(K_mix)))
      param_exprs[[length(param_exprs) + 1]] <- quote(magnitude_theta <- Dim(lower = 0))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- bquote(theta <- Dim(.(K_mix), type = "simplex"))
    }
  }

  if (is.null(magnitude)) param_exprs[[length(param_exprs) + 1]] <- quote(magnitude <- Dim(lower = 0))
  if (is.null(smoothing)) param_exprs[[length(param_exprs) + 1]] <- quote(smoothing <- Dim(lower = 0))
  param_ast <- as.call(param_exprs)

  # --- 3. Model ---
  model_exprs <- list(as.name("{"))

  # Regularization transformations
  if (has_cov_prob) {
    if (regularization == "rhs") {
      model_exprs[[length(model_exprs) + 1]] <- quote(lambda_b <- 1 / sqrt(w_lambda_b))
      model_exprs[[length(model_exprs) + 1]] <- quote(tau_hs_b <- tau0 / sqrt(w_tau_b))
      if (link == "ordered") {
        model_exprs[[length(model_exprs) + 1]] <- quote(lambda_tilde_b <- sqrt(c2_b * lambda_b^2 / (c2_b + tau_hs_b^2 * lambda_b^2)))
        model_exprs[[length(model_exprs) + 1]] <- quote(b <- z_b * tau_hs_b * lambda_tilde_b)
      } else {
        model_exprs[[length(model_exprs) + 1]] <- quote(b <- matrix(0, K_prob, K-1))
        model_exprs[[length(model_exprs) + 1]] <- quote(for (k in 1:(K-1)) {
          lambda_tilde_k <- sqrt(c2_b[k] * lambda_b[, k]^2 / (c2_b[k] + tau_hs_b[k]^2 * lambda_b[, k]^2))
          b[, k] <- z_b[, k] * tau_hs_b[k] * lambda_tilde_k
        })
      }
    } else if (regularization == "ssp") {
      model_exprs[[length(model_exprs) + 1]] <- quote(b <- b_raw * r_b * tau_b)
    }
  }

  model_exprs[[length(model_exprs) + 1]] <- quote(mu <- t(mu_p))

  if (has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- quote(pi_mat <- matrix(mu[1] * 0, N, K))
    if (link == "ordered") {
      model_exprs[[length(model_exprs) + 1]] <- quote(eta <- X_prob %*% b)
      model_exprs[[length(model_exprs) + 1]] <- quote(for (i in 1:N) {
        F_prev <- 0
        for (k in 1:(K-1)) {
          F_k <- inv_logit(cutpoints[k] - eta[i])
          pi_mat[i, k] <- F_k - F_prev
          F_prev <- F_k
        }
        pi_mat[i, K] <- 1 - F_prev
      })
    } else {
      # Sequential with Step-specific b
      model_exprs[[length(model_exprs) + 1]] <- quote(eta_mat <- X_prob %*% b)
      model_exprs[[length(model_exprs) + 1]] <- quote(for (i in 1:N) {
        q_cum <- 1
        for (k in 1:(K-1)) {
          q_k <- inv_logit(alpha[k] + eta_mat[i, k])
          pi_mat[i, k] <- q_cum * (1 - q_k)
          q_cum <- q_cum * q_k
        }
        pi_mat[i, K] <- q_cum
      })
    }
    model_exprs[[length(model_exprs) + 1]] <- quote(log_pi_mat <- log(pi_mat + 1e-15))
  } else if (prob_smoothing) {
    model_exprs[[length(model_exprs) + 1]] <- quote(theta <- softmax(logit_theta))
  }

  model_exprs[[length(model_exprs) + 1]] <- quote(lp <- mu[1] * 0)
  model_exprs[[length(model_exprs) + 1]] <- quote(log_dens_mat <- matrix(mu[1] * 0, N, K))

  # Likelihood Calculation
  if (is_diag) {
    dist_body <- bquote({
      ld <- Y[1] * 0
      m_vec <- mu[k, ]
      s_vec <- if (.(is_sigma_equal)) sigma else sigma[k, ]
      for (p in 1:P) {
        ld <- ld + normal_lpdf(Y[, p], m_vec[p], s_vec[p], sum = FALSE)
      }
      log_dens_mat[, k] <- ld
    })
  } else {
    # Full covariance
    dist_body <- bquote({
      m_vec <- mu[k, ]
      s_vec <- if (.(is_sigma_equal)) sigma else sigma[k, ]
      L_mat <- if (.(covariance == "full")) matrix(L_corr[k, , ], P, P) else matrix(L_corr, P, P)
      log_dens_mat[, k] <- multi_normal_CF_lpdf(Y, mean = as.vector(m_vec), sd = as.vector(s_vec), CF_Omega = L_mat, sum = FALSE)
    })
  }
  model_exprs[[length(model_exprs) + 1]] <- as.call(list(as.name("for"), as.name("k"), quote(1:K), dist_body))

  if (has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- quote(lp <- lp + sum(log_sum_exp(log_pi_mat + log_dens_mat)))
  } else {
    model_exprs[[length(model_exprs) + 1]] <- quote(log_theta <- log(theta))
    model_exprs[[length(model_exprs) + 1]] <- quote(lp <- lp + sum(log_sum_exp(t(t(log_dens_mat) + log_theta))))
  }

  # GP Prior (Items)
  gp_mag <- if (is.null(magnitude)) quote(magnitude) else magnitude
  gp_sm  <- if (is.null(smoothing)) quote(smoothing) else smoothing
  model_exprs[[length(model_exprs) + 1]] <- bquote(for (p in 1:P) {
    lp <- lp + gaussian_process_lpdf(mu[, p], x = rank_coords, magnitude = .(gp_mag), smoothing = .(gp_sm), noise = .(noise))
  })

  if (prob_smoothing && !has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- bquote(lp <- lp + gaussian_process_lpdf(logit_theta, x = rank_coords, magnitude = magnitude_theta, smoothing = .(gp_sm), noise = .(noise)))
  }

  # Priors
  if (use_weak_info) {
    model_exprs[[length(model_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    if (is.null(magnitude)) model_exprs[[length(model_exprs) + 1]] <- bquote(magnitude ~ exponential(.(prior$mag_rate)))
    if (is.null(smoothing)) model_exprs[[length(model_exprs) + 1]] <- bquote(smoothing ~ exponential(.(prior$smooth_rate)))
    if (prob_smoothing && !has_cov_prob) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(magnitude_theta ~ exponential(.(prior$mag_rate)))
    }

    if (multivariate && !is_diag) {
      if (covariance == "full") {
        model_exprs[[length(model_exprs) + 1]] <- bquote(for (k in 1:K) matrix(L_corr[k, , ], P, P) ~ lkj_CF_corr(.(prior$lkj_eta)))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(L_corr ~ lkj_CF_corr(.(prior$lkj_eta)))
      }
    }

    if (has_cov_prob) {
      if (regularization == "rhs") {
        model_exprs[[length(model_exprs) + 1]] <- quote(z_b ~ normal(0, 1))
        model_exprs[[length(model_exprs) + 1]] <- quote(w_lambda_b ~ gamma(0.5, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- quote(w_tau_b ~ gamma(0.5, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- quote(c2_b ~ inverse_gamma(half_slab_df, half_slab_scale2))
      } else if (regularization == "ssp") {
        model_exprs[[length(model_exprs) + 1]] <- quote(b_raw ~ laplace(0, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- bquote(mu_r <- log(.(prior$ssp_ratio) / (1 - .(prior$ssp_ratio))))
        model_exprs[[length(model_exprs) + 1]] <- quote(r_b ~ logit_normal(mu_r, 3))
        model_exprs[[length(model_exprs) + 1]] <- quote(tau_b ~ exponential(1 / tau_scale))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
      }

      if (link == "ordered") {
        model_exprs[[length(model_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(alpha ~ normal(0, .(prior$cutpoint_sd)))
      }
    }
    if (!has_cov_prob && !prob_smoothing) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(theta ~ dirichlet(rep(1, K)))
    }
  } else if (prior_type == "normal") {
    if (!is.null(prior$sigma_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    if (has_cov_prob && !is.null(prior$b_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
    if (is.null(magnitude) && !is.null(prior$mag_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(magnitude ~ exponential(.(prior$mag_rate)))
    if (is.null(smoothing) && !is.null(prior$smooth_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(smoothing ~ exponential(.(prior$smooth_rate)))
    if (has_cov_prob) {
       if (link == "ordered") {
         if (!is.null(prior$cutpoint_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))
       } else {
         if (!is.null(prior$cutpoint_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(alpha ~ normal(0, .(prior$cutpoint_sd)))
       }
    }

    if (!is.null(prior$lkj_eta) && multivariate && !is_diag) {
      if (covariance == "full") {
        model_exprs[[length(model_exprs) + 1]] <- bquote(for (k in 1:K) matrix(L_corr[k, , ], P, P) ~ lkj_CF_corr(.(prior$lkj_eta)))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(L_corr ~ lkj_CF_corr(.(prior$lkj_eta)))
      }
    }
  }
  model_ast <- as.call(model_exprs)

  # --- 4. Generate ---
  generate_exprs <- list(as.name("{"))
  generate_exprs[[length(generate_exprs) + 1]] <- quote(mu <- t(mu_p))
  if (has_cov_prob) {
    if (regularization == "rhs") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(lambda_b <- 1 / sqrt(w_lambda_b))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(tau_hs_b <- tau0 / sqrt(w_tau_b))
      if (link == "ordered") {
        generate_exprs[[length(generate_exprs) + 1]] <- quote(lambda_tilde_b <- sqrt(c2_b * lambda_b^2 / (c2_b + tau_hs_b^2 * lambda_b^2)))
        generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- z_b * tau_hs_b * lambda_tilde_b)
      } else {
        generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- matrix(0, K_prob, K-1))
        generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:(K-1)) {
          lambda_tilde_k <- sqrt(c2_b[k] * lambda_b[, k]^2 / (c2_b[k] + tau_hs_b[k]^2 * lambda_b[, k]^2))
          b[, k] <- z_b[, k] * tau_hs_b[k] * lambda_tilde_k
        })
      }
    } else if (regularization == "ssp") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- b_raw * r_b * tau_b)
    }
    if (link == "ordered") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(eta_mean <- X_means %*% b)
      generate_exprs[[length(generate_exprs) + 1]] <- quote(F_prev <- 0)
      generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- rtmb_vector(0, K))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:(K - 1)) {
        F_k <- inv_logit(cutpoints[k] - eta_mean)
        prob_mean[k] <- F_k - F_prev
        F_prev <- F_k
      })
      generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean[K] <- 1 - F_prev)
    } else {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(eta_mean_mat <- X_means %*% b)
      generate_exprs[[length(generate_exprs) + 1]] <- quote(q_cum <- 1)
      generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- rtmb_vector(0, K))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:(K - 1)) {
        q_k <- inv_logit(alpha[k] + eta_mean_mat[k])
        prob_mean[k] <- q_cum * (1 - q_k)
        q_cum <- q_cum * q_k
      })
      generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean[K] <- q_cum)
    }
  } else {
    if (prob_smoothing) generate_exprs[[length(generate_exprs) + 1]] <- quote(theta <- softmax(logit_theta))
    generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- theta)
  }

  if (!is_diag) {
    generate_exprs[[length(generate_exprs) + 1]] <- bquote(corr <- if (.(covariance == "full")) array(mu[1] * 0, dim = c(K, P, P)) else matrix(mu[1] * 0, P, P))
    if (covariance == "full") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:K) {
        L_mat <- matrix(L_corr[k, , ], P, P)
        corr[k, , ] <- L_mat %*% t(L_mat)
      })
    } else {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(L_mat <- matrix(L_corr, P, P))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(corr <- L_mat %*% t(L_mat))
    }
  }
  generate_base_exprs <- generate_exprs
  if (isTRUE(WAIC)) {
    waic_exprs <- list(as.name("{"))
    waic_exprs[[length(waic_exprs) + 1]] <- quote(log_dens_mat <- matrix(mu[1] * 0, N, K))
    waic_exprs[[length(waic_exprs) + 1]] <- as.call(list(as.name("for"), as.name("k"), quote(1:K), dist_body))
    if (has_cov_prob) {
      waic_exprs[[length(waic_exprs) + 1]] <- quote(pi_mat <- matrix(mu[1] * 0, N, K))
      if (link == "ordered") {
        waic_exprs[[length(waic_exprs) + 1]] <- quote(eta <- X_prob %*% b)
        waic_exprs[[length(waic_exprs) + 1]] <- quote(for (i in 1:N) {
          F_prev <- 0
          for (k in 1:(K - 1)) {
            F_k <- inv_logit(cutpoints[k] - eta[i])
            pi_mat[i, k] <- F_k - F_prev
            F_prev <- F_k
          }
          pi_mat[i, K] <- 1 - F_prev
        })
      } else {
        waic_exprs[[length(waic_exprs) + 1]] <- quote(eta_mat <- X_prob %*% b)
        waic_exprs[[length(waic_exprs) + 1]] <- quote(for (i in 1:N) {
          q_cum <- 1
          for (k in 1:(K - 1)) {
            q_k <- inv_logit(alpha[k] + eta_mat[i, k])
            pi_mat[i, k] <- q_cum * (1 - q_k)
            q_cum <- q_cum * q_k
          }
          pi_mat[i, K] <- q_cum
        })
      }
      waic_exprs[[length(waic_exprs) + 1]] <- quote(log_lik <- log_sum_exp(log(pi_mat + 1e-15) + log_dens_mat))
    } else {
      if (prob_smoothing) {
        waic_exprs[[length(waic_exprs) + 1]] <- quote(theta <- softmax(logit_theta))
      }
      waic_exprs[[length(waic_exprs) + 1]] <- quote(log_lik <- log_sum_exp(t(t(log_dens_mat) + log(theta))))
    }
  }
  generate_ast <- if (isTRUE(WAIC)) {
    .rtmb_waic_generate_ast(as.call(generate_base_exprs), as.call(waic_exprs))
  } else {
    as.call(generate_base_exprs)
  }

  mdl_code <- list(setup = setup_ast, parameters = param_ast, model = model_ast, generate = generate_ast, env = parent.frame())
  class(mdl_code) <- "rtmb_code"

  v_names <- list(mu = list(paste0("Rank", 1:K_mix), colnames(Y_mat)))
  v_names$sigma <- if (is_sigma_equal) colnames(Y_mat) else list(paste0("Rank", 1:K_mix), colnames(Y_mat))

  if (!is_diag) {
    if (covariance == "full") {
      v_names$L_corr <- list(paste0("Rank", 1:K_mix), colnames(Y_mat), colnames(Y_mat))
      v_names$corr <- list(paste0("Rank", 1:K_mix), colnames(Y_mat), colnames(Y_mat))
    } else {
      v_names$L_corr <- list(colnames(Y_mat), colnames(Y_mat))
      v_names$corr <- list(colnames(Y_mat), colnames(Y_mat))
    }
  }

  if (has_cov_prob) {
    if (link == "ordered") {
      v_names$b <- colnames(X_prob)
      v_names$cutpoints <- paste0("C", 1:(K_mix - 1))
    } else {
      step_labels <- paste0(1:(K_mix - 1), "->", 2:K_mix)
      v_names$b <- list(colnames(X_prob), step_labels)
      v_names$alpha <- step_labels
    }
  }

  data_list <- if (setup_from_formula) {
    list(df = setup_df, formula = formula, K = K_mix, rank_coords = rank_coords)
  } else {
    list(Y = Y_mat, N = N_obs, K = K_mix, P = P_dim, rank_coords = rank_coords)
  }
  if (has_cov_prob) {
    if (!setup_from_formula) {
      data_list$X_prob <- X_prob
      data_list$K_prob <- K_prob
      data_list$X_means <- as.vector(colMeans(X_prob))
    }
    if (regularization == "rhs") {
      p0 <- min(prior$expected_vars, K_prob - 1)
      if (p0 < 1) p0 <- 1
      if (!setup_from_formula) {
        data_list$tau0 <- p0 / (K_prob - p0) / sqrt(N_obs)
        data_list$half_slab_df <- prior$slab_df / 2
        data_list$half_slab_scale2 <- prior$slab_scale^2 / 2
      }
    } else if (regularization == "ssp") {
      if (!setup_from_formula) data_list$tau_scale <- prior$max_beta / 1.96
    }
  }
  mdl_code$setup_env <- .rtmb_setup_env(environment(), setup_ast, exclude = names(data_list))

  init_list <- list()
  if (nrow(Y_mat) >= K_mix) {
    # PCA-based ordering for LRT initial values
    # We use PC1 to capture the dominant direction of variation
    pc1 <- try(prcomp(Y_mat, rank. = 1)$x[, 1], silent = TRUE)
    if (!inherits(pc1, "try-error")) {
      # Split by quantiles of PC1
      groups <- cut(pc1, breaks = quantile(pc1, probs = seq(0, 1, length.out = K_mix + 1)), include.lowest = TRUE, labels = FALSE)
      mu_init_p <- matrix(0, P_dim, K_mix)
      for (k in 1:K_mix) {
        idx <- which(groups == k)
        if (length(idx) > 0) {
          mu_init_p[, k] <- colMeans(Y_mat[idx, , drop = FALSE])
        } else {
          mu_init_p[, k] <- colMeans(Y_mat)
        }
      }
      # Initial mu_p: handle potential NAs from colMeans
      mu_init_p[is.na(mu_init_p)] <- 0
      # Ensure it's ordered along the average direction
      if (mean(mu_init_p[, 1]) > mean(mu_init_p[, K_mix])) {
        mu_init_p <- mu_init_p[, K_mix:1]
      }
      # Add small jitter to ensure strict ordering for AD setup
      # AND handle any individual items that might be decreasing
      for (k in 2:K_mix) {
        for (p in 1:P_dim) {
          if (mu_init_p[p, k] <= mu_init_p[p, k-1]) {
            mu_init_p[p, k] <- mu_init_p[p, k-1] + 1e-3
          }
        }
      }
      init_list$mu_p <- mu_init_p

      # Initial b and cutpoints if covariates present
      if (has_cov_prob) {
        # Proxy regression of ranks on X
        fit_b <- try(lm(groups ~ X_prob), silent = TRUE)
        if (!inherits(fit_b, "try-error")) {
          b_est <- coef(fit_b)
          # Exclude intercept if it was included (lm usually does)
          if ("(Intercept)" %in% names(b_est)) b_est <- b_est[names(b_est) != "(Intercept)"]
          b_est[is.na(b_est)] <- 0

          if (regularization == "rhs") {
            init_list$z_b <- if (link == "ordered") as.vector(b_est) else matrix(b_est, ncol(X_prob), K_mix - 1)
            init_list$w_lambda_b <- if (link == "ordered") rep(1, length(init_list$z_b)) else matrix(1, ncol(X_prob), K_mix - 1)
            init_list$w_tau_b <- if (link == "ordered") 1 else rep(1, K_mix - 1)
            init_list$c2_b <- if (link == "ordered") 1 else rep(1, K_mix - 1)
          } else if (regularization == "ssp") {
            init_list$b_raw <- if (link == "ordered") as.vector(b_est) else matrix(b_est, ncol(X_prob), K_mix - 1)
            init_list$r_b <- if (link == "ordered") rep(0.5, length(init_list$b_raw)) else matrix(0.5, ncol(X_prob), K_mix - 1)
            init_list$tau_b <- if (link == "ordered") 1 else rep(1, K_mix - 1)
          } else {
            if (link == "ordered") {
              init_list$b <- as.vector(b_est)
            } else {
              init_list$b <- matrix(b_est, ncol(X_prob), K_mix - 1)
            }
          }
        }

        # Thresholds from group proportions
        props <- table(factor(groups, levels = 1:K_mix)) / length(groups)
        F_k <- cumsum(props)
        jitter <- (1:(K_mix - 1)) * 1e-4
        if (link == "ordered") {
          init_list$cutpoints <- qlogis(pmin(pmax(F_k[1:(K_mix-1)], 0.01), 0.99)) + jitter
        } else {
          # Sequential alphas: q_k = P(Z>k|Z>=k)
          q_vals <- (1 - F_k[1:(K_mix-1)]) / (1 - c(0, F_k)[1:(K_mix-1)])
          init_list$alpha <- qlogis(pmin(pmax(q_vals, 0.01), 0.99)) + jitter
        }
      }
    }
  }

  if (is.null(init_list$mu_p)) {
    mu_init <- matrix(0, P_dim, K_mix)
    for (p in 1:P_dim) {
      r_y <- range(Y_mat[, p])
      mu_init[p, ] <- seq(r_y[1], r_y[2], length.out = K_mix)
    }
    init_list$mu_p <- mu_init
  }

  if (!is.null(Y_mat)) {
    s_init <- apply(Y_mat, 2, sd)
    if (is_sigma_equal) {
      init_list$sigma <- s_init
    } else {
      init_list$sigma <- matrix(s_init, K_mix, P_dim, byrow = TRUE)
    }
  }

  view_order <- c("prob_mean", "mu", "sigma")
  if (has_cov_prob) {
    view_order <- if (link == "ordered") c("b", "cutpoints", view_order) else c("b", "alpha", view_order)
  }
  if (!is_diag) view_order <- c(view_order, "corr")

  mdl <- rtmb_model(data_list, mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order)
  mdl$type <- "lrt"
  mdl$extra <- list(
    source = "wrapper",
    prior_type = prior$type,
    marginal = "mu_p",
    two_stage = isTRUE(two_stage),
    lrt = list(
      formula = formula,
      data = if (setup_from_formula) setup_df else NULL,
      k = K_mix,
      rank_coords = rank_coords,
      covariance = covariance,
      magnitude = magnitude,
      smoothing = smoothing,
      noise = noise,
      prob_smoothing = prob_smoothing,
      link = link,
      prior = prior,
      y_range = y_range,
      optimize_defaults = wrapper_optimize_defaults,
      multivariate = multivariate,
      has_cov_prob = has_cov_prob
    )
  )
  return(mdl)
}

.lrt_response_only_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) < 3L) return(formula)
  stats::as.formula(call("~", formula[[2L]], 1), env = environment(formula))
}

.lrt_log_dmvn <- function(Y, mean, sigma, corr = NULL) {
  Y <- as.matrix(Y)
  P <- ncol(Y)
  if (is.null(corr)) corr <- diag(P)
  S <- diag(as.numeric(sigma), P, P) %*% corr %*% diag(as.numeric(sigma), P, P)
  S <- S + diag(1e-8, P)
  R <- chol(S)
  z <- backsolve(R, t(sweep(Y, 2, mean, "-")), transpose = TRUE)
  quad <- colSums(z^2)
  -0.5 * (P * log(2 * pi) + 2 * sum(log(diag(R))) + quad)
}

.lrt_rank_weights_from_fit <- function(Y, fit, covariance) {
  Y <- as.matrix(Y)
  mu <- t(fit$par$mu_p)
  K <- nrow(mu)
  P <- ncol(mu)
  theta <- fit$par$theta
  if (is.null(theta)) theta <- rep(1 / K, K)
  theta <- pmax(as.numeric(theta), 1e-12)
  theta <- theta / sum(theta)

  sigma <- fit$par$sigma
  is_sigma_equal <- covariance %in% c("diagonal_equal", "full_equal", "full_equal_corr")
  is_diag <- covariance %in% c("diagonal", "diagonal_equal")

  corr_obj <- fit$generate$corr
  log_dens <- matrix(NA_real_, nrow(Y), K)
  for (kk in seq_len(K)) {
    sig_k <- if (is_sigma_equal) as.numeric(sigma) else as.numeric(sigma[kk, ])
    if (is_diag) {
      ld <- rep(0, nrow(Y))
      for (pp in seq_len(P)) {
        ld <- ld + stats::dnorm(Y[, pp], mean = mu[kk, pp], sd = sig_k[pp], log = TRUE)
      }
      log_dens[, kk] <- ld
    } else {
      corr_k <- if (covariance == "full") {
        matrix(corr_obj[kk, , ], P, P)
      } else if (!is.null(corr_obj)) {
        matrix(corr_obj, P, P)
      } else {
        diag(P)
      }
      log_dens[, kk] <- .lrt_log_dmvn(Y, mean = mu[kk, ], sigma = sig_k, corr = corr_k)
    }
  }

  log_w <- sweep(log_dens, 2, log(theta), "+")
  mx <- apply(log_w, 1, max)
  W <- exp(log_w - mx)
  W / rowSums(W)
}

.make_lrt_rank_regression_model <- function(W, X_prob, link, prior, init = NULL) {
  W <- as.matrix(W)
  X_prob <- as.matrix(X_prob)
  K <- ncol(W)
  K_prob <- ncol(X_prob)
  link <- match.arg(link, c("ordered", "sequential"))
  prior_type <- prior$type %||% "flat"

  if (prior_type %in% c("rhs", "ssp")) {
    stop("two_stage = TRUE currently does not support prior_rhs() or prior_ssp().", call. = FALSE)
  }

  setup_ast <- quote({
    N <- nrow(W)
    K <- ncol(W)
    K_prob <- ncol(X_prob)
    X_means <- matrix(colMeans(X_prob), 1, K_prob)
  })

  if (link == "ordered") {
    param_ast <- quote({
      b <- Dim(K_prob)
      cutpoints <- Dim(K - 1, type = "ordered")
    })
    model_exprs <- list(
      as.name("{"),
      quote(lp <- b[1] * 0),
      quote(for (i in 1:N) {
        eta_i <- sum(X_prob[i, ] * b)
        F_prev <- 0
        for (k in 1:(K - 1)) {
          F_k <- inv_logit(cutpoints[k] - eta_i)
          p_k <- F_k - F_prev
          lp <- lp + W[i, k] * log(p_k + 1e-15)
          F_prev <- F_k
        }
        lp <- lp + W[i, K] * log(1 - F_prev + 1e-15)
      })
    )
    generate_ast <- quote({
      eta_mean <- X_means %*% b
      F_prev <- 0
      prob_mean <- rtmb_vector(0, K)
      for (k in 1:(K - 1)) {
        F_k <- inv_logit(cutpoints[k] - eta_mean)
        prob_mean[k] <- F_k - F_prev
        F_prev <- F_k
      }
      prob_mean[K] <- 1 - F_prev
    })
    par_names <- list(b = colnames(X_prob), cutpoints = paste0("C", seq_len(K - 1)))
  } else {
    param_ast <- quote({
      b <- Dim(c(K_prob, K - 1))
      alpha <- Dim(K - 1)
    })
    model_exprs <- list(
      as.name("{"),
      quote(lp <- b[1] * 0),
      quote(for (i in 1:N) {
        eta_i <- X_prob[i, ] %*% b
        q_cum <- 1
        for (k in 1:(K - 1)) {
          q_k <- inv_logit(alpha[k] + eta_i[k])
          p_k <- q_cum * (1 - q_k)
          lp <- lp + W[i, k] * log(p_k + 1e-15)
          q_cum <- q_cum * q_k
        }
        lp <- lp + W[i, K] * log(q_cum + 1e-15)
      })
    )
    generate_ast <- quote({
      eta_mean_mat <- X_means %*% b
      q_cum <- 1
      prob_mean <- rtmb_vector(0, K)
      for (k in 1:(K - 1)) {
        q_k <- inv_logit(alpha[k] + eta_mean_mat[k])
        prob_mean[k] <- q_cum * (1 - q_k)
        q_cum <- q_cum * q_k
      }
      prob_mean[K] <- q_cum
    })
    step_labels <- paste0(seq_len(K - 1), "->", 2:K)
    par_names <- list(b = list(colnames(X_prob), step_labels), alpha = step_labels)
  }

  if (prior_type %in% c("weak", "normal")) {
    if (!is.null(prior$b_sd)) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
    }
    if (!is.null(prior$cutpoint_sd)) {
      if (link == "ordered") {
        model_exprs[[length(model_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(alpha ~ normal(0, .(prior$cutpoint_sd)))
      }
    }
  }

  data_list <- list(W = W, X_prob = X_prob)
  code <- list(
    setup = setup_ast,
    parameters = param_ast,
    model = as.call(model_exprs),
    generate = generate_ast,
    env = parent.frame()
  )
  class(code) <- "rtmb_code"

  if (is.null(init)) {
    init <- list(b = if (link == "ordered") rep(0, K_prob) else matrix(0, K_prob, K - 1))
    avg_w <- pmin(pmax(cumsum(colMeans(W))[seq_len(K - 1)], 0.01), 0.99)
    if (link == "ordered") {
      init$cutpoints <- stats::qlogis(avg_w) + seq_len(K - 1) * 1e-4
    } else {
      q_vals <- (1 - avg_w) / pmax(1 - c(0, avg_w)[seq_len(K - 1)], 1e-6)
      init$alpha <- stats::qlogis(pmin(pmax(q_vals, 0.01), 0.99)) + seq_len(K - 1) * 1e-4
    }
  }

  code$setup_env <- .rtmb_setup_env(environment(), setup_ast, exclude = names(data_list))
  obj <- rtmb_model(data_list, code, par_names = par_names, init = init, view = c("b", "cutpoints", "alpha", "prob_mean"))
  obj$type <- "lrt_two_stage_regression"
  obj$extra <- list(source = "wrapper", prior_type = prior_type, marginal = "b")
  obj
}

.lrt_two_stage_fit_regression <- function(W, X_prob, link, prior, init = NULL,
                                          optimizer = "nlminb", method = "BFGS",
                                          control = list(), se_method = "wald",
                                          num_estimate = 1) {
  reg_model <- .make_lrt_rank_regression_model(W, X_prob, link, prior, init = init)
  reg_model$optimize(
    laplace = FALSE,
    num_estimate = num_estimate,
    optimizer = optimizer,
    method = method,
    control = control,
    se_method = se_method,
    df_method = "none"
  )
}

.lrt_two_stage_optimize <- function(self, laplace = TRUE, init = NULL, num_estimate = 1,
                                    control = list(), optimizer = "nlminb", method = "BFGS",
                                    map = NULL, fixed = NULL, se_method = c("wald", "sampling", "none"),
                                    num_samples = 1000, seed = 123, marginal = "none",
                                    df_method = "auto", view = NULL, ...) {
  se_method <- match.arg(se_method)
  if (identical(se_method, "sampling")) {
    stop("two_stage = TRUE currently supports se_method = 'wald' or 'none' only.", call. = FALSE)
  }
  if (!is.null(fixed) || !is.null(map)) {
    stop("two_stage = TRUE currently does not support fixed or map in optimize().", call. = FALSE)
  }

  cfg <- self$extra$lrt
  if (is.null(cfg) || !isTRUE(cfg$has_cov_prob)) {
    stop("two_stage estimation requires an rtmb_lrt model with covariates.", call. = FALSE)
  }
  if (identical(num_estimate, 1) && !is.null(cfg$optimize_defaults$num_estimate)) {
    num_estimate <- cfg$optimize_defaults$num_estimate
  }

  dot_args <- list(...)
  delta_eps <- dot_args$two_stage_delta_eps %||% 1e-4
  delta_method <- isTRUE(dot_args$two_stage_delta %||% TRUE)

  stage1_formula <- .lrt_response_only_formula(cfg$formula)
  if (!getOption("BayesRTMB.silent", FALSE)) {
    cat("Starting two-stage LRT optimization...\n")
    cat("Stage 1/2: estimating latent ranks without covariates...\n")
  }
  stage1_model <- rtmb_lrt(
    formula = stage1_formula,
    k = cfg$k,
    data = cfg$data,
    rank_coords = cfg$rank_coords,
    covariance = cfg$covariance,
    magnitude = cfg$magnitude,
    smoothing = cfg$smoothing,
    noise = cfg$noise,
    prob_smoothing = FALSE,
    link = cfg$link,
    prior = cfg$prior,
    y_range = cfg$y_range,
    two_stage = FALSE
  )
  stage1_fit <- stage1_model$optimize(
    laplace = laplace,
    init = init,
    num_estimate = num_estimate,
    control = control,
    optimizer = optimizer,
    method = method,
    se_method = se_method,
    num_samples = num_samples,
    seed = seed,
    marginal = marginal,
    df_method = df_method,
    view = view
  )

  Y <- self$data$Y
  X_prob <- self$data$X_prob
  W <- .lrt_rank_weights_from_fit(Y, stage1_fit, cfg$covariance)

  if (!getOption("BayesRTMB.silent", FALSE)) {
    cat("Stage 2/2: estimating rank regression with fixed soft rank weights...\n")
  }
  stage2_fit <- .lrt_two_stage_fit_regression(
    W = W,
    X_prob = X_prob,
    link = cfg$link,
    prior = cfg$prior,
    optimizer = optimizer,
    method = method,
    control = control,
    se_method = se_method,
    num_estimate = num_estimate
  )

  delta_applied <- FALSE
  delta_se_cond <- delta_se_extra <- NULL
  if (delta_method && identical(se_method, "wald") && !is.null(stage1_fit$vcov_unc)) {
    if (!getOption("BayesRTMB.silent", FALSE)) {
      cat("Applying delta-method correction for Stage 1 uncertainty...\n")
    }
    base_est <- stage2_fit$df_fixed$Estimate
    active <- stage1_fit$idx_fix_active %||% seq_along(stage1_fit$par_unc)
    V_eta <- stage1_fit$vcov_unc[active, active, drop = FALSE]
    J <- matrix(0, length(base_est), length(active))
    init_stage2 <- stage2_fit$par
    old_silent <- options(BayesRTMB.silent = TRUE)
    on.exit(options(old_silent), add = TRUE)

    for (jj in seq_along(active)) {
      u_tmp <- stage1_fit$par_unc
      step <- delta_eps * (abs(u_tmp[active[jj]]) + 1)
      u_tmp[active[jj]] <- u_tmp[active[jj]] + step
      par_tmp <- stage1_model$to_constrained(unconstrained_vector_to_list(u_tmp, stage1_model$par_list))
      gen_tmp <- if (!is.null(stage1_model$generate)) {
        tryCatch(stage1_model$generate(stage1_model$data, par_tmp), error = function(e) stage1_fit$generate)
      } else {
        stage1_fit$generate
      }
      fit_tmp <- list(par = par_tmp, generate = gen_tmp)
      W_tmp <- .lrt_rank_weights_from_fit(Y, fit_tmp, cfg$covariance)
      reg_tmp <- tryCatch(
        .lrt_two_stage_fit_regression(
          W = W_tmp,
          X_prob = X_prob,
          link = cfg$link,
          prior = cfg$prior,
          init = init_stage2,
          optimizer = optimizer,
          method = method,
          control = control,
          se_method = "none",
          num_estimate = 1
        ),
        error = function(e) NULL
      )
      if (!is.null(reg_tmp) && !is.null(reg_tmp$df_fixed)) {
        J[, jj] <- (reg_tmp$df_fixed$Estimate - base_est) / step
      }
    }

    V_extra <- J %*% V_eta %*% t(J)
    se_cond <- stage2_fit$df_fixed$`Std. Error`
    delta_se_cond <- se_cond
    delta_se_extra <- sqrt(pmax(diag(V_extra), 0))
    se_total <- sqrt(pmax(se_cond^2 + diag(V_extra), 0))
    stage2_fit$df_fixed$`Std. Error` <- se_total
    stage2_fit$df_fixed$`Lower 95%` <- base_est - stats::qnorm(0.975) * se_total
    stage2_fit$df_fixed$`Upper 95%` <- base_est + stats::qnorm(0.975) * se_total
    delta_applied <- TRUE
  }

  stage2_fit$model$type <- "lrt"
  stage2_fit$model$extra$two_stage <- TRUE
  stage2_fit$model$extra$two_stage_details <- list(
    stage1_fit = stage1_fit,
    weights = W,
    stage1_objective = stage1_fit$objective,
    stage2_objective = stage2_fit$objective,
    correction = if (delta_applied) "delta" else "none",
    se_conditional = delta_se_cond,
    se_stage1 = delta_se_extra
  )
  class(stage2_fit) <- unique(c("rtmb_lrt_two_stage", class(stage2_fit)))
  stage2_fit
}

#' @export
print.rtmb_lrt_two_stage <- function(x, ...) {
  x$print(...)
  details <- x$model$extra$two_stage_details %||% list()
  cat("\nTwo-stage LRT note:\n")
  cat("  The displayed objective is the Stage 2 rank-regression objective.\n")
  cat("  It is not directly comparable to the joint LRT objective.\n")
  if (!is.null(details$correction)) {
    cat(sprintf("  Stage 1 uncertainty correction: %s\n", details$correction))
  }
  if (!is.null(details$stage1_objective)) {
    cat(sprintf("  Stage 1 latent-rank objective: %.2f\n", details$stage1_objective))
  }
  if (!is.null(details$stage2_objective)) {
    cat(sprintf("  Stage 2 rank-regression objective: %.2f\n", details$stage2_objective))
  }
  invisible(x)
}
