#' Mixture Model Wrapper for RTMB
#'
#' @description
#' Provides a user-friendly interface for fitting Gaussian mixture models
#' with optional covariates on class membership probabilities and various
#' covariance structures.
#'
#' @param formula A formula specifying the response variable(s). For multivariate, use \code{cbind(y1, y2) ~ 1}.
#' @param k Number of mixture components.
#' @param data A data frame containing the variables in the model.
#' @param covariance Covariance structure: "diagonal" (default), "diagonal_equal", "full", "full_equal", or "full_equal_corr".
#' @param prior Prior configuration: `prior_flat()`, `prior_normal()`,
#'   `prior_weak()`, `prior_rhs()`, or `prior_ssp()`. Default is
#'   `prior_flat()`. If `y_range` is supplied with the default flat prior,
#'   the wrapper automatically switches to `prior_weak()`.
#' @param y_range Optional numeric vector or matrix defining the theoretical range (min, max) of response variables.
#'   Specifying this automatically enables weakly informative priors if `prior` is `prior_flat()`.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Additional arguments passed to `rtmb_model`.
#' @return A \code{RTMB_Model} object.
#' @example inst/examples/ex_mixture.R
#' @export
rtmb_mixture <- function(formula, k = 2, data = NULL,
                         covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
                         prior = prior_flat(), y_range = NULL, fixed = NULL,
                         WAIC = FALSE, ...) {

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2 || k != as.integer(k)) {
    stop("'k' must be an integer greater than or equal to 2.", call. = FALSE)
  }
  k <- as.integer(k)

  if (is.null(prior)) {
    prior <- prior_flat()
  }

  # Automatically switch to prior_weak() if y_range is provided and prior is default flat
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && identical(prior$type, "flat")) {
    prior <- prior_weak()
  }

  if (!inherits(prior, "rtmb_prior")) {
    stop(
      "prior must be an object of class 'rtmb_prior'. ",
      "Use prior_flat(), prior_normal(), prior_weak(), prior_rhs(), or prior_ssp().",
      call. = FALSE
    )
  }

  # NSE for formula: handle case where formula is just a variable name in data
  formula_expr <- substitute(formula)
  formula_val <- try(formula, silent = TRUE)

  if (inherits(formula_val, "try-error") ||
      (!inherits(formula_val, "formula") && !is.matrix(formula_val) && !is.vector(formula_val))) {
    # Check if variables in formula_expr are in data
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
      dirichlet_alpha = NULL
    )
    prior <- .merge_prior(normal_defaults, prior)
  } else {
    default_prior <- list(Intercept_sd = 10, mu_sd = 10, b_sd = 10, sigma_rate = 1, lkj_eta = 1.0, dirichlet_alpha = 1.0)
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
    if (prior_type == "weak" && is.null(y_range)) {
      warning(
        "prior_weak() is used without 'y_range'. ",
        "Default hyperparameters will be applied, which may not be well-scaled for your data. ",
        "Consider specifying 'y_range' for better prior calibration.",
        call. = FALSE
      )
    }
  }

  setup_from_formula <- inherits(formula, "formula") && !is.null(data)

  # Parse formulas and prepare data
  if (is.character(formula) && !is.null(data) && all(formula %in% names(data))) {
    Y_mat <- as.matrix(data[, formula, drop = FALSE])
    X_prob <- NULL
  } else if (is.matrix(formula) || (is.vector(formula) && !is.character(formula))) {
    Y_mat <- as.matrix(formula)
    X_prob <- NULL
  } else if (!inherits(formula, "formula")) {
    # Case where a single variable name is passed
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
  if (!is.numeric(Y_mat) && !is.logical(Y_mat)) {
    stop("The mixture response must be numeric.", call. = FALSE)
  }
  storage.mode(Y_mat) <- "double"
  if (anyNA(Y_mat)) {
    stop("rtmb_mixture() does not currently support missing response values.", call. = FALSE)
  }
  N_obs <- nrow(Y_mat)
  P_dim <- ncol(Y_mat)
  if (is.null(N_obs) || N_obs < 1L) {
    stop("The mixture response must contain at least one observation.", call. = FALSE)
  }
  if (is.null(P_dim) || P_dim < 1L) {
    stop("The mixture response must contain at least one variable.", call. = FALSE)
  }

  if (is.null(colnames(Y_mat))) {
    colnames(Y_mat) <- paste0("Y", 1:P_dim)
  }

  multivariate <- P_dim > 1
  has_cov_prob <- !is.null(X_prob)
  K_prob <- if (has_cov_prob) ncol(X_prob) else 0
  is_sigma_equal <- covariance %in% c("diagonal_equal", "full_equal", "full_equal_corr")

  # --- 1. Dynamic AST Construction: setup ---
  if (setup_from_formula) {
    setup_exprs <- list(
      as.name("{"),
      quote(mf <- model.frame(formula, df)),
      quote(Y <- model.response(mf)),
      quote(if (is.null(Y)) Y <- as.matrix(mf)),
      quote(Y <- as.matrix(Y)),
      quote(N <- nrow(Y)),
      quote(P <- ncol(Y)),
      quote(K <- K)
    )
  } else {
    setup_exprs <- list(
      as.name("{"),
      quote(N <- N),
      quote(K <- K),
      quote(P <- P)
    )
  }
  if (has_cov_prob) {
    if (setup_from_formula) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_prob <- model.matrix(formula, mf))
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

  # --- 2. Dynamic AST Construction: parameters ---
  param_exprs <- list(as.name("{"))
  if (has_cov_prob) {
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
  } else {
    param_exprs[[length(param_exprs) + 1]] <- bquote(theta <- Dim(.(K_mix), type = "simplex"))
  }

  if (multivariate) {
    param_exprs[[length(param_exprs) + 1]] <- bquote(mu <- Dim(c(K, P)))
    s_dim <- if (is_sigma_equal && covariance != "diagonal_equal") P_dim else c(K_mix, P_dim)
    if (covariance == "diagonal_equal") s_dim <- P_dim
    param_exprs[[length(param_exprs) + 1]] <- bquote(sigma <- Dim(.(s_dim), lower = 0))
  } else {
    param_exprs[[length(param_exprs) + 1]] <- bquote(mu <- Dim(K))
    s_dim <- if (is_sigma_equal) NULL else K_mix
    param_exprs[[length(param_exprs) + 1]] <- if (is.null(s_dim)) quote(sigma <- Dim(lower = 0)) else bquote(sigma <- Dim(.(s_dim), lower = 0))
  }

  if (multivariate && covariance %in% c("full", "full_equal", "full_equal_corr")) {
    if (covariance == "full") {
      param_exprs[[length(param_exprs) + 1]] <- quote(L_corr <- Dim(c(K, P, P), type = "CF_corr"))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- quote(L_corr <- Dim(c(P, P), type = "CF_corr"))
    }
  }
  param_ast <- as.call(param_exprs)

  # --- 3. Dynamic AST Construction: model ---
  model_exprs <- list(as.name("{"))

  # Regularization transformations
  if (has_cov_prob) {
    if (regularization == "rhs") {
      model_exprs[[length(model_exprs) + 1]] <- quote(lambda_b <- 1 / sqrt(w_lambda_b))
      model_exprs[[length(model_exprs) + 1]] <- quote(tau_hs_b <- tau0 / sqrt(w_tau_b))
      model_exprs[[length(model_exprs) + 1]] <- quote(b <- matrix(0, K_prob, K-1)) 
      model_exprs[[length(model_exprs) + 1]] <- quote(for (k in 1:(K-1)) {
        lambda_tilde_k <- sqrt(c2_b[k] * lambda_b[, k]^2 / (c2_b[k] + tau_hs_b[k]^2 * lambda_b[, k]^2))
        b[, k] <- z_b[, k] * tau_hs_b[k] * lambda_tilde_k
      })
    } else if (regularization == "ssp") {
      model_exprs[[length(model_exprs) + 1]] <- quote(b <- b_raw * r_b * tau_b)
    }
  }

  model_exprs[[length(model_exprs) + 1]] <- quote(lp <- mu[1] * 0)
  model_exprs[[length(model_exprs) + 1]] <- quote(log_dens_mat <- matrix(mu[1] * 0, N, K))

  is_diag <- covariance %in% c("diagonal", "diagonal_equal")
  s_k_expr <- if (is_sigma_equal) quote(sigma) else quote(sigma[k, ])
  if (!multivariate && !is_sigma_equal) s_k_expr <- quote(sigma[k])

  dist_body <- if (multivariate) {
    if (is_diag) {
      bquote({
        ld <- Y[1] * 0
        s_vec <- .(s_k_expr)
        m_vec <- mu[k, ]
        for (p in 1:P) {
          ld <- ld + normal_lpdf(Y[, p], m_vec[p], s_vec[p], sum = FALSE)
        }
        log_dens_mat[, k] <- ld
      })
    } else {
      L_corr_k_expr <- if (covariance == "full") quote(matrix(L_corr[k, , ], P, P)) else quote(L_corr)
      bquote({
        log_dens_mat[, k] <- multi_normal_CF_lpdf(Y, mean = mu[k, ], sd = .(s_k_expr), CF_Omega = .(L_corr_k_expr), sum = FALSE)
      })
    }
  } else {
    bquote({
      log_dens_mat[, k] <- normal_lpdf(Y, mu[k], .(s_k_expr), sum = FALSE)
    })
  }

  model_exprs[[length(model_exprs) + 1]] <- as.call(list(as.name("for"), as.name("k"), quote(1:K), dist_body))

  if (has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- quote(eta_prob <- X_prob %*% b)
    model_exprs[[length(model_exprs) + 1]] <- quote(log_pi_mat <- matrix(eta_prob[1] * 0, N, K))
    model_exprs[[length(model_exprs) + 1]] <- quote(for (i in 1:N) {
      log_pi_mat[i, ] <- log_softmax(c(0, eta_prob[i, ]))
    })
    model_exprs[[length(model_exprs) + 1]] <- quote(lp <- lp + sum(log_sum_exp(log_pi_mat + log_dens_mat)))
  } else {
    model_exprs[[length(model_exprs) + 1]] <- quote(log_theta <- log(theta))
    model_exprs[[length(model_exprs) + 1]] <- quote(lp <- lp + sum(log_sum_exp(t(t(log_dens_mat) + log_theta))))
  }

  # Priors
  if (use_weak_info) {
    # Component priors
    model_exprs[[length(model_exprs) + 1]] <- bquote(mu ~ normal(0, .(prior$Intercept_sd)))
    model_exprs[[length(model_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))

    if (multivariate && covariance %in% c("full", "full_equal", "full_equal_corr")) {
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
    } else {
      model_exprs[[length(model_exprs) + 1]] <- bquote(theta ~ dirichlet(rep(.(prior$dirichlet_alpha), K)))
    }
    } else if (prior_type == "normal") {
      if (!is.null(prior$Intercept_sd) || !is.null(prior$mu_sd)) {
        mu_sd_val <- if (!is.null(prior$mu_sd)) prior$mu_sd else prior$Intercept_sd
        model_exprs[[length(model_exprs) + 1]] <- bquote(mu ~ normal(0, .(mu_sd_val)))
      }

      if (!is.null(prior$sigma_rate)) {
        model_exprs[[length(model_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
      }

      if (has_cov_prob && !is.null(prior$b_sd)) {
        model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
      }

      if (!is.null(prior$lkj_eta) && multivariate && !is_diag) {
        if (covariance == "full") {
          model_exprs[[length(model_exprs) + 1]] <-
            bquote(for (k in 1:K) matrix(L_corr[k, , ], P, P) ~ lkj_CF_corr(.(prior$lkj_eta)))
        } else {
          model_exprs[[length(model_exprs) + 1]] <-
            bquote(L_corr ~ lkj_CF_corr(.(prior$lkj_eta)))
        }
      }
    }

  model_ast <- as.call(model_exprs)

  # --- 4. Dynamic AST Construction: generate ---
  generate_exprs <- list(as.name("{"))
  if (has_cov_prob) {
    if (regularization == "rhs") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(lambda_b <- 1 / sqrt(w_lambda_b))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(tau_hs_b <- tau0 / sqrt(w_tau_b))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- matrix(0, P, K-1))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:(K-1)) {
        lambda_tilde_k <- sqrt(c2_b[k] * lambda_b[, k]^2 / (c2_b[k] + tau_hs_b[k]^2 * lambda_b[, k]^2))
        b[, k] <- z_b[, k] * tau_hs_b[k] * lambda_tilde_k
      })
    } else if (regularization == "ssp") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- b_raw * r_b * tau_b)
    }
    generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- softmax(c(0, X_means %*% b)))
  } else {
    generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- theta)
  }

  if (multivariate && !is_diag) {
    if (covariance == "full") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(corr <- array(mu[1] * 0, dim = c(K, P, P)))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:K) {
        corr[k, , ] <- matrix(L_corr[k, , ], P, P) %*% t(matrix(L_corr[k, , ], P, P))
      })
    } else {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(corr <- L_corr %*% t(L_corr))
    }
  }
  generate_base_exprs <- generate_exprs
  if (isTRUE(WAIC)) {
    waic_exprs <- list(as.name("{"))
    waic_exprs[[length(waic_exprs) + 1]] <- quote(log_dens_mat <- matrix(mu[1] * 0, N, K))
    waic_exprs[[length(waic_exprs) + 1]] <- as.call(list(as.name("for"), as.name("k"), quote(1:K), dist_body))
    if (has_cov_prob) {
      waic_exprs[[length(waic_exprs) + 1]] <- quote(eta_prob <- X_prob %*% b)
      waic_exprs[[length(waic_exprs) + 1]] <- quote(log_pi_mat <- matrix(eta_prob[1] * 0, N, K))
      waic_exprs[[length(waic_exprs) + 1]] <- quote(for (i in 1:N) {
        log_pi_mat[i, ] <- log_softmax(c(0, eta_prob[i, ]))
      })
      waic_exprs[[length(waic_exprs) + 1]] <- quote(log_lik <- log_sum_exp(log_pi_mat + log_dens_mat))
    } else {
      waic_exprs[[length(waic_exprs) + 1]] <- quote(log_lik <- log_sum_exp(t(t(log_dens_mat) + log(theta))))
    }
  }
  generate_ast <- if (isTRUE(WAIC)) {
    .rtmb_waic_generate_ast(as.call(generate_base_exprs), as.call(waic_exprs))
  } else {
    as.call(generate_base_exprs)
  }

  mdl_code <- list(
    setup = setup_ast,
    parameters = param_ast,
    model = model_ast,
    generate = generate_ast,
    env = parent.frame()
  )
  class(mdl_code) <- "rtmb_code"

  v_names <- list()
  if (multivariate) {
    v_names$mu <- list(paste0("C", 1:K_mix), colnames(Y_mat))
    v_names$sigma <- if (is_sigma_equal && covariance != "diagonal_equal") colnames(Y_mat) else list(paste0("C", 1:K_mix), colnames(Y_mat))
    if (covariance == "diagonal_equal") v_names$sigma <- colnames(Y_mat)
  } else {
    v_names$mu <- paste0("C", 1:K_mix)
    v_names$sigma <- if (is_sigma_equal) "sigma" else paste0("C", 1:K_mix)
  }

  if (has_cov_prob) {
    colnames(X_prob)[colnames(X_prob) == "(Intercept)"] <- "Intercept"
    v_names$b <- list(colnames(X_prob), paste0("C", 2:K_mix))
  }
  if (multivariate && !is_diag) {
    if (covariance == "full") {
      v_names$L_corr <- list(paste0("C", 1:K_mix), colnames(Y_mat), colnames(Y_mat))
      v_names$corr <- list(paste0("C", 1:K_mix), colnames(Y_mat), colnames(Y_mat))
    } else {
      v_names$L_corr <- list(colnames(Y_mat), colnames(Y_mat))
      v_names$corr <- list(colnames(Y_mat), colnames(Y_mat))
    }
  }

  data_list <- if (setup_from_formula) {
    list(df = setup_df, formula = formula, K = K_mix)
  } else {
    list(Y = Y_mat, N = N_obs, K = K_mix, P = P_dim)
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

  init_list <- list()
  if (nrow(Y_mat) >= K_mix) {
    km <- try(kmeans(Y_mat, centers = K_mix, nstart = 5), silent = TRUE)
    if (!inherits(km, "try-error")) {
      if (multivariate) {
        init_list$mu <- km$centers
      } else {
        init_list$mu <- as.vector(km$centers)
      }
      props <- as.vector(table(factor(km$cluster, levels = 1:K_mix))) / nrow(Y_mat)
      init_list$theta <- props

      if (has_cov_prob) {
        b_init <- matrix(0, ncol(X_prob), K_mix - 1)
        for (k in 2:K_mix) {
          z_k <- as.numeric(km$cluster == k)
          fit_k <- try(suppressWarnings(glm(z_k ~ X_prob - 1 + offset(rep(qlogis(max(0.01, props[k])), length(z_k))), family = binomial)), silent = TRUE)
          if (!inherits(fit_k, "try-error")) {
            b_init[, k-1] <- coef(fit_k)
          }
        }
        if (regularization == "rhs") {
          init_list$z_b <- b_init
          init_list$w_lambda_b <- matrix(1, nrow(b_init), ncol(b_init))
          init_list$w_tau_b <- rep(1, ncol(b_init))
          init_list$c2_b <- rep(1, ncol(b_init))
        } else if (regularization == "ssp") {
          init_list$b_raw <- b_init
          init_list$r_b <- matrix(0.5, nrow(b_init), ncol(b_init))
          init_list$tau_b <- rep(1, ncol(b_init))
        } else {
          init_list$b <- b_init
        }
      }
    }
  }

  if (is.null(init_list$mu)) {
    if (multivariate) {
      mu_init <- matrix(0, K_mix, P_dim)
      for (p in 1:P_dim) {
        r_y <- range(Y_mat[, p])
        mu_init[, p] <- seq(r_y[1], r_y[2], length.out = K_mix + 2)[2:(K_mix + 1)]
      }
      init_list$mu <- mu_init
    } else {
      r_y <- range(Y_mat)
      init_list$mu <- seq(r_y[1], r_y[2], length.out = K_mix + 2)[2:(K_mix + 1)]
    }
  }

  if (!is.null(Y_mat)) {
     init_list$sigma <- if (is_sigma_equal) apply(Y_mat, 2, sd) else matrix(apply(Y_mat, 2, sd), K_mix, P_dim, byrow = TRUE)
  }

  view_order <- if (has_cov_prob) c("b", "prob_mean", "mu", "sigma") else c("prob_mean", "mu", "sigma")
  if (multivariate && !is_diag) view_order <- c(view_order, "corr")

  mdl <- rtmb_model(data_list, mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order)
  mdl$type <- "mixture"
  mdl$extra <- list(
    source = "wrapper",
    prior_type = prior$type,
    marginal = "mu"
  )
  return(mdl)
}
