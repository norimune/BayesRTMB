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
#' \strong{2. Regularization and Variable Selection:}
#' For rtmb_lm, rtmb_glm, and rtmb_glmer, you can handle high-dimensional predictors (where the number of variables is large relative to the sample size) using the penalty argument:
#' \itemize{
#'   \item \code{"none"}: Standard flat or weakly informative priors.
#'   \item \code{"rhs"}: Regularized Horseshoe prior for continuous shrinkage.
#'   \item \code{"ssp"}: Spike-and-Slab prior for sparse variable selection.
#' }
#' \emph{Note: When using regularization, you must specify \code{y_range = c(min, max)}
#' to let the model set appropriate global scales for the priors.}
#'
#' \strong{3. Weakly Informative Priors (\code{y_range}):}
#' By providing the theoretical range of your response variable via \code{y_range},
#' the wrappers automatically construct "Weakly Informative Priors". These priors
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


#' Specify a uniform or manual prior
#'
#' @param Intercept_sd Standard deviation for the intercept prior (Normal). Default is NULL (flat).
#' @param b_sd Standard deviation for the coefficients prior (Normal). Default is NULL (flat).
#' @param sigma_rate Rate for the residual standard deviation prior (Exponential). Default is NULL (flat).
#' @param tau_rate Rate for the random effects standard deviation prior (Exponential). Default is NULL (flat).
#' @param ... Optional hyperparameters
#' @return A list with class "rtmb_prior"
#' @export
prior_uniform <- function(Intercept_sd = NULL, b_sd = NULL, sigma_rate = NULL, tau_rate = NULL, ...) {
  res <- list(type = "uniform", Intercept_sd = Intercept_sd, b_sd = b_sd, sigma_rate = sigma_rate, tau_rate = tau_rate, ...)
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

#' RTMB-based GLMM wrapper function
#'
#' @param formula lme4-style formula (e.g., Y ~ X + (1 | GID))
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param laplace Logical; whether to marginalize random effects using Laplace approximation
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Use prior_weak(), prior_rhs(), or prior_ssp(). Default is NULL (flat prior).
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Required when using weakly informative or regularized priors with continuous models.
#' @param init List of initial values (generated automatically based on glm if omitted)
#' @param null Character string specifying the target parameter for the null model.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glmer <- function(formula, data, family = "gaussian", laplace = FALSE,
                       prior = prior_uniform(),
                       y_range = NULL,
                       init = NULL, null = NULL) {

  if (is.null(prior)) {
    prior <- prior_uniform()
  }

  # Automatically switch to prior_weak() if y_range is provided and prior is default uniform
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && prior$type == "uniform" &&
      is.null(prior$Intercept_sd) && is.null(prior$b_sd) &&
      is.null(prior$sigma_rate) && is.null(prior$tau_rate)) {
    prior <- prior_weak()
  }

  if (!inherits(prior, "rtmb_prior")) {
    stop("prior must be an object of class 'rtmb_prior'. Use prior_weak(), prior_rhs(), or prior_ssp().")
  }

  prior_type <- prior$type
  regularization <- if (prior_type %in% c("rhs", "ssp")) prior_type else "none"
  use_weak_info <- prior_type %in% c("weak", "rhs", "ssp")
  has_random <- !is.null(findbars(formula))

  default_prior <- list(Intercept_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, nu_rate = 0.1, cutpoint_sd = 2.5, shape_rate = 1.0, phi_rate = 1.0, lkj_eta = 1.0)
  prior <- modifyList(default_prior, prior)

  if (use_weak_info) {
    if (is.null(prior$max_beta)) prior$max_beta <- 1.0
    if (is.null(prior$sd_ratio)) prior$sd_ratio <- 0.5
    if (is.null(prior$expected_vars)) prior$expected_vars <- 3
    if (is.null(prior$slab_scale)) prior$slab_scale <- 2.0
    if (is.null(prior$slab_df)) prior$slab_df <- 4.0
    if (is.null(prior$ssp_ratio)) prior$ssp_ratio <- 0.25
  }

  families_requiring_yrange <- c("gaussian", "lognormal", "student_t")
  is_continuous <- family %in% families_requiring_yrange

  if (use_weak_info && is_continuous && is.null(y_range)) {
    stop(sprintf("Specifying 'y_range' is required when using prior_%s() with family = '%s'.", prior_type, family))
  }

  if (has_random) {
    all_vars_form <- subbars(formula)
    fixed_form <- nobars(formula)
    bars <- findbars(formula)

    # Expand bars with * or : using terms() expansion
    expanded_bars <- list()
    if (length(bars) > 0) {
      for (bar in bars) {
        grp_expr <- bar[[3]]
        # Use terms() to expand a*b into a, b, a:b
        term_labels <- attr(terms(as.formula(paste("~", deparse(grp_expr)))), "term.labels")
        for (label in term_labels) {
          new_bar <- bar
          new_bar[[3]] <- parse(text = label)[[1]]
          expanded_bars[[length(expanded_bars) + 1]] <- new_bar
        }
      }
    }
    bars <- expanded_bars
    num_bars <- length(bars)
    suffix <- function(b) if (num_bars > 1) paste0("_", b) else ""

    mf <- model.frame(all_vars_form, data = data)
    Y <- model.response(mf)
    X <- model.matrix(fixed_form, mf)
    offset <- model.offset(mf)
    N <- nrow(mf)

    dat_ranef <- list()
    num_ranef_list <- list()
    num_groups_list <- list()
    ranef_names_list <- list()
    group_labels_list <- list()

    for (b in 1:num_bars) {
      bar <- bars[[b]]
      re_form <- as.formula(paste("~", deparse(bar[[2]])))

      Z_mf <- model.matrix(re_form, mf)

      # Handle grouping factor (potentially with interactions)
      group_var_label <- deparse(bar[[3]])
      vars <- trimws(unlist(strsplit(group_var_label, ":")))

      if (length(vars) > 1) {
        flist <- interaction(mf[vars], sep = ":", drop = TRUE)
      } else {
        flist <- as.factor(mf[[vars]])
      }
      group_idx <- as.integer(flist)

      if (length(group_idx) == 0) {
        stop("Grouping factor '", group_var_label, "' resulted in an empty vector.")
      }
      if (any(is.na(group_idx))) {
        stop("NA values found in grouping factor '", group_var_label, "'. Check your data.")
      }

      num_groups <- length(levels(flist))
      num_ranef <- ncol(Z_mf)

      Zt <- matrix(0, nrow = num_groups * num_ranef, ncol = N)
      for (i in 1:N) {
        g <- group_idx[i]
        Zt[((g - 1) * num_ranef + 1):(g * num_ranef), i] <- Z_mf[i, ]
      }

      dat_ranef[[paste0("Zt", suffix(b))]] <- Zt
      dat_ranef[[paste0("group_idx", suffix(b))]] <- group_idx

      num_ranef_list[[b]] <- num_ranef
      num_groups_list[[b]] <- num_groups

      r_names <- colnames(Z_mf)
      r_names[r_names == "(Intercept)"] <- "Int"
      ranef_names_list[[b]] <- r_names
      group_labels_list[[b]] <- group_var_label
    }
  } else {
    mf <- model.frame(formula, data)
    Y <- model.response(mf)
    X <- model.matrix(formula, mf)
    offset <- model.offset(mf)
    num_bars <- 0
    dat_ranef <- list()
    num_ranef_list <- list()
    num_groups_list <- list()
    ranef_names_list <- list()
    group_labels_list <- list()
  }

  has_intercept <- "(Intercept)" %in% colnames(X)
  if (family == "ordered") has_intercept <- FALSE
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  fixed_names <- colnames(X)
  K_tmp <- ncol(X)

  if (is.matrix(Y)) {
    if (ncol(Y) == 2 && family == "binomial") {
      trials <- as.numeric(Y[, 1] + Y[, 2])
      Y <- as.numeric(Y[, 1])
    } else if (ncol(Y) == 1) {
      Y <- as.numeric(Y[, 1])
      trials <- rep(1, length(Y))
    } else stop("Invalid matrix format for the response variable.")
  } else {
    trials <- rep(1, length(Y))
    if (is.factor(Y)) {
      Y <- if (family %in% c("bernoulli", "binomial")) as.numeric(Y) - 1 else as.numeric(Y)
    } else Y <- as.numeric(Y)
  }

  # --- Automatic generation of initial values for SSP ---
  if (regularization == "ssp" && is.null(init) && K_tmp > 0) {
    rhs_mod <- NULL
    suppressMessages(capture.output({
      rhs_prior_obj <- prior_rhs(expected_vars = prior$expected_vars, slab_scale = prior$slab_scale, slab_df = prior$slab_df)
      rhs_mod <- rtmb_glmer(formula, data, family, laplace, prior = rhs_prior_obj, y_range = y_range, init = NULL)
    }))

    jac_target <- if (laplace) "random" else "none"

    ad_setup <- NULL
    suppressMessages(capture.output({
      ad_setup <- rhs_mod$build_ad_obj(laplace = laplace, jacobian_target = jac_target)
    }))

    ad_obj <- ad_setup$ad_obj
    best_obj <- Inf
    best_par <- ad_obj$par

    for (i in 1:4) {
      start_par <- if (i == 1) ad_obj$par else ad_obj$par + rnorm(length(ad_obj$par), 0, 0.5)
      opt <- suppressWarnings(try(nlminb(start = start_par, objective = ad_obj$fn, gradient = ad_obj$gr), silent = TRUE))
      if (!inherits(opt, "try-error") && !is.na(opt$objective) && opt$objective < best_obj) {
        best_obj <- opt$objective
        best_par <- opt$par
      }
    }

    ad_obj$fn(best_par)
    full_par <- if (!is.null(ad_obj$env$last.par.best)) ad_obj$env$last.par.best else ad_obj$env$last.par

    unc_est_list <- unconstrained_vector_to_list(full_par, rhs_mod$par_list)
    con_est_list <- to_constrained(unc_est_list, rhs_mod$par_list)
    if (!is.null(rhs_mod$transform)) {
      tran_res <- rhs_mod$transform(rhs_mod$data, con_est_list)
      con_est_list <- c(con_est_list, tran_res)
    }

    b_est <- as.numeric(con_est_list$b)
    if (length(b_est) == 0 || any(is.na(b_est))) b_est <- rep(0, K_tmp)

    r_init <- ifelse(abs(b_est) > mean(abs(b_est)), 0.9, 0.1)
    tau_init <- pmax(abs(b_est), 0.01)
    beta_raw_init <- b_est / (r_init * tau_init)

    init <- list(beta_raw = beta_raw_init, r = r_init, tau = tau_init)
    if (has_intercept && !is.null(con_est_list$Intercept_c)) {
      init$Intercept_c <- as.numeric(con_est_list$Intercept_c)[1]
    }
  } else if (is.null(init)) {
    tryCatch({
      fixed_formula <- nobars(formula)
      glm_fam <- switch(family, "gaussian" = gaussian(), "student_t" = gaussian(), "lognormal" = gaussian(link="log"), "bernoulli" = binomial(), "binomial" = binomial(), "poisson" = poisson(), "neg_binomial" = poisson(), "gamma" = Gamma(link="log"), gaussian())
      init_fit <- suppressWarnings(glm(fixed_formula, data = data, family = glm_fam))
      init_coef <- unname(coef(init_fit))
      init_coef[is.na(init_coef)] <- 0

      init <- list()
      if (has_intercept) {
        if (prior_type %in% c("flat", "uniform")) {
          init$Intercept <- init_coef[1]
        } else {
          init$Intercept_c <- init_coef[1]
        }
      }

      if (regularization == "none" && K_tmp > 0) {
        if (has_intercept) init$b <- init_coef[-1]
        else init$b <- init_coef
      }
    }, error = function(e) init <- NULL)
  }

  dat <- list(Y = Y, trials = trials, X = X)
  if (has_random) {
    dat <- c(dat, dat_ranef)
  }
  if (!is.null(offset)) dat$offset <- offset
  if (family == "ordered") dat$num_categories <- length(unique(Y))

  # --- Setup AST ---
  setup_exprs <- list()
  setup_exprs[[1]] <- quote(N <- length(Y))
  setup_exprs[[2]] <- quote(K <- ncol(X))
  if (has_random) {
    for (b in 1:num_bars) {
      group_idx_name <- as.name(paste0("group_idx", suffix(b)))
      Zt_name <- as.name(paste0("Zt", suffix(b)))
      Z_mat_name <- as.name(paste0("Z_mat", suffix(b)))
      num_groups_name <- as.name(paste0("num_groups", suffix(b)))
      num_ranef_name <- as.name(paste0("num_ranef", suffix(b)))

      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(num_groups_name) <- max(.(group_idx_name)))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(num_ranef_name) <- nrow(.(Zt_name)) / .(num_groups_name))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(Z_mat_name) <- matrix(0, nrow = N, ncol = .(num_ranef_name)))

      loop_body <- bquote({
        g <- .(group_idx_name)[i]
        .(Z_mat_name)[i, ] <- .(Zt_name)[((g - 1) * .(num_ranef_name) + 1):(g * .(num_ranef_name)), i]
      })
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(for (i in 1:N) .(loop_body))
    }
  }

  if (use_weak_info) {
    if (family %in% c("bernoulli", "binomial", "ordered")) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- pi / sqrt(3))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(alpha_prior_sd <- base_scale * .(prior$max_beta))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- 0)
    } else if (family %in% c("poisson", "neg_binomial", "gamma")) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- 1.0)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- base_scale)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- 0)
    } else {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(half_d_y <- diff(y_range) / 2)
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(base_scale <- half_d_y * .(prior$sd_ratio))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- half_d_y)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mean(y_range))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate <- 1.0 / base_scale)
    }

    setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_rate <- 1.0 / base_scale)

    if (K_tmp > 0) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_mean <- apply(X, 2, mean))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_sd <- apply(X, 2, sd))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(beta_prior_sd <- .(prior$max_beta) * base_scale / X_sd)
      if (has_intercept && !prior_type %in% c("flat", "uniform")) {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(X_c <- X - rep(1, N) %*% t(X_mean))
      }

      if (regularization == "rhs") {
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(p0 <- min(.(prior$expected_vars), K - 1))
        setup_exprs[[length(setup_exprs) + 1]] <- quote(if (p0 < 1) p0 <- 1)
        setup_exprs[[length(setup_exprs) + 1]] <- quote(tau0 <- (p0 / (K - p0)) * (base_scale / sqrt(N)))
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(half_slab_df <- .(prior$slab_df) / 2.0)
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(half_slab_scale2 <- (.(prior$slab_df) * .(prior$slab_scale)^2) / 2.0)
      } else if (regularization == "ssp") {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_scale <- base_scale / X_sd)
      }
    }
  } else {
    # Minimal calculations for fixed prior distributions
    # No centering or scaling for flat/uniform priors
    NULL
  }
  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(dat)
  eval(setup_ast, tmp_env)
  N <- tmp_env$N; K <- tmp_env$K

  # --- Parameters AST ---
  param_exprs <- list()
  if (has_intercept) {
    if (prior_type %in% c("flat", "uniform")) {
      param_exprs[[length(param_exprs) + 1]] <- quote(Intercept <- Dim(1))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- quote(Intercept_c <- Dim(1))
    }
  }
  if (K > 0) {
    if (regularization == "none") {
      param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(K))
    } else if (regularization == "rhs") {
      param_exprs[[length(param_exprs) + 1]] <- quote(z <- Dim(K))
      param_exprs[[length(param_exprs) + 1]] <- quote(lambda <- Dim(K, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(w_lambda <- Dim(K, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(tau_hs <- Dim(1, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(w_tau <- Dim(1, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(c2 <- Dim(1, lower = 0))
    } else if (regularization == "ssp") {
      param_exprs[[length(param_exprs) + 1]] <- quote(beta_raw <- Dim(K))
      param_exprs[[length(param_exprs) + 1]] <- quote(r <- Dim(K, lower = 0.001, upper = 0.999))
      param_exprs[[length(param_exprs) + 1]] <- quote(tau <- Dim(K, lower = 0))
    }
  }

  if (has_random) {
    for (b in 1:num_bars) {
      num_ranef_val <- num_ranef_list[[b]]
      num_groups_val <- num_groups_list[[b]]
      sd_name <- as.name(paste0("sd", suffix(b)))
      r_re_name <- as.name(paste0("r_re", suffix(b)))
      CF_corr_name <- as.name(paste0("CF_corr", suffix(b)))

      if (num_ranef_val == 1) {
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(sd_name) <- Dim(.(num_ranef_val), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(r_re_name) <- Dim(.(num_groups_val), random = TRUE))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(sd_name) <- Dim(.(num_ranef_val), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(CF_corr_name) <- Dim(c(.(num_ranef_val), .(num_ranef_val)), type = "CF_corr"))
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(r_re_name) <- Dim(c(.(num_groups_val), .(num_ranef_val)), random = TRUE))
      }
    }
  }

  if (family %in% c("gaussian", "lognormal", "student_t")) param_exprs[[length(param_exprs) + 1]] <- quote(sigma <- Dim(1, lower = 0))
  if (family == "student_t") param_exprs[[length(param_exprs) + 1]] <- quote(nu <- Dim(1, lower = 2))
  if (family == "gamma") param_exprs[[length(param_exprs) + 1]] <- quote(shape <- Dim(1, lower = 0))
  if (family == "neg_binomial") param_exprs[[length(param_exprs) + 1]] <- quote(phi <- Dim(1, lower = 0))
  if (family == "ordered") param_exprs[[length(param_exprs) + 1]] <- quote(cutpoints <- Dim(num_categories - 1, type = "ordered"))
  param_ast <- as.call(c(list(as.name("{")), param_exprs))

  # --- Transform AST ---
  tran_exprs <- list()
  if (K > 0) {
    if (regularization == "rhs") {
      tran_exprs[[length(tran_exprs) + 1]] <- quote(lambda_sq <- lambda^2)
      tran_exprs[[length(tran_exprs) + 1]] <- quote(tau_sq <- tau_hs^2)
      tran_exprs[[length(tran_exprs) + 1]] <- quote(lambda_tilde <- sqrt((c2 * lambda_sq) / (c2 + tau_sq * lambda_sq)))
      tran_exprs[[length(tran_exprs) + 1]] <- quote(b <- z * lambda_tilde * tau_hs)
    } else if (regularization == "ssp") {
      tran_exprs[[length(tran_exprs) + 1]] <- quote(b <- beta_raw * r * tau)
    }
  }
  if (has_intercept) {
    if (!prior_type %in% c("flat", "uniform")) {
      if (K > 0) tran_exprs[[length(tran_exprs) + 1]] <- quote(Intercept <- Intercept_c - sum(X_mean * b))
      else tran_exprs[[length(tran_exprs) + 1]] <- quote(Intercept <- Intercept_c)
    }
  }
  if (has_random) {
    for (b in 1:num_bars) {
      if (num_ranef_list[[b]] > 1) {
        CF_corr_name <- as.name(paste0("CF_corr", suffix(b)))
        corr_name <- as.name(paste0("corr", suffix(b)))
        tran_exprs[[length(tran_exprs) + 1]] <- bquote(.(corr_name) <- .(CF_corr_name) %*% t(.(CF_corr_name)))
      }
    }
  }
  tran_ast <- if (length(tran_exprs) > 0) as.call(c(list(as.name("{")), tran_exprs)) else NULL

  # --- Model AST ---
  transform_exprs <- list()
  if (has_intercept) {
    if (prior_type %in% c("flat", "uniform")) {
      if (K > 0) transform_exprs[[1]] <- quote(eta <- Intercept + X %*% b)
      else transform_exprs[[1]] <- quote(eta <- rep(Intercept, N))
    } else {
      if (K > 0) transform_exprs[[1]] <- quote(eta <- Intercept_c + X_c %*% b)
      else transform_exprs[[1]] <- quote(eta <- rep(Intercept_c, N))
    }
  } else {
    if (K > 0) transform_exprs[[1]] <- quote(eta <- X %*% b)
    else transform_exprs[[1]] <- quote(eta <- rep(0, N))
  }

  if (has_random) {
    for (b in 1:num_bars) {
      Z_mat_name <- as.name(paste0("Z_mat", suffix(b)))
      r_re_name <- as.name(paste0("r_re", suffix(b)))
      group_idx_name <- as.name(paste0("group_idx", suffix(b)))
      sd_name <- as.name(paste0("sd", suffix(b)))

      if (num_ranef_list[[b]] > 1) {
        transform_exprs[[length(transform_exprs) + 1]] <- bquote(for (i in 1:N) eta[i] <- eta[i] + sum(.(Z_mat_name)[i, ] * .(r_re_name)[.(group_idx_name)[i], ] * .(sd_name)))
      } else {
        transform_exprs[[length(transform_exprs) + 1]] <- bquote(eta <- eta + .(Z_mat_name)[,1] * .(r_re_name)[.(group_idx_name)] * .(sd_name))
      }
    }
  }
  if (!is.null(offset)) transform_exprs[[length(transform_exprs) + 1]] <- quote(eta <- eta + offset)

  ll_data_exprs <- list()
  ll_data_exprs[[1]] <- switch(family,
                               "gaussian" = quote(Y ~ normal(eta, sigma)),
                               "lognormal" = quote(Y ~ lognormal(eta, sigma)),
                               "student_t" = quote(Y ~ student_t(nu, eta, sigma)),
                               "gamma" = quote(Y ~ gamma(shape, shape / exp(eta))),
                               "bernoulli" = quote(Y ~ bernoulli_logit(eta)),
                               "binomial" = quote(Y ~ binomial_logit(trials, eta)),
                               "poisson" = quote(Y ~ poisson(exp(eta))),
                               "neg_binomial" = quote(Y ~ neg_binomial_2(exp(eta), phi)),
                               "ordered" = quote(Y ~ ordered_logistic(eta, cutpoints))
  )

  ll_random_exprs <- list()
  if (has_random) {
    for (b in 1:num_bars) {
      r_re_name <- as.name(paste0("r_re", suffix(b)))
      if (num_ranef_list[[b]] > 1) {
        CF_corr_name <- as.name(paste0("CF_corr", suffix(b)))
        num_ranef_name <- as.name(paste0("num_ranef", suffix(b)))
        num_groups_name <- as.name(paste0("num_groups", suffix(b)))

        ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(.(CF_corr_name) ~ lkj_CF_corr(.(prior$lkj_eta)))
        ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(for (j in 1:.(num_groups_name)) .(r_re_name)[j, ] ~ multi_normal_CF(rep(0, .(num_ranef_name)), rep(1, .(num_ranef_name)), .(CF_corr_name)))
      } else {
        ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(.(r_re_name) ~ normal(0, 1))
      }
    }
  }

  prior_exprs <- list()
  if (prior_type %in% c("weak", "rhs", "ssp")) {
    if (family %in% c("gaussian", "lognormal", "student_t")) {
      if (use_weak_info && prior_type == "weak") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(sigma ~ exponential(sigma_rate))
      } else if (!is.null(prior$sigma_rate)) {
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
      }
    }

    if (family == "student_t" && !is.null(prior$nu_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(nu ~ exponential(.(prior$nu_rate)))
    if (family == "gamma" && !is.null(prior$shape_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(shape ~ exponential(.(prior$shape_rate)))
    if (family == "neg_binomial" && !is.null(prior$phi_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(phi ~ exponential(.(prior$phi_rate)))
    if (family == "ordered" && !is.null(prior$cutpoint_sd)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))

    if (has_intercept) {
      if (use_weak_info && prior_type == "weak") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(Intercept_c ~ normal(mid_y, alpha_prior_sd))
      } else if (!is.null(prior$Intercept_sd)) {
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(Intercept_c ~ normal(0, .(prior$Intercept_sd)))
      }
    }

    if (K > 0) {
      if (regularization == "none") {
        if (use_weak_info && prior_type == "weak") {
          prior_exprs[[length(prior_exprs) + 1]] <- quote(b ~ normal(0, beta_prior_sd))
        } else if (!is.null(prior$b_sd)) {
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
        }
      } else if (regularization == "rhs") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(z ~ normal(0, 1))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(w_lambda ~ gamma(0.5, 0.5))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(lambda ~ half_normal(1 / sqrt(w_lambda)))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(w_tau ~ gamma(0.5, 0.5))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(tau_hs ~ half_normal(tau0 / sqrt(w_tau)))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(c2 ~ inverse_gamma(half_slab_df, half_slab_scale2))
      } else if (regularization == "ssp") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(beta_raw ~ laplace(0, 0.5))
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(mu_r <- log(.(prior$ssp_ratio) / (1 - .(prior$ssp_ratio))))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(r ~ logit_normal(mu_r, 3))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(tau ~ exponential(1 / tau_scale))
      }
    }

    if (has_random) {
      for (b in 1:num_bars) {
        sd_name <- as.name(paste0("sd", suffix(b)))
        if (use_weak_info && prior_type == "weak") {
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(.(sd_name) ~ exponential(tau_rate))
        } else if (!is.null(prior$sd_rate)) {
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(.(sd_name) ~ exponential(.(prior$sd_rate)))
        }
      }
    }
  } else if (prior_type == "uniform") {
    if (!is.null(prior$sigma_rate)) {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    }
    if (has_intercept && !is.null(prior$Intercept_sd)) {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(Intercept ~ normal(0, .(prior$Intercept_sd)))
    }
    if (K > 0 && !is.null(prior$b_sd)) {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
    }
    if (has_random && !is.null(prior$tau_rate)) {
      for (b in 1:num_bars) {
        sd_name <- as.name(paste0("sd", suffix(b)))
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(.(sd_name) ~ exponential(.(prior$tau_rate)))
      }
    }
  }

  model_exprs <- list()
  model_exprs[[length(model_exprs) + 1]] <- "# Transform"
  model_exprs <- c(model_exprs, transform_exprs)
  model_exprs[[length(model_exprs) + 1]] <- "# Likelihood (Data)"
  model_exprs <- c(model_exprs, ll_data_exprs)
  if (has_random) {
    model_exprs[[length(model_exprs) + 1]] <- "# Likelihood (Random)"
    model_exprs <- c(model_exprs, ll_random_exprs)
  }
  model_exprs[[length(model_exprs) + 1]] <- "# Priors"
  model_exprs <- c(model_exprs, prior_exprs)

  model_ast <- as.call(c(list(as.name("{")), model_exprs))

  code_obj <- list(setup = setup_ast, parameters = param_ast)
  if (!is.null(tran_ast)) code_obj$transform <- tran_ast
  code_obj$model <- model_ast

  par_names_list <- list()
  if (K > 0) {
    par_names_list$b <- fixed_names
    if (regularization == "rhs") {
      par_names_list$z <- fixed_names; par_names_list$lambda <- fixed_names; par_names_list$w_lambda <- fixed_names
    } else if (regularization == "ssp") {
      par_names_list$beta_raw <- fixed_names; par_names_list$r <- fixed_names; par_names_list$tau <- fixed_names
    }
  }
  if (has_random) {
    for (b in 1:num_bars) {
      r_names <- ranef_names_list[[b]]
      if (num_bars > 1) {
        r_names <- paste0(r_names, ":", group_labels_list[[b]])
      }
      par_names_list[[paste0("sd", suffix(b))]] <- r_names
      if (num_ranef_list[[b]] > 1) par_names_list[[paste0("corr", suffix(b))]] <- r_names
    }
  }

  view_vars <- c()
  if (has_intercept) view_vars <- c("Intercept")
  if (K > 0) view_vars <- c(view_vars, "b")
  view_vars <- c(view_vars, "sigma")
  if (has_random) {
    for (b in 1:num_bars) {
      view_vars <- c(view_vars, paste0("sd", suffix(b)))
      if (num_ranef_list[[b]] > 1) view_vars <- c(view_vars, paste0("corr", suffix(b)))
    }
  }

  ordered_data <- env_to_ordered_list(tmp_env, dat, setup_ast)
  obj <- rtmb_model(data = ordered_data, code = code_obj, par_names = par_names_list, init = init, view = view_vars)
  obj$formula <- formula
  obj$raw_data <- data
  obj$family <- family

  if (!is.null(null)) {
    obj <- obj$null_model(pars = null)
  }

  return(obj)
}

#' RTMB-based Linear Mixed Model (LMM) wrapper function
#'
#' @param formula Formula
#' @param data Data frame
#' @param laplace Logical; whether to marginalize random effects using Laplace approximation
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Default is NULL (flat prior).
#' @param y_range Theoretical minimum and maximum values of the response variable
#' @param init Initial values
#' @param null Null model parameters
#' @return RTMB_Model object
#' @export
#' @example inst/examples/ex_lm.R
rtmb_lmer <- function(formula, data, laplace = TRUE,
                      prior = prior_uniform(),
                      y_range = NULL,
                      init = NULL, null = NULL) {
  rtmb_glmer(formula = formula, data = data, family = "gaussian",
             laplace = laplace,
             prior = prior,
             y_range = y_range,
             init = init,
             null = null)
}

#' RTMB-based GLM wrapper function (no random effects)
#'
#' @param formula Formula
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Use prior_weak(), prior_rhs(), or prior_ssp(). Default is NULL (flat prior).
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param init List of initial values
#' @param null Character string specifying the target parameter for the null model.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glm <- function(formula, data, family = "gaussian",
                     prior = prior_uniform(),
                     y_range = NULL,
                     init = NULL, null = NULL) {
  rtmb_glmer(
    formula = formula,
    data = data,
    family = family,
    laplace = FALSE,
    prior = prior,
    y_range = y_range,
    init = init,
    null = null
  )
}

#' RTMB-based Linear Regression wrapper function
#'
#' @param formula Formula (e.g., Y ~ X1 + X2)
#' @param data Data frame
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Use prior_weak(), prior_rhs(), or prior_ssp(). Default is NULL (flat prior).
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param init List of initial values.
#' @param null Character string specifying the target parameter for the null model.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_lm <- function(formula, data,
                    prior = prior_uniform(),
                    y_range = NULL,
                    init = NULL, null = NULL) {

  rtmb_glm(
    formula = formula,
    data = data,
    family = "gaussian",
    prior = prior,
    y_range = y_range,
    init = init,
    null = null
  )
}



#' RTMB-based Factor Analysis Wrapper
#'
#' @description
#' Performs exploratory factor analysis (EFA) using RTMB. It supports standard
#' rotation methods (e.g., varimax, promax) as well as regularized factor analysis
#' using a Spike-and-Slab Prior (SSP) for estimating sparse loading matrices.
#'
#' @param data Observation data frame or matrix (N x J).
#' @param nfactors Number of factors (K).
#' @param rotate String specifying the rotation method (e.g., "varimax", "promax", "ssp"). If NULL, no rotation is applied. Specifying "ssp" performs regularized factor analysis.
#' @param score Logical; if TRUE, factor scores are calculated in the generate block (default is FALSE).
#' @param prior List of hyperparameters for prior distributions. `ssp_ratio` represents the proportion of non-zero loadings per factor when "ssp" is specified.
#' @param init List of initial values. If not provided, initial values are automatically generated based on PCA or the psych package.
#' @example inst/examples/ex_fa.R
#' @export
rtmb_fa <- function(data, nfactors = 1, rotate = NULL, score = FALSE,
                    prior = list(mean_sd = 10, loadings_sd = 1, sd_rate = 10, ssp_ratio = 0.25),
                    init = NULL) {

  Y <- as.matrix(data)
  K <- nfactors

  default_prior <- list(mean_sd = 10, loadings_sd = 1, sd_rate = 10, ssp_ratio = 0.25)
  if (!is.null(prior)) {
    prior <- modifyList(default_prior, prior)
  } else {
    prior <- default_prior
  }

  var_names <- colnames(data)
  if (is.null(var_names)) var_names <- paste0("V", 1:ncol(Y))

  if (K >= ncol(Y)) stop("The number of factors (K) must be less than the number of observed variables (J).")

  # Determine if SSP model is used
  is_ssp <- !is.null(rotate) && rotate == "ssp"

  # Common setup block
  setup_ast <- quote({
    N <- nrow(Y)
    J <- ncol(Y)

    # Check for missing values and calculate sufficient statistics accordingly
    if (any(is.na(Y))) {
      warning("Missing values (NAs) detected in the data. Using pairwise complete observations for sufficient statistics. Results may be approximate.")
      Y_bar <- apply(Y, 2, mean, na.rm = TRUE)
      S_Y <- cov(Y, use = "pairwise.complete.obs") * (N - 1)
    } else {
      Y_bar <- apply(Y, 2, mean)
      S_Y <- cov(Y) * (N - 1)
    }
  })

  if (is_ssp) {
    if (is.null(prior$ssp_ratio)) prior$ssp_ratio <- 0.25

    dat_fa <- list(
      Y = Y, K_factors = K,
      prior_mean_sd = prior$mean_sd, prior_sd_rate = prior$sd_rate,
      ssp_ratio = prior$ssp_ratio
    )

    tmp_env <- list2env(dat_fa)
    eval(setup_ast, tmp_env)
    N <- tmp_env$N; J <- tmp_env$J; Y_bar <- tmp_env$Y_bar; S_Y <- tmp_env$S_Y

    if (score) dat_fa$Y <- Y

    # Automatic generation of initial values for SSP (using PCA and promax from Base R)
    if (is.null(init)) {
      tryCatch({
        eig <- eigen(S_Y / (N - 1))
        L_pca <- eig$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(eig$values[1:K]), nrow = K, ncol = K)
        var_Y <- diag(S_Y / (N - 1))

        if (K > 1) {
          pm <- stats::promax(L_pca, m = 4)
          init_Lambda <- unclass(pm$loadings)
          T_mat <- pm$rotmat
          Phi <- cov2cor(solve(t(T_mat) %*% T_mat))
          init_CF_Omega <- t(chol(Phi))
        } else {
          init_Lambda <- L_pca
          Phi <- matrix(1, nrow = 1, ncol = 1)
          init_CF_Omega <- matrix(1, nrow = 1, ncol = 1)
        }

        h2_raw <- rowSums(init_Lambda * (init_Lambda %*% Phi))
        init_sd <- sqrt(pmax(var_Y - h2_raw, 0.01 * var_Y))

        L_std <- init_Lambda / sqrt(var_Y)
        init_r <- ifelse(abs(L_std) > 0.2, 0.9, 0.1)

        init <- list(
          mean = Y_bar,
          Lambda_star = init_Lambda,
          sd = init_sd,
          CF_Omega = init_CF_Omega,
          r = init_r,
          tau = matrix(1.0, nrow = J, ncol = K)
        )
      }, error = function(e) {
        stop("Failed to automatically generate initial values. Please provide initial values manually via the 'init' argument. Details: ", conditionMessage(e))
      })
    }

    param_ast <- quote({
      mean = Dim(dim = J)
      Lambda_star = Dim(c(J, K_factors))
      r = Dim(c(J, K_factors), lower = 0.001, upper = 0.999)
      tau = Dim(c(J, K_factors), lower = 0)
      sd = Dim(J, lower = 0)
      CF_Omega = Dim(c(K_factors, K_factors), type = "CF_corr")
    })

    tran_ast <- quote({
      L_raw <- Lambda_star * r * tau
      h2 <- rowSums(L_raw * (L_raw %*% CF_Omega))
      var_Y <- h2 + sd^2
      sd_Y <- sqrt(var_Y)
      L <- L_raw / sd_Y
      fa_cor <- CF_Omega %*% t(CF_Omega)
    })

    model_ast <- quote({
      S_Y ~ sufficient_multi_normal_fa(N, Y_bar, mean, sd, L_raw %*% CF_Omega)
      mean ~ normal(0, prior_mean_sd)
      sd ~ exponential(prior_sd_rate)
      CF_Omega ~ lkj_CF_corr(1)

      mu_r <- log(ssp_ratio / (1 - ssp_ratio))
      r ~ logit_normal(mu_r, 3)
      tau ~ exponential(1)
      Lambda_star ~ laplace(0, 1)
    })

    base_gq <- quote({
      Sigma <- L_raw %*% fa_cor %*% t(L_raw) + diag(sd^2)
      var_total <- diag(Sigma)
      var_common <- rowSums(L_raw * (L_raw %*% CF_Omega))
      communality <- var_common / var_total
      out <- list(communality = communality)
    })

    score_expr <- if (score) {
      quote({
        Y_c <- matrix(0, nrow = N, ncol = J)
        for (i in 1:N) Y_c[i, ] <- Y[i, ] - mean
        out$score <- Y_c %*% solve(Sigma, L_raw %*% fa_cor)
      })
    } else quote({})

    ret_expr <- quote({ return(out) })
    gq_ast <- as.call(c(list(as.name("{")), as.list(base_gq)[-1], as.list(score_expr)[-1], as.list(ret_expr)[-1]))

    code_args <- list(setup = setup_ast, parameters = param_ast, transform = tran_ast, model = model_ast, generate = gq_ast)
    code_obj <- eval(as.call(c(list(as.name("rtmb_code")), code_args)))

    p_names <- list(
      mean = var_names,
      Lambda_star = var_names,
      r = var_names,
      tau = var_names,
      sd = var_names,
      CF_Omega = paste0("Factor", 1:K),
      L_raw = var_names,
      L = var_names,
      fa_cor = paste0("Factor", 1:K),
      communality = var_names
    )
    if (score) {
      ind_names <- rownames(data)
      if (is.null(ind_names)) ind_names <- paste0("Id", 1:N)
      p_names[["score"]] <- list(ind_names, paste0("Factor", 1:K))
    }

    target_view <- c("L", "sd", "fa_cor")

    obj <- rtmb_model(
      data = dat_fa,
      code = code_obj,
      par_names = p_names,
      init = init,
      view = target_view
    )

    return(obj)

  } else {
    # --- Existing rotation logic (varimax, promax, etc.) ---
    dat_fa <- list(
      Y = Y, K_factors = K,
      prior_mean_sd = prior$mean_sd, prior_loadings_sd = prior$loadings_sd, prior_sd_rate = prior$sd_rate
    )

    tmp_env <- list2env(dat_fa)
    eval(setup_ast, tmp_env)
    N <- tmp_env$N; J <- tmp_env$J; Y_bar <- tmp_env$Y_bar; S_Y <- tmp_env$S_Y

    if (score) dat_fa$Y <- Y

    if (is.null(init)) {
      tryCatch({
        eig <- eigen(S_Y / (N - 1))
        L_pca <- eig$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(eig$values[1:K]), nrow=K, ncol=K)
        L_top <- L_pca[1:K, , drop = FALSE]
        qr_res <- qr(t(L_top))
        Q_mat <- qr.Q(qr_res)
        init_loadings <- L_pca %*% Q_mat
        init_sd <- pmax(diag(S_Y / (N - 1)) - rowSums(init_loadings^2), 0.01)^0.5

        for (j in 1:J) {
          for (k in 1:K) {
            if (j < k) init_loadings[j, k] <- 0
          }
        }
        init <- list(mean = Y_bar, L_raw = init_loadings, sd = init_sd)
      }, error = function(e) { init <- NULL })
    }

    param_ast <- quote({
      mean <- Dim(dim = J)
      L_raw <- Dim(dim = c(J, K_factors), type = "lower_tri")
      sd <- Dim(dim = J, lower = 0)
    })

    tran_ast <- quote({
      h2 <- rowSums(L_raw^2)
      var_Y <- h2 + sd^2
      sd_Y <- sqrt(var_Y)
      L <- L_raw / sd_Y
    })

    model_ast <- quote({
      S_Y ~ sufficient_multi_normal_fa(N, Y_bar, mean, sd, L_raw)
      mean ~ normal(0, prior_mean_sd)
      sd ~ exponential(prior_sd_rate)
      L_raw ~ lower_tri_normal(0, prior_loadings_sd)
    })

    base_gq <- quote({
      Sigma <- L_raw %*% t(L_raw) + diag(sd^2)
      var_total <- diag(Sigma)
      var_common <- rowSums(L_raw^2)
      communality <- var_common / var_total
      out <- list(communality = communality)
    })

    if (!is.null(rotate)) {
      rot_loadings_name <- paste0("L_", rotate)

      if (exists(rotate, mode = "function")) {
        rot_fn <- match.fun(rotate)
        fn_call <- as.name(rotate)
      } else if (requireNamespace("GPArotation", quietly = TRUE) &&
                 exists(rotate, where = asNamespace("GPArotation"), mode = "function")) {
        rot_fn <- getFromNamespace(rotate, "GPArotation")
        fn_call <- call("::", as.name("GPArotation"), as.name(rotate))
      } else {
        stop("Rotation function not found: ", rotate, ". If this is from GPArotation, please install it using install.packages('GPArotation').")
      }

      dummy_L <- matrix(rnorm(J * K), J, K)
      test_rot <- rot_fn(dummy_L)

      is_matrix_rot <- is.matrix(test_rot)
      has_phi <- !is_matrix_rot && !is.null(test_rot$Phi)

      if (is_matrix_rot) {
        rot_expr <- bquote({
          rot_obj <- .(fn_call)(L)
          rot_mat <- unclass(rot_obj)
          out[[.(rot_loadings_name)]] <- rot_mat
        })
        score_expr <- if (score) bquote({
          Y_c <- matrix(0, nrow = N, ncol = J)
          for (i in 1:N) Y_c[i, ] <- Y[i, ] - mean
          rot_raw <- unclass(.(fn_call)(L_raw))
          if (!is.matrix(rot_raw)) rot_raw <- unclass(rot_raw$loadings)
          out$score <- Y_c %*% solve(Sigma, rot_raw)
        }) else quote({})
      } else {
        if (has_phi) {
          rot_expr <- bquote({
            rot_obj <- .(fn_call)(L)
            rot_mat <- unclass(rot_obj$loadings)
            out$fa_cor <- rot_obj$Phi
            out[[.(rot_loadings_name)]] <- rot_mat
          })
          score_expr <- if (score) bquote({
            Y_c <- matrix(0, nrow = N, ncol = J)
            for (i in 1:N) Y_c[i, ] <- Y[i, ] - mean
            rot_raw_obj <- .(fn_call)(L_raw)
            rot_raw_mat <- unclass(rot_raw_obj$loadings)
            out$score <- Y_c %*% solve(Sigma, rot_raw_mat %*% rot_raw_obj$Phi)
          }) else quote({})
        } else {
          rot_expr <- bquote({
            rot_obj <- .(fn_call)(L)
            rot_mat <- unclass(rot_obj$loadings)
            out[[.(rot_loadings_name)]] <- rot_mat
          })
          score_expr <- if (score) bquote({
            Y_c <- matrix(0, nrow = N, ncol = J)
            for (i in 1:N) Y_c[i, ] <- Y[i, ] - mean
            rot_raw <- unclass(.(fn_call)(L_raw)$loadings)
            out$score <- Y_c %*% solve(Sigma, rot_raw)
          }) else quote({})
        }
      }
    } else {
      has_phi <- FALSE
      rot_expr <- quote({})
      score_expr <- if (score) {
        quote({
          Y_c <- matrix(0, nrow = N, ncol = J)
          for (i in 1:N) Y_c[i, ] <- Y[i, ] - mean
          out$score <- Y_c %*% solve(Sigma, L_raw)
        })
      } else quote({})
    }

    ret_expr <- quote({ return(out) })

    gq_ast <- as.call(c(list(as.name("{")), as.list(base_gq)[-1], as.list(rot_expr)[-1], as.list(score_expr)[-1], as.list(ret_expr)[-1]))

    code_args <- list(setup = setup_ast, parameters = param_ast, transform = tran_ast, model = model_ast, generate = gq_ast)
    code_obj <- eval(as.call(c(list(as.name("rtmb_code")), code_args)))

    p_names <- list(
      mean = var_names,
      L_raw = var_names,
      sd = var_names,
      L = var_names,
      communality = var_names
    )
    if (!is.null(rotate)) {
      p_names[[paste0("L_", rotate)]] <- var_names
      if (has_phi) p_names[["fa_cor"]] <- paste0("Factor", 1:K)
    }
    if (score) {
      ind_names <- rownames(data)
      if (is.null(ind_names)) ind_names <- paste0("Id", 1:N)
      p_names[["score"]] <- list(ind_names, paste0("Factor", 1:K))
    }

    target_view <- if (!is.null(rotate)) c(paste0("L_", rotate), "sd", "fa_cor") else c("L", "sd", "fa_cor")

    obj <- rtmb_model(
      data = dat_fa,
      code = code_obj,
      par_names = p_names,
      init = init,
      view = target_view
    )

    return(obj)
  }
}

#' RTMB-based IRT (Item Response Theory) Wrapper
#'
#' @description
#' Performs Item Response Theory modeling (1PL, 2PL, 3PL, and Graded Response Model) using RTMB.
#' Missing values (NA) in the data are automatically removed internally for efficient computation.
#'
#' @param data A data frame or matrix of item responses (N persons x J items).
#' @param model Character string for the model type: "1PL", "2PL", or "3PL".
#' @param type Character string for the data type: "binary" or "ordered".
#' @param prior List of hyperparameters for prior distributions.
#' @param init List of initial values.
#' @example inst/examples/ex_irt.R
#' @export
rtmb_irt <- function(data, model = c("2PL", "1PL", "3PL"), type = c("binary", "ordered"),
                     prior = list(a_log_mean = 0, a_log_sd = 0.5, b_mean = 0, b_sd = 2.5, c_alpha = 1, c_beta = 4, theta_sd = 1),
                     init = NULL) {

  model <- match.arg(model)
  type <- match.arg(type)

  if (model == "3PL" && type == "ordered") {
    stop("The 3PL model only supports 'binary' data.")
  }

  Y <- as.matrix(data)

  # Get row and column names
  item_names <- colnames(Y)
  if (is.null(item_names)) item_names <- paste0("Item", 1:ncol(Y))
  person_names <- rownames(Y)
  if (is.null(person_names)) person_names <- paste0("Person", 1:nrow(Y))

  # Exclude NA and convert to long format (for speedup and reducing unnecessary computation graph)
  obs_data <- which(!is.na(Y), arr.ind = TRUE)
  person_idx <- as.integer(obs_data[, "row"])
  item_idx <- as.integer(obs_data[, "col"])
  Y_obs <- Y[obs_data]

  if (type == "ordered") {
    # Ordered scale assumes categories start from 1 (specification of ordered_logistic)
    min_y <- min(Y_obs)
    if (min_y == 0) {
      Y_obs <- Y_obs + 1
      message("Since the minimum value of the ordered data was 0, it was automatically converted internally to a 1-based index.")
    }
  }

  default_prior <- list(a_log_mean = 0, a_log_sd = 0.5, b_mean = 0, b_sd = 2.5, c_alpha = 1, c_beta = 4, theta_sd = 1)
  if (!is.null(prior)) prior <- modifyList(default_prior, prior) else prior <- default_prior

  # Remove preprocessing from the wrapper function side and pass only the raw matrix Y
  dat <- list(
    Y = Y,
    prior_a_log_mean = prior$a_log_mean,
    prior_a_log_sd = prior$a_log_sd,
    prior_b_mean = prior$b_mean,
    prior_b_sd = prior$b_sd,
    prior_c_alpha = prior$c_alpha,
    prior_c_beta = prior$c_beta,
    prior_theta_sd = prior$theta_sd
  )

  # --- Construction of Setup Block ---
  setup_exprs <- list()
  setup_exprs[[1]] <- quote(obs_data <- which(!is.na(Y), arr.ind = TRUE))
  setup_exprs[[2]] <- quote(person_idx <- as.integer(obs_data[, "row"]))
  setup_exprs[[3]] <- quote(item_idx <- as.integer(obs_data[, "col"]))
  setup_exprs[[4]] <- quote(Y_obs <- Y[obs_data])
  setup_exprs[[5]] <- quote(N_persons <- nrow(Y)) # or max(person_idx)
  setup_exprs[[6]] <- quote(N_items <- ncol(Y))   # or max(item_idx)
  setup_exprs[[7]] <- quote(N_obs <- length(Y_obs))

  if (type == "ordered") {
    # The process of correcting 0-based ordered scale to 1-based is also done in setup
    setup_exprs[[8]] <- quote(if (min(Y_obs) == 0) Y_obs <- Y_obs + 1)
    setup_exprs[[9]] <- quote(K_cat <- max(Y_obs))
  }

  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(dat)
  eval(setup_ast, tmp_env)

  # --- Construction of Parameters Block ---
  param_exprs <- list()
  param_exprs[[length(param_exprs) + 1]] <- quote(theta <- Dim(N_persons, random = TRUE))

  if (type == "binary") {
    param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(N_items))
  } else {
    param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(c(N_items, K_cat - 1), type = "ordered"))
  }

  if (model %in% c("2PL", "3PL")) {
    param_exprs[[length(param_exprs) + 1]] <- quote(a <- Dim(N_items, lower = 0))
  }
  if (model == "3PL") {
    param_exprs[[length(param_exprs) + 1]] <- quote(c <- Dim(N_items, lower = 0, upper = 1))
  }
  param_ast <- as.call(c(list(as.name("{")), param_exprs))

  # --- Construction of Model Block ---
  loop_body <- list()
  loop_body[[1]] <- quote(p <- person_idx[i])
  loop_body[[2]] <- quote(j <- item_idx[i])
  loop_body[[3]] <- quote(y <- Y_obs[i])

  if (type == "binary") {
    if (model == "1PL") {
      loop_body[[4]] <- quote(eta <- theta[p] - b[j])
    } else {
      loop_body[[4]] <- quote(eta <- a[j] * (theta[p] - b[j]))
    }

    if (model == "3PL") {
      loop_body[[5]] <- quote(prob <- c[j] + (1 - c[j]) * plogis(eta))
      loop_body[[6]] <- quote(y ~ bernoulli(prob))
    } else {
      loop_body[[5]] <- quote(y ~ bernoulli_logit(eta))
    }
  } else if (type == "ordered") {
    if (model == "1PL") {
      loop_body[[4]] <- quote(eta <- theta[p])
    } else {
      loop_body[[4]] <- quote(eta <- a[j] * theta[p])
    }
    loop_body[[5]] <- quote(y ~ ordered_logistic(eta, b[j, ]))
  }

  # Construct the entire loop
  loop_ast <- quote(for (i in 1:N_obs) {})
  loop_ast[[4]] <- as.call(c(list(as.name("{")), loop_body))

  model_exprs <- list()
  model_exprs[[length(model_exprs) + 1]] <- "# Likelihood"
  model_exprs[[length(model_exprs) + 1]] <- loop_ast

  model_exprs[[length(model_exprs) + 1]] <- "# Priors"
  if (model %in% c("2PL", "3PL")) {
    model_exprs[[length(model_exprs) + 1]] <- quote(a ~ lognormal(prior_a_log_mean, prior_a_log_sd))
  }
  if (type == "binary") {
    model_exprs[[length(model_exprs) + 1]] <- quote(b ~ normal(prior_b_mean, prior_b_sd))
  } else {
    model_exprs[[length(model_exprs) + 1]] <- quote(for (j in 1:N_items) b[j, ] ~ normal(prior_b_mean, prior_b_sd))
  }
  if (model == "3PL") {
    model_exprs[[length(model_exprs) + 1]] <- quote(c ~ beta(prior_c_alpha, prior_c_beta))
  }
  model_exprs[[length(model_exprs) + 1]] <- quote(theta ~ normal(0, prior_theta_sd))

  model_ast <- as.call(c(list(as.name("{")), model_exprs))

  # --- Mapping of Parameter Names ---
  code_obj <- list(setup = setup_ast, parameters = param_ast, model = model_ast)

  par_names_list <- list()
  par_names_list$theta <- person_names
  if (type == "binary") {
    par_names_list$b <- item_names
  } else {
    par_names_list$b <- list(item_names, paste0("Threshold", 1:(tmp_env$K_cat - 1)))
  }
  if (model %in% c("2PL", "3PL")) par_names_list$a <- item_names
  if (model == "3PL") par_names_list$c <- item_names

  view_vars <- c()
  if (model %in% c("2PL", "3PL")) view_vars <- c(view_vars, "a")
  view_vars <- c(view_vars, "b")
  if (model == "3PL") view_vars <- c(view_vars, "c")

  # --- Automatic Setting of Initial Values ---
  if (is.null(init)) {
    init <- list()
    if (type == "binary") {
      init$b <- rep(0, length(item_names))
    }
    if (model %in% c("2PL", "3PL")) init$a <- rep(1.0, length(item_names))
    if (model == "3PL") init$c <- rep(0.1, length(item_names))
  }

  obj <- rtmb_model(data = as.list(tmp_env), code = code_obj, par_names = par_names_list, init = init, view = view_vars)

  return(obj)
}


#' Wrapper for estimating correlation matrix (multivariate normal distribution)
#'
#' @description
#' Estimates a correlation matrix (along with means and standard deviations) assuming a multivariate normal distribution from observation data.
#' If there are 2 observed variables, it automatically switches to estimate the scalar correlation coefficient (`corr`) directly.
#'
#' @param data Observation data frame or matrix (N x P).
#' @param prior List of hyperparameters for prior distributions. Default is \code{list(lkj_eta = 1.0, mu_sd = 10, sigma_rate = 1.0)}.
#' @param init List of initial values (optional).
#' @param null Target when creating a null model (e.g., \code{"corr"}). Optional.
#' @return An instance of the \code{RTMB_Model} class.
#' @example inst/examples/ex_corr.R
#' @export
rtmb_corr <- function(data, prior = list(lkj_eta = 1.0, mu_sd = 10, sigma_rate = 1.0), init = NULL, null = NULL) {

  Y <- as.matrix(data)
  N <- nrow(Y)
  P <- ncol(Y)

  default_prior <- list(
    lkj_eta = 1.0,
    mu_sd = 10,
    sigma_rate = 1.0
  )
  if (!is.null(prior)) {
    prior <- modifyList(default_prior, prior)
  } else {
    prior <- default_prior
  }

  var_names <- colnames(data)
  if (is.null(var_names)) var_names <- paste0("V", 1:P)

  dat <- list(
    N = N, P = P, Y = Y,
    Y_bar = apply(Y, 2, mean), S_Y = cov(Y) * (N - 1),
    prior_lkj_eta = prior$lkj_eta, prior_sigma_rate = prior$sigma_rate, prior_mu_sd = prior$mu_sd
  )

  if (P == 2) {
    # --- For 2 variables: Estimate correlation coefficient (scalar) directly ---
    code_obj <- rtmb_code(
      parameters = {
        mean = Dim(P)
        sd   = Dim(P, lower = 0)
        corr = Dim(lower = -1, upper = 1)
      },
      model = {
        mean ~ normal(0, prior_mu_sd)
        sd ~ exponential(prior_sigma_rate)
        corr ~ lkj_corr(prior_lkj_eta)

        # Safe matrix initialization to prevent AD type destruction
        CF_corr <- matrix(corr * 0, nrow = 2, ncol = 2)
        CF_corr[1, 1] <- 1
        CF_corr[2, 1] <- corr
        CF_corr[2, 2] <- sqrt(1 - corr^2)

        S_Y ~ sufficient_multi_normal_CF(N, Y_bar, mean, sd, CF_corr)
      }
    )
    par_names_list <- list(mean = var_names, sd = var_names, corr = "rho")

  } else {
    # --- For 3 or more variables: Estimate via Cholesky factor (matrix) ---
    code_obj <- rtmb_code(
      parameters = {
        mean    = Dim(P)
        sd      = Dim(P, lower = 0)
        CF_corr = Dim(c(P, P), type = "CF_corr")
      },
      model = {
        mean ~ normal(0, prior_mu_sd)
        sd ~ exponential(prior_sigma_rate)
        CF_corr ~ lkj_CF_corr(prior_lkj_eta)
        S_Y ~ sufficient_multi_normal_CF(N, Y_bar, mean, sd, CF_corr)
      },
      transform = {
        corr <- CF_corr %*% t(CF_corr)
      }
    )
    par_names_list <- list(mean = var_names, sd = var_names, corr = var_names)
  }

  # Create model instance
  obj <- rtmb_model(
    data = dat,
    code = code_obj,
    par_names = par_names_list,
    init = init,
    view = c("corr")
  )

  # Apply null_model if specified
  if (!is.null(null)) {
    obj <- obj$null_model(target = null)
  }

  return(obj)
}

#' RTMB-based Bayesian two-sample t-test wrapper function
#'
#' @description
#' Performs a Bayesian two-sample t-test using RTMB.
#' It estimates the effect size (delta) with a Cauchy prior, allowing for robust inference
#' and calculation of Bayes factors.
#'
#' @param x Numeric vector of responses for group 1, or a formula (e.g., `y ~ group`).
#' @param y Numeric vector of responses for group 2. Required if `x` is not a formula.
#' @param data Data frame containing the variables in the formula.
#' @param r Numeric; Cauchy prior scale for the effect size (delta). Default is 0.707.
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param use_weak_info Logical; whether to explicitly use weakly informative priors.
#' @param prior List of hyperparameters for the default fixed priors.
#' @param weak_info_prior List of hyperparameters for the weakly informative priors.
#' @param init List of initial values.
#' @param null Character string specifying the target parameter for the null model (e.g., "delta" or "delta ~ cauchy(0, r)").
#' @param y1 Deprecated. Use `x` instead.
#' @param y2 Deprecated. Use `y` instead.
#' @return An \code{RTMB_Model} object.
#' @example inst/examples/ex_ttest.R
#' @export
rtmb_ttest <- function(x, y = NULL, data = NULL, r = 0.707,
                       y_range = NULL,
                       use_weak_info = FALSE,
                       prior = list(mean_sd = 10, sd_rate = 0.1),
                       weak_info_prior = list(sd_ratio = 0.5),
                       init = NULL, null = NULL,
                       y1 = NULL, y2 = NULL) {

  if (!is.null(y1)) x <- y1
  if (!is.null(y2)) y <- y2

  if (inherits(x, "formula")) {
    if (is.null(data)) {
      mf <- model.frame(x, parent.frame())
    } else {
      mf <- model.frame(x, data)
    }
    response <- mf[[1]]
    group <- as.factor(mf[[2]])

    if (length(levels(group)) != 2) {
      stop("The grouping variable must have exactly 2 levels.")
    }

    Y1 <- as.numeric(na.omit(response[group == levels(group)[1]]))
    Y2 <- as.numeric(na.omit(response[group == levels(group)[2]]))
  } else {
    if (is.null(y)) stop("y must be provided if x is not a formula.")
    Y1 <- as.numeric(na.omit(x))
    Y2 <- as.numeric(na.omit(y))
  }

  if (!is.null(y_range)) {
    use_weak_info <- TRUE
  }

  if (use_weak_info && is.null(y_range)) {
    stop("Specifying 'y_range' is required when using weakly informative priors.")
  }

  default_prior <- list(mean_sd = 10, sd_rate = 0.1)
  if (!is.null(prior)) {
    prior <- modifyList(default_prior, prior)
  } else {
    prior <- default_prior
  }

  default_weak_prior <- list(sd_ratio = 0.5)
  if (!is.null(weak_info_prior)) {
    weak_info_prior <- modifyList(default_weak_prior, weak_info_prior)
  } else {
    weak_info_prior <- default_weak_prior
  }


  dat <- list(Y1 = Y1, Y2 = Y2, r = r)
  tmp_env <- list2env(dat, parent = environment())

  # --- Setup AST ---
  if (use_weak_info) {
    setup_exprs <- list()
    setup_exprs[[length(setup_exprs) + 1]] <- quote(half_d_y <- diff(y_range) / 2)
    setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- half_d_y * weak_info_prior$sd_ratio)
    setup_exprs[[length(setup_exprs) + 1]] <- quote(mu_prior_sd <- half_d_y)
    setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mean(y_range))
    setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate <- 1.0 / base_scale)

    setup_ast <- as.call(c(list(as.name("{")), setup_exprs))
    eval(setup_ast, tmp_env)
  } else {
    setup_ast <- NULL
  }

  # --- Parameters AST ---
  param_ast <- quote({
    mean = Dim(1)
    sd = Dim(1, lower = 0)
    delta = Dim(1)
  })

  # --- Transform AST ---
  tran_ast <- quote({
    diff <- delta * sd
    mean0 <- mean - diff / 2
    mean1 <- mean + diff / 2
  })

  # --- Model AST ---
  model_exprs <- list()
  model_exprs[[length(model_exprs) + 1]] <- quote(Y1 ~ normal(mean0, sd))
  model_exprs[[length(model_exprs) + 1]] <- quote(Y2 ~ normal(mean1, sd))
  model_exprs[[length(model_exprs) + 1]] <- quote(delta ~ cauchy(0, r))

  if (use_weak_info) {
    model_exprs[[length(model_exprs) + 1]] <- quote(mean ~ normal(mid_y, mu_prior_sd))
    model_exprs[[length(model_exprs) + 1]] <- quote(sd ~ exponential(sigma_rate))
  } else {
    model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(0, .(prior$mean_sd)))
    model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(.(prior$sd_rate)))
  }

  model_ast <- as.call(c(list(as.name("{")), model_exprs))

  # --- Construction of model object ---
  code_args <- list(parameters = param_ast, transform = tran_ast, model = model_ast)
  if (!is.null(setup_ast)) {
    code_args <- c(list(setup = setup_ast), code_args)
  }

  code_obj <- eval(as.call(c(list(as.name("rtmb_code")), code_args)))

  ordered_data <- env_to_ordered_list(tmp_env, dat, setup_ast)

  view_vars <- c("delta", "mean", "sd")

  obj <- rtmb_model(
    data = ordered_data,
    code = code_obj,
    par_names = list(),
    init = init,
    view = view_vars
  )

  if (!is.null(null)) {
    obj <- obj$null_model(target = null)
  }

  return(obj)
}

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
#' @param prob_formula Optional formula for latent class regression (covariates for probabilities).
#' @param prior_type Prior type: "uniform" (default) or "weakly_informative".
#' @param ... Additional arguments passed to \code{rtmb_model}.
#' @return A \code{RTMB_Model} object.
#' @export
rtmb_mixture <- function(formula, k = 2, data = NULL,
                         covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
                         prior_type = c("uniform", "weakly_informative"), ...) {

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
  prior_type <- match.arg(prior_type)
  K_mix <- k

  # Parse formulas and prepare data
  if (is.matrix(formula) || is.vector(formula)) {
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
    # RHS is formula[[3]] for two-sided formula
    rhs <- formula[[3]]
    is_empty_rhs <- (is.numeric(rhs) && rhs == 1) || (is.symbol(rhs) && rhs == "1")
    
    if (!is_empty_rhs) {
      X_prob <- model.matrix(formula, data = data)
    } else {
      X_prob <- NULL
    }
  }

  Y_mat <- as.matrix(Y_mat)
  N_obs <- nrow(Y_mat)
  P_dim <- ncol(Y_mat)
  
  if (is.null(colnames(Y_mat))) {
    colnames(Y_mat) <- paste0("Y", 1:P_dim)
  }

  multivariate <- P_dim > 1
  has_cov_prob <- !is.null(X_prob)
  is_sigma_equal <- covariance %in% c("diagonal_equal", "full_equal", "full_equal_corr")

  # --- 1. Dynamic AST Construction: setup ---
  setup_exprs <- list(
    as.name("{"),
    bquote(N <- .(N_obs)),
    bquote(K <- .(K_mix)),
    bquote(P <- .(P_dim))
  )
  if (has_cov_prob) {
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(X_means <- matrix(.(colMeans(X_prob)), 1, .(ncol(X_prob))))
  }
  setup_ast <- as.call(setup_exprs)

  # --- 2. Dynamic AST Construction: parameters ---
  param_exprs <- list(as.name("{"))
  if (has_cov_prob) {
    param_exprs[[length(param_exprs) + 1]] <- bquote(b <- Dim(c(.(ncol(X_prob)), K - 1)))
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

  if (prior_type != "uniform") {
    model_exprs[[length(model_exprs) + 1]] <- quote(mu ~ normal(0, mu_sd))
    model_exprs[[length(model_exprs) + 1]] <- quote(sigma ~ exponential(sigma_rate))
    if (multivariate && covariance %in% c("full", "full_equal", "full_equal_corr")) {
      if (covariance == "full") {
        model_exprs[[length(model_exprs) + 1]] <- quote(for (k in 1:K) {
          matrix(L_corr[k, , ], P, P) ~ lkj_CF_corr(2)
        })
      } else {
        model_exprs[[length(model_exprs) + 1]] <- quote(L_corr ~ lkj_CF_corr(2))
      }
    }
    if (has_cov_prob) {
      model_exprs[[length(model_exprs) + 1]] <- quote(b ~ normal(0, beta_sd))
    } else {
      model_exprs[[length(model_exprs) + 1]] <- quote(theta ~ dirichlet(rep(1, K)))
    }
  }
  model_ast <- as.call(model_exprs)

  # --- 4. Dynamic AST Construction: generate ---
  generate_exprs <- list(as.name("{"))
  if (has_cov_prob) {
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
  generate_ast <- as.call(generate_exprs)

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

  data_list <- list(Y = Y_mat)
  if (has_cov_prob) data_list$X_prob <- X_prob

  init_list <- list()
  if (nrow(Y_mat) >= K_mix) {
    # Use kmeans for better initial values
    km <- try(kmeans(Y_mat, centers = K_mix, nstart = 5), silent = TRUE)
    if (!inherits(km, "try-error")) {
      if (multivariate) {
        init_list$mu <- km$centers
      } else {
        init_list$mu <- as.vector(km$centers)
      }
      # Initial theta / probabilities
      props <- as.vector(table(factor(km$cluster, levels = 1:K_mix))) / nrow(Y_mat)
      init_list$theta <- props
      
      # Initial b if covariates present
      if (has_cov_prob) {
        b_init <- matrix(0, ncol(X_prob), K_mix - 1)
        for (k in 2:K_mix) {
          z_k <- as.numeric(km$cluster == k)
          fit_k <- try(suppressWarnings(glm(z_k ~ X_prob - 1 + offset(rep(qlogis(max(0.01, props[k])), length(z_k))), family = binomial)), silent = TRUE)
          if (!inherits(fit_k, "try-error")) {
            b_init[, k-1] <- coef(fit_k)
          }
        }
        init_list$b <- b_init
      }
    }
  }

  if (is.null(init_list$mu)) {
    # Fallback to even spread
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

  mdl <- BayesRTMB::rtmb_model(data_list, mdl_code, par_names = v_names, init = init_list, view = view_order)
  return(mdl)
}

#' Fit a Latent Rank Theory (LRT) Model
#'
#' @description
#' Fits a Latent Rank Theory model, which is a mixture model with ordered ranks
#' and Gaussian Process smoothing on the mean profiles.
#'
#' @param formula A formula specifying the response variable(s).
#' @param k Number of ranks (mixture components).
#' @param data A data frame containing the variables.
#' @param magnitude Signal standard deviation for the GP prior. If NULL, it is estimated.
#' @param smoothing Length-scale for the GP prior. If NULL, it is estimated.
#' @param noise Measurement noise for the GP prior (default is 0.01).
#' @param prob_formula Optional formula for latent class regression.
#' @param prior_type Prior type: "weakly_informative" or "uniform".
#' @param ... Additional arguments passed to \code{rtmb_model}.
#' @return A \code{RTMB_Model} object.
#' @export
rtmb_lrt <- function(formula, k = 3, data = NULL,
                     covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
                     magnitude = NULL, smoothing = NULL, noise = 0.01,
                     prob_smoothing = FALSE,
                     link = c("ordered", "sequential"),
                     prior_type = c("weakly_informative", "uniform"), ...) {

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
  prior_type <- match.arg(prior_type)
  link <- match.arg(link)
  K_mix <- k

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

  Y_mat <- as.matrix(Y_mat)
  N_obs <- nrow(Y_mat)
  P_dim <- ncol(Y_mat)
  if (is.null(colnames(Y_mat))) colnames(Y_mat) <- paste0("Y", 1:P_dim)
  has_cov_prob <- !is.null(X_prob)

  multivariate <- P_dim > 1
  is_diag <- covariance %in% c("diagonal", "diagonal_equal")
  is_sigma_equal <- covariance %in% c("diagonal_equal", "full_equal", "full_equal_corr")

  # --- 1. Setup ---
  setup_exprs <- list(
    as.name("{"),
    bquote(N <- .(N_obs)),
    bquote(K <- .(K_mix)),
    bquote(P <- .(P_dim)),
    bquote(rank_coords <- 1:.(K_mix))
  )
  if (has_cov_prob) {
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(X_means <- matrix(.(colMeans(X_prob)), 1, .(ncol(X_prob))))
  }
  setup_ast <- as.call(setup_exprs)

  # --- 2. Parameters ---
  param_exprs <- list(as.name("{"))
  param_exprs[[length(param_exprs) + 1]] <- bquote(mu_p <- Dim(c(P, K), type = "ordered"))
  
  # Sigma parameterization based on covariance
  if (is_sigma_equal) {
    param_exprs[[length(param_exprs) + 1]] <- bquote(sigma <- Dim(P, lower = 0))
  } else {
    param_exprs[[length(param_exprs) + 1]] <- bquote(sigma <- Dim(c(K, P), lower = 0))
  }

  if (multivariate && !is_diag) {
    if (covariance == "full") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr <- Dim(c(K, P, P), type = "CF_corr"))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr <- Dim(c(P, P), type = "CF_corr"))
    }
  }
  
  if (has_cov_prob) {
    if (link == "ordered") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(b <- Dim(.(ncol(X_prob))))
      param_exprs[[length(param_exprs) + 1]] <- bquote(cutpoints <- Dim(.(K_mix - 1), type = "ordered"))
    } else {
      # Sequential: Separate b for each step
      param_exprs[[length(param_exprs) + 1]] <- bquote(b <- Dim(c(.(ncol(X_prob)), .(K_mix - 1))))
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
  model_exprs[[length(model_exprs) + 1]] <- bquote({
    for (p in 1:P) {
      lp <- lp + gaussian_process_lpdf(mu[, p], x = rank_coords, magnitude = .(gp_mag), smoothing = .(gp_sm), noise = .(noise))
    }
  })

  if (prob_smoothing && !has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- bquote({
      lp <- lp + gaussian_process_lpdf(logit_theta, x = rank_coords, magnitude = magnitude_theta, smoothing = .(gp_sm), noise = .(noise))
    })
  }

  if (prior_type != "uniform") {
    model_exprs[[length(model_exprs) + 1]] <- quote(sigma ~ exponential(1))
    if (is.null(magnitude)) model_exprs[[length(model_exprs) + 1]] <- quote(magnitude ~ exponential(1))
    if (is.null(smoothing)) model_exprs[[length(model_exprs) + 1]] <- quote(smoothing ~ exponential(1))
    if (prob_smoothing && !has_cov_prob) {
      model_exprs[[length(model_exprs) + 1]] <- quote(magnitude_theta ~ exponential(1))
    }
    if (has_cov_prob) {
      model_exprs[[length(model_exprs) + 1]] <- quote(b ~ normal(0, 5))
      if (link == "ordered") model_exprs[[length(model_exprs) + 1]] <- quote(cutpoints ~ normal(0, 5))
      else model_exprs[[length(model_exprs) + 1]] <- quote(alpha ~ normal(0, 5))
    } else if (!prob_smoothing) {
      model_exprs[[length(model_exprs) + 1]] <- quote(theta ~ dirichlet(rep(1, K)))
    }
  }
  model_ast <- as.call(model_exprs)

  # --- 4. Generate ---
  generate_exprs <- list(as.name("{"))
  generate_exprs[[length(generate_exprs) + 1]] <- quote(mu <- t(mu_p))
  if (has_cov_prob) {
    if (link == "ordered") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote({
        eta_mean <- X_means %*% b
        F_prev <- 0
        prob_mean <- numeric(K)
        for (k in 1:(K-1)) {
          F_k <- inv_logit(cutpoints[k] - eta_mean)
          prob_mean[k] <- F_k - F_prev
          F_prev <- F_k
        }
        prob_mean[K] <- 1 - F_prev
      })
    } else {
      generate_exprs[[length(generate_exprs) + 1]] <- quote({
        eta_mean_mat <- X_means %*% b
        q_cum <- 1
        prob_mean <- numeric(K)
        for (k in 1:(K-1)) {
          q_k <- inv_logit(alpha[k] + eta_mean_mat[k])
          prob_mean[k] <- q_cum * (1 - q_k)
          q_cum <- q_cum * q_k
        }
        prob_mean[K] <- q_cum
      })
    }
  } else {
    if (prob_smoothing) generate_exprs[[length(generate_exprs) + 1]] <- quote(theta <- softmax(logit_theta))
    generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- theta)
  }
  
  if (!is_diag) {
    generate_exprs[[length(generate_exprs) + 1]] <- bquote({
      corr <- if (.(covariance == "full")) array(mu[1] * 0, dim = c(K, P, P)) else matrix(mu[1] * 0, P, P)
      if (.(covariance == "full")) {
        for (k in 1:K) {
          L_mat <- matrix(L_corr[k, , ], P, P)
          corr[k, , ] <- L_mat %*% t(L_mat)
        }
      } else {
        L_mat <- matrix(L_corr, P, P)
        corr <- L_mat %*% t(L_mat)
      }
    })
  }
  generate_ast <- as.call(generate_exprs)

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

  data_list <- list(Y = Y_mat)
  if (has_cov_prob) data_list$X_prob <- X_prob

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
          
          if (link == "ordered") {
            init_list$b <- as.vector(b_est)
          } else {
            init_list$b <- matrix(b_est, ncol(X_prob), K_mix - 1)
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
  
  mdl <- BayesRTMB::rtmb_model(data_list, mdl_code, par_names = v_names, init = init_list, view = view_order)
  return(mdl)
}
