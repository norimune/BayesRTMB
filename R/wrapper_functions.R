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

#' RTMB-based GLMM wrapper function
#'
#' @param formula lme4-style formula (e.g., Y ~ X + (1 | GID))
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param laplace Logical; whether to marginalize random effects using Laplace approximation
#' @param penalty Type of regularization for fixed effects: "none", "rhs" (Regularized Horseshoe), or "ssp" (Spike and Slab Prior). Default is "none".
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param use_weak_info Logical; whether to explicitly use weakly informative priors (requires y_range for continuous models).
#' @param prior List of hyperparameters for the default fixed priors.
#' @param weak_info_prior List of hyperparameters for the weakly informative priors and regularization.
#' @param init List of initial values (generated automatically based on glm if omitted)
#' @param null Character string specifying the target parameter for the null model.
#' @import reformulas
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glmer <- function(formula, data, family = "gaussian", laplace = FALSE,
                       penalty = c("none", "rhs", "ssp"),
                       y_range = NULL,
                       use_weak_info = FALSE,
                       prior = list(Intercept_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, nu_rate = 0.1, cutpoint_sd = 2.5, shape_rate = 1.0, phi_rate = 1.0, lkj_eta = 1.0),
                       weak_info_prior = list(max_beta = 1, sd_ratio = 0.5, expected_vars = 3, slab_scale = 2.0, slab_df = 4.0, ssp_ratio = 0.25),
                       init = NULL, null = NULL) {

  regularization <- match.arg(penalty)
  if (!requireNamespace("reformulas", quietly = TRUE)) stop("The 'reformulas' package is required to parse the formula.")
  has_random <- !is.null(reformulas::findbars(formula))

  default_prior <- eval(formals(rtmb_glmer)$prior)
  prior <- modifyList(default_prior, prior)
  default_weak_prior <- eval(formals(rtmb_glmer)$weak_info_prior)
  weak_info_prior <- modifyList(default_weak_prior, weak_info_prior)

  if (!is.null(y_range) || regularization %in% c("rhs", "ssp")) {
    use_weak_info <- TRUE
  }

  families_requiring_yrange <- c("gaussian", "lognormal", "student_t")
  is_continuous <- family %in% families_requiring_yrange

  if (use_weak_info && is_continuous && is.null(y_range)) {
    stop("Specifying 'y_range' is required when using weakly informative priors (or regularization) with family = '", family, "'.")
  }

  if (has_random) {
    parsed <- lme4::lFormula(formula, data = data, control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE  = "ignore"))
    Y <- model.response(parsed$fr)
    X <- parsed$X
    reTrms <- parsed$reTrms

    if (length(reTrms$Ztlist) > 1) stop("Currently, only a single grouping variable is supported.")

    Zt <- as.matrix(reTrms$Ztlist[[1]])
    group_idx <- as.integer(reTrms$flist[[1]])
    offset <- model.offset(parsed$fr)
    ranef_names <- parsed$reTrms$cnms[[1]]
    ranef_names[ranef_names == "(Intercept)"] <- "Int"
  } else {
    mf <- model.frame(formula, data)
    Y <- model.response(mf)
    X <- model.matrix(formula, mf)
    offset <- model.offset(mf)
    Zt <- NULL
    group_idx <- NULL
    ranef_names <- NULL
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
      rhs_mod <- rtmb_glmer(formula, data, family, laplace, penalty = "rhs", y_range = y_range, use_weak_info = use_weak_info, prior = prior, weak_info_prior = weak_info_prior, init = NULL)
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
      fixed_formula <- if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars(formula) else suppressWarnings(reformulas::nobars(formula))
      glm_fam <- switch(family, "gaussian" = gaussian(), "student_t" = gaussian(), "lognormal" = gaussian(link="log"), "bernoulli" = binomial(), "binomial" = binomial(), "poisson" = poisson(), "neg_binomial" = poisson(), "gamma" = Gamma(link="log"), gaussian())
      init_fit <- suppressWarnings(glm(fixed_formula, data = data, family = glm_fam))
      init_coef <- unname(coef(init_fit))
      init_coef[is.na(init_coef)] <- 0

      init <- list()
      if (has_intercept) init$Intercept_c <- init_coef[1]

      if (regularization == "none" && K_tmp > 0) {
        if (has_intercept) init$b <- init_coef[-1]
        else init$b <- init_coef
      }
    }, error = function(e) init <- NULL)
  }

  dat <- list(Y = Y, trials = trials, X = X)
  if (has_random) {
    dat$Zt <- Zt
    dat$group_idx <- group_idx
  }
  if (!is.null(offset)) dat$offset <- offset
  if (family == "ordered") dat$num_categories <- length(unique(Y))

  # --- Setup AST ---
  setup_exprs <- list()
  setup_exprs[[1]] <- quote(N <- length(Y))
  setup_exprs[[2]] <- quote(K <- ncol(X))
  if (has_random) {
    setup_exprs[[3]] <- quote(num_groups <- max(group_idx))
    setup_exprs[[4]] <- quote(num_ranef <- nrow(Zt) / num_groups)
    setup_exprs[[5]] <- quote(Z_mat <- matrix(0, nrow = N, ncol = num_ranef))
    setup_exprs[[6]] <- quote(
      for (i in 1:N) {
        g <- group_idx[i]
        Z_mat[i, ] <- Zt[((g - 1) * num_ranef + 1):(g * num_ranef), i]
      }
    )
  }

  if (use_weak_info) {
    if (family %in% c("bernoulli", "binomial", "ordered")) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- pi / sqrt(3))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- base_scale * weak_info_prior$max_beta)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- 0)
    } else if (family %in% c("poisson", "neg_binomial", "gamma")) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- 1.0)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- base_scale)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- 0)
    } else {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(half_d_y <- diff(y_range) / 2)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- half_d_y * weak_info_prior$sd_ratio)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- half_d_y)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mean(y_range))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate <- 1.0 / base_scale)
    }

    setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_rate <- 1.0 / base_scale)

    if (K_tmp > 0) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_mean <- apply(X, 2, mean))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_sd <- apply(X, 2, sd))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(beta_prior_sd <- weak_info_prior$max_beta * base_scale / X_sd)
      if (has_intercept) {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(X_c <- X - rep(1, N) %*% t(X_mean))
      }

      if (regularization == "rhs") {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(p0 <- min(weak_info_prior$expected_vars, K - 1))
        setup_exprs[[length(setup_exprs) + 1]] <- quote(if (p0 < 1) p0 <- 1)
        setup_exprs[[length(setup_exprs) + 1]] <- quote(tau0 <- (p0 / (K - p0)) * (base_scale / sqrt(N)))
        setup_exprs[[length(setup_exprs) + 1]] <- quote(half_slab_df <- weak_info_prior$slab_df / 2.0)
        setup_exprs[[length(setup_exprs) + 1]] <- quote(half_slab_scale2 <- (weak_info_prior$slab_df * weak_info_prior$slab_scale^2) / 2.0)
      } else if (regularization == "ssp") {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_scale <- base_scale / X_sd)
      }
    }
  } else {
    # Minimal calculations for fixed prior distributions
    if (K_tmp > 0) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_mean <- apply(X, 2, mean))
      if (has_intercept) {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(X_c <- X - rep(1, N) %*% t(X_mean))
      }
    }
  }
  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(dat)
  eval(setup_ast, tmp_env)
  N <- tmp_env$N; K <- tmp_env$K
  if (has_random) {
    num_groups <- tmp_env$num_groups; num_ranef <- tmp_env$num_ranef
  }

  # --- Parameters AST ---
  param_exprs <- list()
  if (has_intercept) param_exprs[[length(param_exprs) + 1]] <- quote(Intercept_c <- Dim(1))
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
    if (num_ranef == 1) {
      param_exprs[[length(param_exprs) + 1]] <- quote(sd <- Dim(num_ranef, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(r_re <- Dim(num_groups, random = TRUE))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- quote(sd <- Dim(num_ranef, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(CF_corr <- Dim(c(num_ranef, num_ranef), type = "CF_corr"))
      param_exprs[[length(param_exprs) + 1]] <- quote(r_re <- Dim(c(num_groups, num_ranef), random = TRUE))
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
    if (K > 0) tran_exprs[[length(tran_exprs) + 1]] <- quote(Intercept <- Intercept_c - sum(X_mean * b))
    else tran_exprs[[length(tran_exprs) + 1]] <- quote(Intercept <- Intercept_c)
  }
  if (has_random && num_ranef > 1) tran_exprs[[length(tran_exprs) + 1]] <- quote(corr <- CF_corr %*% t(CF_corr))
  tran_ast <- if (length(tran_exprs) > 0) as.call(c(list(as.name("{")), tran_exprs)) else NULL

  # --- Model AST ---
  transform_exprs <- list()
  if (has_intercept) {
    if (K > 0) transform_exprs[[length(transform_exprs) + 1]] <- quote(eta <- as.vector(Intercept_c + X_c %*% b))
    else transform_exprs[[length(transform_exprs) + 1]] <- quote(eta <- rep(Intercept_c, N))
  } else {
    if (K > 0) transform_exprs[[length(transform_exprs) + 1]] <- quote(eta <- as.vector(X %*% b))
    else transform_exprs[[length(transform_exprs) + 1]] <- quote(eta <- rep(0, N))
  }
  if (has_random) {
    if (num_ranef > 1) {
      transform_exprs[[length(transform_exprs) + 1]] <- quote(for (i in 1:N) eta[i] <- eta[i] + sum(Z_mat[i, ] * r_re[group_idx[i], ] * sd))
    } else {
      transform_exprs[[length(transform_exprs) + 1]] <- quote(eta <- eta + Z_mat[,1] * r_re[group_idx] * sd)
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
    if (num_ranef > 1) {
      ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(CF_corr ~ lkj_CF_corr(.(prior$lkj_eta)))
      ll_random_exprs[[length(ll_random_exprs) + 1]] <- quote(for (j in 1:num_groups) r_re[j, ] ~ multi_normal_CF(rep(0, num_ranef), rep(1, num_ranef), CF_corr))
    } else {
      ll_random_exprs[[length(ll_random_exprs) + 1]] <- quote(r_re ~ normal(0, 1))
    }
  }

  prior_exprs <- list()
  if (family %in% c("gaussian", "lognormal", "student_t")) {
    if (use_weak_info) {
      prior_exprs[[length(prior_exprs) + 1]] <- quote(sigma ~ exponential(sigma_rate))
    } else {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    }
  }

  if (family == "student_t") prior_exprs[[length(prior_exprs) + 1]] <- bquote(nu ~ exponential(.(prior$nu_rate)))
  if (family == "gamma") prior_exprs[[length(prior_exprs) + 1]] <- bquote(shape ~ exponential(.(prior$shape_rate)))
  if (family == "neg_binomial") prior_exprs[[length(prior_exprs) + 1]] <- bquote(phi ~ exponential(.(prior$phi_rate)))
  if (family == "ordered") prior_exprs[[length(prior_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))

  if (has_intercept) {
    if (use_weak_info) {
      prior_exprs[[length(prior_exprs) + 1]] <- quote(Intercept_c ~ normal(mid_y, alpha_prior_sd))
    } else {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(Intercept_c ~ normal(0, .(prior$Intercept_sd)))
    }
  }

  if (K > 0) {
    if (regularization == "none") {
      if (use_weak_info) {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(b ~ normal(0, beta_prior_sd))
      } else {
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
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(mu_r <- log(.(weak_info_prior$ssp_ratio) / (1 - .(weak_info_prior$ssp_ratio))))
      prior_exprs[[length(prior_exprs) + 1]] <- quote(r ~ logit_normal(mu_r, 3))
      prior_exprs[[length(prior_exprs) + 1]] <- quote(tau ~ exponential(1 / tau_scale))
    }
  }

  if (has_random) {
    if (use_weak_info) {
      prior_exprs[[length(prior_exprs) + 1]] <- quote(sd ~ exponential(tau_rate))
    } else {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(sd ~ exponential(.(prior$sd_rate)))
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
    par_names_list$sd <- ranef_names
    if (num_ranef > 1) par_names_list$corr <- ranef_names
  }

  view_vars <- c()
  if (has_intercept) view_vars <- c("Intercept")
  if (K > 0) view_vars <- c(view_vars, "b")
  view_vars <- c(view_vars, "sigma")
  if (has_random) {
    view_vars <- c(view_vars, "sd", "corr")
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

#' RTMB-based GLM wrapper function (no random effects)
#'
#' @param formula Formula
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param penalty Type of regularization for fixed effects: "none", "rhs" (Regularized Horseshoe), or "ssp" (Spike and Slab Prior). Default is "none".
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param use_weak_info Logical; whether to explicitly use weakly informative priors (requires y_range for continuous models).
#' @param prior List of hyperparameters for the default fixed priors.
#' @param weak_info_prior List of hyperparameters for the weakly informative priors and regularization.
#' @param init List of initial values
#' @param null Character string specifying the target parameter for the null model.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glm <- function(formula, data, family = "gaussian",
                     penalty = c("none", "rhs", "ssp"),
                     y_range = NULL,
                     use_weak_info = FALSE,
                     prior = list(Intercept_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, nu_rate = 0.1, cutpoint_sd = 2.5, shape_rate = 1.0, phi_rate = 1.0, lkj_eta = 1.0),
                     weak_info_prior = list(max_beta = 1, sd_ratio = 0.5, expected_vars = 3, slab_scale = 2.0, slab_df = 4.0, ssp_ratio = 0.25),
                     init = NULL, null = NULL) {
  regularization <- match.arg(penalty)
  rtmb_glmer(
    formula = formula,
    data = data,
    family = family,
    laplace = FALSE,
    penalty = regularization,
    y_range = y_range,
    use_weak_info = use_weak_info,
    prior = prior,
    weak_info_prior = weak_info_prior,
    init = init,
    null = null
  )
}

#' RTMB-based Linear Regression wrapper function
#'
#' @param formula Formula (e.g., Y ~ X1 + X2)
#' @param data Data frame
#' @param penalty Type of regularization for fixed effects: "none", "rhs" (Regularized Horseshoe), or "ssp" (Spike and Slab Prior). Default is "none".
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param use_weak_info Logical; whether to explicitly use weakly informative priors.
#' @param prior List of hyperparameters for the default fixed priors.
#' @param weak_info_prior List of hyperparameters for the weakly informative priors and regularization.
#' @param init List of initial values.
#' @param null Character string specifying the target parameter for the null model.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_lm <- function(formula, data,
                    penalty = c("none", "rhs", "ssp"),
                    y_range = NULL,
                    use_weak_info = FALSE,
                    prior = list(Intercept_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, nu_rate = 0.1, cutpoint_sd = 2.5, shape_rate = 1.0, phi_rate = 1.0),
                    weak_info_prior = list(max_beta = 1, sd_ratio = 0.5, expected_vars = 3, slab_scale = 2.0, slab_df = 4.0, ssp_ratio = 0.25),
                    init = NULL, null = NULL) {

  regularization <- match.arg(penalty)
  rtmb_glm(
    formula = formula,
    data = data,
    family = "gaussian",
    y_range = y_range,
    use_weak_info = use_weak_info,
    penalty = regularization,
    prior = prior,
    weak_info_prior = weak_info_prior,
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
    Y_bar <- apply(Y, 2, mean)
    S_Y <- cov(Y) * (N - 1)
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
        stop("Rotation function not found: ", rotate)
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
