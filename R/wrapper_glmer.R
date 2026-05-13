#' RTMB-based GLMM wrapper function
#'
#' @param formula lme4-style formula (e.g., Y ~ X + (1 | GID))
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param laplace Logical; whether to marginalize random effects using Laplace approximation
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Use prior_weak(), prior_rhs(), or prior_ssp(). Default is `prior_flat()`.
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Required when using weakly informative or regularized priors with continuous models.
#' @param init List of initial values (generated automatically based on glm if omitted)

#' @param gmc Character vector of variable names for Grand Mean Centering (GMC). If "all", all numeric variables are centered.
#' @param cwc List for Centering Within Cluster (CWC). Should contain \code{cluster} (group variable) and \code{pars} (variable names to center).
#' @param view Character vector of parameter names to prioritize in summary.
#' @param factors Character vector of variable names to be treated as factors.
#' @param contrasts Character string specifying the contrast type ("treatment" or "sum").
#' @param sigma_by Character vector specifying variables to group residual variance by (heteroscedasticity).
#' @param resid_corr Residual correlation structure: "ar1" (Autoregressive), "cs" (Compound Symmetry), "toep" (Toeplitz), or "un" (Unstructured).
#' @param resid_time Variable name for time points in residual correlation.
#' @param resid_group Variable name for grouping in residual correlation.
#' @param within Optional list for wide-to-long conversion.
#' @param generate Optional expression for generated quantities.
#' @param .force_sum Logical; internal use only.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glmer <- function(formula, data, family = "gaussian", laplace = FALSE,
                       prior = prior_flat(),
                       y_range = NULL,
                       init = NULL,
                       fixed = NULL,
                       gmc = NULL,
                       cwc = NULL,
                       view = NULL,
                       within = NULL,
                       factors = NULL,
                       contrasts = "treatment",
                       sigma_by = NULL,
                       resid_corr = NULL,
                       resid_time = NULL,
                       resid_group = NULL,
                       generate = NULL,
                       .force_sum = FALSE) {

  # --- Handle Wide to Long conversion if needed ---
  if (family == "gaussian" && is.call(formula[[2]]) && identical(formula[[2]][[1]], as.name("cbind"))) {
    processed <- .handle_wide_to_long(formula, data, within, factors)
    formula <- processed$formula
    data <- processed$data
    factors <- processed$factors
  }

  # --- 0. Contrast Management (Automatic sum-to-zero) ---


   if (!is.null(resid_corr)) {
     valid_resid_corr <- c("ar1", "cs", "un", "toep")
     if (!(resid_corr %in% valid_resid_corr)) {
       stop(sprintf("Invalid 'resid_corr' value: '%s'. Valid options are: %s",
            resid_corr, paste(valid_resid_corr, collapse = ", ")))
     }
   }

  # --- 0. Data Preparation (Factors) ---
  if (!is.null(factors)) {
    data <- as.data.frame(data)
    for (f in factors) {
      if (f %in% names(data)) {
        data[[f]] <- as.factor(data[[f]])
      }
    }
  }

  # --- 0. Data Centering (GMC / CWC) ---
  if (!is.null(gmc) || !is.null(cwc)) {
    data_centered <- as.data.frame(data)

    # 1. Grand Mean Centering
    if (!is.null(gmc)) {
      target_gmc <- if (identical(gmc, "all")) {
        names(data_centered)[sapply(data_centered, is.numeric)]
      } else {
        gmc
      }
      for (v in target_gmc) {
        if (v %in% names(data_centered)) {
          data_centered[[v]] <- data_centered[[v]] - mean(data_centered[[v]], na.rm = TRUE)
        } else {
          warning(sprintf("Variable '%s' for GMC not found in data.", v))
        }
      }
    }

    # 2. Centering Within Cluster
    if (!is.null(cwc)) {
      # Robustly extract cluster and pars
      if ("pars" %in% names(cwc)) {
        cluster_var <- cwc$cluster
        target_pars <- cwc$pars
      } else if ("cluster" %in% names(cwc)) {
        cluster_var <- cwc$cluster
        # Any elements not named "cluster" are considered target_pars
        target_pars <- unlist(cwc[names(cwc) != "cluster" | names(cwc) == ""])
      } else if (length(cwc) >= 2) {
        # Assume first is cluster, rest are pars
        cluster_var <- cwc[[1]]
        target_pars <- unlist(cwc[-1])
      } else {
        stop("Invalid 'cwc' format. Use list(cluster = ID, pars = c('var1', 'var2')) or list(ID, 'var1').")
      }

      # Resolve group ID (character or symbol)
      group_id <- if (is.character(cluster_var)) {
        data_centered[[cluster_var]]
      } else {
        # cluster_var might be a symbol (e.g., group)
        # Try to evaluate it in the context of data_centered
        tryCatch(eval(cluster_var, data_centered, parent.frame()),
                 error = function(e) NULL)
      }

      if (is.null(group_id)) {
        # If still NULL, try to deparse and check if it's a column name
        group_var_name <- deparse(cluster_var)
        if (group_var_name %in% names(data_centered)) {
          group_id <- data_centered[[group_var_name]]
        } else {
           stop("Cluster variable for CWC not found.")
        }
      }

      for (v in target_pars) {
        if (v %in% names(data_centered)) {
          group_means <- tapply(data_centered[[v]], group_id, mean, na.rm = TRUE)
          data_centered[[v]] <- data_centered[[v]] - group_means[as.character(group_id)]
        } else {
          warning(sprintf("Variable '%s' for CWC not found in data.", v))
        }
      }
    }
    data <- data_centered
  }

  # --- 0.1 Handle sigma_by (Heteroscedasticity) ---
  sigma_idx <- NULL
  num_sigma_groups <- 1
  sigma_group_names <- "sigma"
  if (!is.null(sigma_by)) {
    if (!family %in% c("gaussian", "lognormal", "student_t")) {
      stop("'sigma_by' is only supported for continuous families (gaussian, lognormal, student_t).")
    } else {
      # Identify categorical factors in sigma_by
      if (identical(sigma_by, "all")) {
        # Extract all variables in fixed part of formula and check if they are factors
        vars_in_form <- all.vars(nobars(formula))
        # Exclude response (the first variable in all.vars is the response)
        vars_in_form <- vars_in_form[-1]
        is_cat <- sapply(data[vars_in_form], function(x) is.factor(x) || is.character(x))
        target_vars <- vars_in_form[is_cat]
      } else {
        target_vars <- sigma_by
      }

      if (length(target_vars) > 0) {
        sigma_factor <- if (length(target_vars) > 1) {
          interaction(data[target_vars], drop = TRUE, sep = ":")
        } else {
          as.factor(data[[target_vars]])
        }
        sigma_idx <- as.integer(sigma_factor)
        num_sigma_groups <- nlevels(sigma_factor)
        sigma_group_names <- levels(sigma_factor)
      }
    }
  }

  # --- 0.2 Default Initial Values ---
  if (is.null(init)) {
    init <- list()
  }
  if (is.list(init)) {
    # Set safe default for residual correlation (e.g. CS, AR1) if not specified by user.
    # This avoids non-positive definite matrices during the initial check.
    if (!is.null(resid_corr) && !("rho_resid" %in% names(init)) && resid_corr %in% c("ar1", "cs", "toep")) {
      init$rho_resid <- 0.1
    }
  }

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

  prior_type <- prior$type

  is_flat <- identical(prior_type, "flat")
  is_normal <- identical(prior_type, "normal")
  is_weak <- identical(prior_type, "weak")
  is_jzs <- identical(prior_type, "jzs")

  regularization <- if (prior_type %in% c("rhs", "ssp")) prior_type else "none"

  use_centering <- prior_type %in% c("normal", "weak", "rhs", "ssp", "jzs")
  use_weak_info <- prior_type %in% c("weak", "rhs", "ssp", "jzs")
  use_normal_prior <- identical(prior_type, "normal")

  if (use_normal_prior) {
    if (is.null(prior$mu_sd) && !is.null(prior$Intercept_sd)) {
      prior$mu_sd <- prior$Intercept_sd
    }
  }

  has_random <- !is.null(findbars(formula))

  # If classic mode and random effects are present (or forced), we internally
  # force 'sum' contrasts for stable Type III ANOVA/DFs.
  # Otherwise, we use the user-requested contrast (defaulting to 'treatment').
  actual_contrasts <- if (.force_sum) "sum" else contrasts

  if (!is.null(actual_contrasts)) {
    if (actual_contrasts == "sum") {
      old_opts <- options(contrasts = c("contr.sum", "contr.poly"))
      on.exit(options(old_opts), add = TRUE)
    } else if (actual_contrasts == "treatment") {
      old_opts <- options(contrasts = c("contr.treatment", "contr.poly"))
      on.exit(options(old_opts), add = TRUE)
    }
  }

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

  if (use_weak_info && is_continuous && is.null(y_range) && prior_type != "jzs") {
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

    # Add residual correlation variables to formula for model.frame if needed
    mf_form <- all_vars_form
    if (!is.null(resid_group)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_group)))
    if (!is.null(resid_time)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_time)))

    mf <- model.frame(mf_form, data = data)
    Y <- model.response(mf)
    X <- model.matrix(fixed_form, mf)
    X_assign <- attr(X, "assign")
    X_terms <- attr(terms(fixed_form), "term.labels")
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
    # Add residual correlation variables to formula for model.frame if needed
    mf_form <- formula
    if (!is.null(resid_group)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_group)))
    if (!is.null(resid_time)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_time)))

    mf <- model.frame(mf_form, data)
    Y <- model.response(mf)
    X <- model.matrix(formula, mf)
    X_assign <- attr(X, "assign")
    X_terms <- attr(terms(formula), "term.labels")
    offset <- model.offset(mf)
    N <- nrow(mf)
    num_bars <- 0
    dat_ranef <- list()
    num_ranef_list <- list()
    num_groups_list <- list()
    ranef_names_list <- list()
    group_labels_list <- list()
  }
  
  jzs_V_val <- NULL

  has_intercept <- "(Intercept)" %in% colnames(X)
  if (family == "ordered") has_intercept <- FALSE
  fixed_colnames <- colnames(X)
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
        if (use_centering) {
          init$Intercept_c <- init_coef[1]
        } else {
          init$Intercept <- init_coef[1]
        }
      }

      if (regularization == "none" && K_tmp > 0) {
        if (has_intercept) init$b <- init_coef[-1]
        else init$b <- init_coef
      }

      # Initial values for residual correlation
      if (!is.null(resid_corr)) {
        if (resid_corr %in% c("ar1", "cs")) init$rho_resid <- 0.1
        else if (resid_corr == "toep") init$rho_resid <- rep(0.1, dat$max_T_resid - 1)
        else if (resid_corr == "un") init$L_resid <- diag(dat$max_T_resid)
      }

      if (family %in% c("gaussian", "lognormal", "student_t")) {
        init$sigma <- if (num_sigma_groups > 1) rep(stats::sd(Y) * 0.5, num_sigma_groups) else stats::sd(Y) * 0.5
      }

    }, error = function(e) init <- NULL)
  }

  dat <- list(Y = Y, trials = trials, X = X,
              sigma_idx = sigma_idx, num_sigma_groups = num_sigma_groups)

  if (!is.null(resid_corr)) {
    if (family != "gaussian") stop("Residual correlation structures are currently only supported for Gaussian models.")
    resid_corr <- tolower(resid_corr)
    if (!resid_corr %in% c("ar1", "cs", "un", "toep")) stop("Currently 'ar1', 'cs', 'un', and 'toep' are supported for resid_corr.")

    # Identify group for residuals
    if (is.null(resid_group)) {
      if (has_random) {
        resid_group_var <- as.character(bars[[1]][[3]])
      } else {
        stop("Residual correlation requires 'resid_group' (e.g., resid_group = 'ID') if no random effects are specified.")
      }
    } else {
      resid_group_var <- resid_group
    }

    if (!(resid_group_var %in% names(mf))) stop(sprintf("Residual grouping variable '%s' not found in data.", resid_group_var))

    group_resid <- as.numeric(as.factor(mf[[resid_group_var]]))
    dat$group_resid <- group_resid
    dat$num_groups_resid <- max(group_resid)
    if (!is.null(resid_time)) {
      dat$time_resid <- as.numeric(as.factor(mf[[resid_time]])) - 1
    } else {
      # Use row order within group
      dat$time_resid <- ave(rep(1, N), group_resid, FUN = seq_along) - 1
    }
    dat$max_T_resid <- max(table(group_resid))
  }

  if (has_random) {
    dat <- c(dat, dat_ranef)
  }
  if (!is.null(sigma_idx)) {
    dat$sigma_idx <- sigma_idx
    dat$num_sigma_groups <- num_sigma_groups
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
      # For Poisson/log-link models, Stan often uses Normal(log(mean(y)), 2.5 or 5) for intercept
      # and Normal(0, 2.5) for coefficients.
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- log(mean(Y) + 0.5))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(alpha_prior_sd <- .(prior$sd_ratio))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- 1.0)
    } else {
      if (prior_type == "jzs") {
         # For JZS, we use defaults if y_range is missing
         if (!is.null(y_range)) {
           setup_exprs[[length(setup_exprs) + 1]] <- quote(half_d_y <- diff(y_range) / 2)
           setup_exprs[[length(setup_exprs) + 1]] <- bquote(base_scale <- half_d_y * .(prior$sd_ratio))
           setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- half_d_y)
           setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mean(y_range))
         } else {
           # JZS defaults for alpha/sigma part if y_range is missing (though we use flat priors anyway)
           setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- 1.0)
           setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- 10.0)
           setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- 0)
         }
         setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate <- 1.0 / base_scale)
         setup_exprs[[length(setup_exprs) + 1]] <- bquote(jzs_r <- .(prior$r))
      } else {
         setup_exprs[[length(setup_exprs) + 1]] <- quote(half_d_y <- diff(y_range) / 2)
         setup_exprs[[length(setup_exprs) + 1]] <- bquote(base_scale <- half_d_y * .(prior$sd_ratio))
         setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- half_d_y)
         setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mean(y_range))
         setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate <- 1.0 / base_scale)
      }
    }

    setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_rate <- 1.0 / base_scale)

    if (K_tmp > 0) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_sd <- apply(X, 2, sd))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(beta_prior_sd <- .(prior$max_beta) * base_scale / X_sd)
      
      if (prior_type == "jzs") {
        # Calculate (X_c' X_c / N)^-1 for multivariate JZS
        # X_c is centered X.
        X_mean_val <- apply(X, 2, mean)
        X_c_val <- X - rep(1, N) %*% t(X_mean_val)
        
        X_for_jzs <- X_c_val
        
        if (ncol(X_for_jzs) > 0) {
          # Calculate (X'X/N)^-1 with numerical safety
          XtX_inv <- solve(t(X_for_jzs) %*% X_for_jzs / N)
          # Force symmetry and add a tiny jitter to ensure PD
          jzs_V_val <- 0.5 * (XtX_inv + t(XtX_inv))
          diag(jzs_V_val) <- diag(jzs_V_val) + 1e-10
          
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(jzs_V <- .(jzs_V_val))
        } else {
          # Null model case: no predictors for JZS
          setup_exprs[[length(setup_exprs) + 1]] <- quote(jzs_V <- matrix(0, 0, 0))
        }
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
  }

  if (use_centering && K_tmp > 0) {
    setup_exprs[[length(setup_exprs) + 1]] <- quote(X_mean <- apply(X, 2, mean))
    if (has_intercept) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_c <- X - rep(1, N) %*% t(X_mean))
    }
  }
  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(dat)
  eval(setup_ast, tmp_env)
  N <- tmp_env$N; K <- tmp_env$K
  num_sigma_groups <- if (!is.null(tmp_env$num_sigma_groups)) tmp_env$num_sigma_groups else 1

  # --- Parameters AST ---
  param_exprs <- list()
  if (has_intercept) {
    if (use_centering) {
      param_exprs[[length(param_exprs) + 1]] <- quote(Intercept_c <- Dim(1))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- quote(Intercept <- Dim(1))
    }
  }
  if (K > 0) {
    if (prior_type == "jzs") {
      param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(K))
    } else if (regularization == "none") {
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

  if (family %in% c("gaussian", "lognormal", "student_t")) param_exprs[[length(param_exprs) + 1]] <- quote(sigma <- Dim(num_sigma_groups, lower = 0))
  if (!is.null(resid_corr)) {
    if (resid_corr %in% c("ar1", "cs")) {
      param_exprs[[length(param_exprs) + 1]] <- quote(rho_resid <- Dim(type = "interval", lower = -0.99, upper = 0.99))
    } else if (resid_corr == "un") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(L_resid <- Dim(c(.(dat$max_T_resid), .(dat$max_T_resid)), type = "CF_corr"))
    } else if (resid_corr == "toep") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(rho_resid <- Dim(.(dat$max_T_resid - 1), type = "interval", lower = -0.99, upper = 0.99))
    }
  }
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
    if (use_centering) {
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
    if (use_centering) {
      if (K > 0) transform_exprs[[1]] <- quote(eta <- Intercept_c + X_c %*% b)
      else transform_exprs[[1]] <- quote(eta <- rep(Intercept_c, N))
    } else {
      if (K > 0) transform_exprs[[1]] <- quote(eta <- Intercept + X %*% b)
      else transform_exprs[[1]] <- quote(eta <- rep(Intercept, N))
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
  if (!is.null(resid_corr)) {
    vi_inner_expr <- switch(resid_corr,
      "ar1" = quote(Vi[r,c] <- s_r * s_c * (rho_resid^abs(ti[r] - ti[c]))),
      "cs"  = quote(if (r == c) Vi[r,c] <- s_r^2 else Vi[r,c] <- s_r * s_c * rho_resid),
      "toep"= quote({
               lag <- abs(ti[r] - ti[c])
               if (lag == 0) Vi[r,c] <- s_r * s_c
               else if (lag <= length(rho_resid)) Vi[r,c] <- s_r * s_c * rho_resid[lag]
               else Vi[r,c] <- 0
            }),
      "un" = quote({
               if (r == c) Vi[r,c] <- s_r^2
               else {
                 idx_r <- ti[r] + 1
                 idx_c <- ti[c] + 1
                 Vi[r,c] <- s_r * s_c * R_mat_resid[idx_r, idx_c]
               }
            })
    )

    # Build model components list
    ll_data_exprs <- list()
    ll_data_exprs[[length(ll_data_exprs) + 1]] <- quote(has_sig_idx <- !is.null(sigma_idx))
    if (resid_corr == "un") {
      ll_data_exprs[[length(ll_data_exprs) + 1]] <- quote(R_mat_resid <- L_resid %*% t(L_resid))
    }

    ll_data_exprs[[length(ll_data_exprs) + 1]] <- bquote(
      for (i in 1:num_groups_resid) {
        idx <- which(group_resid == i)
        ni <- length(idx)
        ti <- time_resid[idx]
        mu_i <- eta[idx]

        Vi <- matrix(0, ni, ni)
        for (r in 1:ni) {
          s_r <- if (has_sig_idx) sigma[sigma_idx[idx[r]]] else sigma
          for (c in 1:ni) {
            s_c <- if (has_sig_idx) sigma[sigma_idx[idx[c]]] else sigma
            .(vi_inner_expr)
          }
        }
        diag(Vi) <- diag(Vi) + 1e-3
        Y[idx] ~ multi_normal(mu_i, Vi)
      }
    )
  } else {
    ll_data_exprs[[1]] <- switch(family,
                                 "gaussian" = if (!is.null(sigma_idx)) quote(Y ~ normal(eta, sigma[sigma_idx])) else quote(Y ~ normal(eta, sigma)),
                                 "lognormal" = if (!is.null(sigma_idx)) quote(Y ~ lognormal(eta, sigma[sigma_idx])) else quote(Y ~ lognormal(eta, sigma)),
                                 "student_t" = if (!is.null(sigma_idx)) quote(Y ~ student_t(nu, eta, sigma[sigma_idx])) else quote(Y ~ student_t(nu, eta, sigma)),
                                 "gamma" = quote(Y ~ gamma(shape, shape / exp(eta))),
                                 "bernoulli" = quote(Y ~ bernoulli_logit(eta)),
                                 "binomial" = quote(Y ~ binomial_logit(trials, eta)),
                                 "poisson" = quote(Y ~ poisson(exp(eta))),
                                 "neg_binomial" = quote(Y ~ neg_binomial_2(exp(eta), phi)),
                                 "ordered" = quote(Y ~ ordered_logistic(eta, cutpoints))
    )
  }

  ll_random_exprs <- list()
  if (has_random) {
    for (b in 1:num_bars) {
      r_re_name <- as.name(paste0("r_re", suffix(b)))
      if (num_ranef_list[[b]] > 1) {
        CF_corr_name <- as.name(paste0("CF_corr", suffix(b)))
        num_ranef_name <- as.name(paste0("num_ranef", suffix(b)))
        num_groups_name <- as.name(paste0("num_groups", suffix(b)))

        if (!is_flat && !is.null(prior$lkj_eta)) {
          ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(.(CF_corr_name) ~ lkj_CF_corr(.(prior$lkj_eta)))
        }
        ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(for (j in 1:.(num_groups_name)) .(r_re_name)[j, ] ~ multi_normal_CF(rep(0, .(num_ranef_name)), rep(1, .(num_ranef_name)), .(CF_corr_name)))
      } else {
        ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(.(r_re_name) ~ normal(0, 1))
      }
    }
  }

  prior_exprs <- list()
  if (prior_type %in% c("weak", "rhs", "ssp", "jzs")) {
    if (family %in% c("gaussian", "lognormal", "student_t")) {
      if (use_weak_info && prior_type == "weak") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(sigma ~ exponential(sigma_rate))
      } else if (prior_type == "jzs") {
        # Jeffreys prior for variance (p(sigma) proportional to 1/sigma)
        if (num_sigma_groups > 1) {
          prior_exprs[[length(prior_exprs) + 1]] <- quote(lp <- lp - sum(log(sigma)))
        } else {
          prior_exprs[[length(prior_exprs) + 1]] <- quote(lp <- lp - log(sigma))
        }
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
      } else if (prior_type == "jzs") {
        # Flat prior for intercept (p(alpha) proportional to 1) - no term added
        prior_exprs[[length(prior_exprs) + 1]] <- quote("# Flat prior for Intercept (JZS)")
      } else if (!is.null(prior$Intercept_sd)) {
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(Intercept_c ~ normal(0, .(prior$Intercept_sd)))
      }
    }

    if (K > 0) {
      if (prior_type == "jzs") {
        # Multivariate JZS (Single-g prior): b ~ Multivariate-Cauchy(0, r^2 * sigma^2 * (X'X/N)^-1)
        # This matches the BayesFactor package and handles correlations correctly.
        # We use bquote to embed the pre-calculated dimension.
        if (!is.null(jzs_V_val) && ncol(jzs_V_val) > 0) {
          K_jzs <- ncol(jzs_V_val)
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(b ~ multi_cauchy(rep(0, .(K_jzs)), (jzs_r * sigma)^2 * jzs_V))
        }
      } else if (regularization == "none") {
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
        } else if (!is_flat && !is.null(prior$sd_rate)) {
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(.(sd_name) ~ exponential(.(prior$sd_rate)))
        }
      }
    }
  } else if (prior_type == "normal") {
    if (family %in% c("gaussian", "lognormal", "student_t") && !is.null(prior$sigma_rate)) {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    }

    if (has_intercept && !is.null(prior$mu_sd)) {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(Intercept_c ~ normal(0, .(prior$mu_sd)))
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

    if (family == "student_t" && !is.null(prior$nu_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(nu ~ exponential(.(prior$nu_rate)))
    if (family == "gamma" && !is.null(prior$shape_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(shape ~ exponential(.(prior$shape_rate)))
    if (family == "neg_binomial" && !is.null(prior$phi_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(phi ~ exponential(.(prior$phi_rate)))
    if (family == "ordered" && !is.null(prior$cutpoint_sd)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))
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
  if (!is.null(generate)) code_obj$generate <- generate

  par_names_list <- list()
  if (num_sigma_groups > 1) {
    par_names_list$sigma <- sigma_group_names
  }
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
      par_names_list[[paste0("sd", suffix(b))]] <- paste0(group_labels_list[[b]], ":", r_names)
      if (num_ranef_list[[b]] > 1) par_names_list[[paste0("corr", suffix(b))]] <- r_names
    }
  }

  view_vars <- c()
  if (has_intercept) view_vars <- c(if (use_centering) "Intercept_c" else "Intercept")
  if (K > 0) view_vars <- c(view_vars, "b")
  view_vars <- c(view_vars, "sigma")
  if (has_random) {
    for (b in 1:num_bars) {
      view_vars <- c(view_vars, paste0("sd", suffix(b)))
      if (num_ranef_list[[b]] > 1) view_vars <- c(view_vars, paste0("corr", suffix(b)))
    }
  }
  if (!is.null(resid_corr)) {
    if (resid_corr %in% c("ar1", "cs", "toep")) view_vars <- c(view_vars, "rho_resid")
    if (resid_corr == "un") view_vars <- c(view_vars, "L_resid")
  }

  ordered_data <- env_to_ordered_list(tmp_env, dat, setup_ast)
  obj <- rtmb_model(data = ordered_data, code = code_obj, par_names = par_names_list, init = init,
                    fixed = fixed,
                    view = if (!is.null(view)) view else view_vars, silent = FALSE)
  obj$formula <- formula
  obj$raw_data <- data
  obj$family <- family
  # Store preferences
  obj$contrasts <- actual_contrasts
  obj$requested_contrasts <- contrasts


  # --- Metadata for classic inference and other methods ---
  obj$type <- if (has_random || !is.null(resid_corr)) {
    if (family %in% c("gaussian", "lognormal", "student_t")) "lmer" else "glmer"
  } else {
    if (family %in% c("gaussian", "lognormal", "student_t")) "lm" else "glm"
  }
  obj$extra <- list(
    source = "wrapper",
    prior_type = prior$type,
    X_assign = X_assign, 
    X_terms = X_terms, 
    X_colnames = fixed_colnames,
    factors = factors,
    within = within
  )

  fixed_effects <- if (K > 0) "b" else character(0)
  if (has_intercept) {
    fixed_effects <- c(if (use_centering) "Intercept_c" else "Intercept", fixed_effects)
  }
  obj$extra$marginal <- fixed_effects

  return(obj)
}
