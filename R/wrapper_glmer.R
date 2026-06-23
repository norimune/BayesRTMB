#' Reconstruct an Observation-Level Random-Effect Design Matrix
#'
#' @description
#' Converts the transposed block random-effect design matrix used internally by
#' `rtmb_glmer()` into an observation-level matrix. This is kept as a small
#' utility for users who want to work from the block matrix returned by
#' `make_glmer_re_terms()`.
#'
#' @param Zt Transposed block random-effect design matrix.
#' @param group_idx Integer group index for each observation.
#' @param N Number of observations. Defaults to `length(group_idx)`.
#'
#' @return A matrix with `N` rows and one column per random-effect coefficient.
#' @export
make_glmer_Z_matrix <- function(Zt, group_idx, N = length(group_idx)) {
  num_groups <- max(group_idx)
  num_ranef <- nrow(Zt) / num_groups
  Z_mat <- matrix(0, nrow = N, ncol = num_ranef)
  for (i in seq_len(N)) {
    g <- group_idx[i]
    Z_mat[i, ] <- Zt[((g - 1) * num_ranef + 1):(g * num_ranef), i]
  }
  Z_mat
}

.deparse_glmer_term <- function(x) {
  paste(deparse(x), collapse = " ")
}

.expand_glmer_bars <- function(bars) {
  expanded_bars <- list()
  if (length(bars) > 0) {
    for (bar in bars) {
      grp_expr <- bar[[3]]
      term_labels <- attr(terms(as.formula(paste("~", .deparse_glmer_term(grp_expr)))), "term.labels")
      for (label in term_labels) {
        new_bar <- bar
        new_bar[[3]] <- parse(text = label)[[1]]
        expanded_bars[[length(expanded_bars) + 1]] <- new_bar
      }
    }
  }
  expanded_bars
}

.make_glmer_Zt <- function(Z_mat, group_idx, num_groups) {
  N <- nrow(Z_mat)
  num_ranef <- ncol(Z_mat)
  Zt <- matrix(0, nrow = num_groups * num_ranef, ncol = N)
  for (i in seq_len(N)) {
    g <- group_idx[i]
    Zt[((g - 1) * num_ranef + 1):(g * num_ranef), i] <- Z_mat[i, ]
  }
  Zt
}

.select_glmer_setup_data <- function(formula, data, resid_group = NULL, resid_time = NULL) {
  has_random <- !is.null(findbars(formula))
  mf_form <- if (has_random) subbars(formula) else formula
  if (!inherits(mf_form, "formula")) mf_form <- as.formula(mf_form, env = environment(formula))
  if (!is.null(resid_group)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_group)))
  if (!is.null(resid_time)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_time)))

  vars <- unique(all.vars(mf_form))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  as.data.frame(data)[, vars, drop = FALSE]
}

.resolve_glmer_cwc_value <- function(expr, env, data_names = character()) {
  if (is.character(expr)) {
    return(expr)
  }
  if (is.symbol(expr)) {
    name <- as.character(expr)
    if (name %in% data_names) {
      return(name)
    }
    value <- tryCatch(eval(expr, envir = env), error = function(e) NULL)
    if (is.null(value)) name else value
  } else {
    eval(expr, envir = env)
  }
}

.resolve_glmer_cwc_spec <- function(expr, env, data_names = character()) {
  if (identical(expr, quote(NULL))) {
    return(NULL)
  }

  spec <- if (is.call(expr) && identical(expr[[1]], as.name("list"))) {
    elems <- as.list(expr)[-1]
    values <- lapply(elems, .resolve_glmer_cwc_value, env = env, data_names = data_names)
    names(values) <- names(elems)
    values
  } else {
    tryCatch(
      eval(expr, envir = env),
      error = function(e) {
        stop(
          "Invalid 'cwc' format. Use list(cluster = ID, pars = c('var1', 'var2')) ",
          "or list(ID, 'var1').",
          call. = FALSE
        )
      }
    )
  }

  if (!is.null(spec) && !is.list(spec)) {
    stop(
      "Invalid 'cwc' format. Use list(cluster = ID, pars = c('var1', 'var2')) ",
      "or list(ID, 'var1').",
      call. = FALSE
    )
  }
  spec
}

.glmer_fixed_numeric_vars <- function(formula, data, exclude = character()) {
  vars_in_fixed <- all.vars(nobars(formula))
  if (length(vars_in_fixed) > 0L) {
    vars_in_fixed <- vars_in_fixed[-1]
  }
  vars_in_fixed <- unique(vars_in_fixed)
  vars_in_fixed <- setdiff(vars_in_fixed, exclude)
  vars_in_fixed <- vars_in_fixed[vars_in_fixed %in% names(data)]
  vars_in_fixed[vapply(data[vars_in_fixed], is.numeric, logical(1))]
}

#' Prepare GLMM Formula Components
#'
#' @description
#' Builds the response, fixed-effect model matrix, and lme4-style random-effect
#' design components used by `rtmb_glmer()`. This is the user-facing helper that
#' reproduces what the wrapper places in the generated `setup` block.
#'
#' @param formula lme4-style formula.
#' @param data Data frame.
#' @param family Character string of the distribution family.
#' @param resid_group Optional residual-correlation grouping variable.
#' @param resid_time Optional residual-correlation time variable.
#' @param within Optional list for wide-to-long conversion when the response is
#'   written as `cbind(...)`.
#' @param factors Optional character vector of variables to treat as factors.
#' @param missing Missing value handling strategy: "listwise".
#'
#' @return A list containing `Y`, `X`, `trials`, `offset`, `N`, fixed-effect
#' metadata, and random-effect terms.
#' @export
make_glmer_re_terms <- function(formula, data, family = "gaussian",
                                resid_group = NULL, resid_time = NULL,
                                within = NULL, factors = NULL,
                                missing = "listwise") {
  # Expand dot (.) in formula if present
  if (!is.null(data) && inherits(formula, "formula")) {
    terms_expanded <- try(terms(formula, data = data), silent = TRUE)
    if (!inherits(terms_expanded, "try-error")) {
      formula <- formula(terms_expanded)
    }
  }

  if (family == "gaussian" && is.call(formula[[2]]) && identical(formula[[2]][[1]], as.name("cbind"))) {
    processed <- .handle_wide_to_long(formula, data, within, factors)
    formula <- processed$formula
    data <- processed$data
    factors <- processed$factors
  }

  if (!is.null(factors)) {
    data <- as.data.frame(data)
    for (f in factors) {
      if (f %in% names(data)) data[[f]] <- as.factor(data[[f]])
    }
  }

  data <- .select_glmer_setup_data(formula, data, resid_group = resid_group, resid_time = resid_time)
  has_random <- !is.null(findbars(formula))
  fixed_form <- if (has_random) nobars(formula) else formula
  mf_form <- if (has_random) subbars(formula) else formula
  if (!inherits(fixed_form, "formula")) fixed_form <- as.formula(fixed_form, env = environment(formula))
  if (!inherits(mf_form, "formula")) mf_form <- as.formula(mf_form, env = environment(formula))

  if (!is.null(resid_group)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_group)))
  if (!is.null(resid_time)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_time)))

  na_act <- if (missing == "listwise") na.omit else na.pass
  mf <- model.frame(mf_form, data = data, na.action = na_act)
  Y <- model.response(mf)
  X_full <- model.matrix(fixed_form, mf)
  X_assign <- attr(X_full, "assign")
  X_terms <- attr(terms(fixed_form), "term.labels")
  offset <- model.offset(mf)
  N <- nrow(mf)

  num_categories <- NULL
  if (is.matrix(Y)) {
    if (ncol(Y) == 2 && family == "binomial") {
      trials <- as.numeric(Y[, 1] + Y[, 2])
      Y <- as.numeric(Y[, 1])
    } else if (ncol(Y) == 1) {
      Y <- as.numeric(Y[, 1])
      trials <- rep(1, length(Y))
    } else {
      stop("Invalid matrix format for the response variable.")
    }
  } else {
    trials <- rep(1, length(Y))
    if (is.factor(Y)) {
      Y <- if (family %in% c("bernoulli", "binomial")) as.numeric(Y) - 1 else as.numeric(Y)
    } else {
      Y <- as.numeric(Y)
    }
  }

  if (family %in% c("ordered", "sequential")) {
    if (any(!is.finite(Y)) || any(abs(Y - round(Y)) > .Machine$double.eps^0.5) || any(Y < 1)) {
      stop(
        "Ordinal responses must be category codes 1, 2, ..., K or an ordered/factor response. ",
        "Check that the response variable is not being centered or otherwise transformed.",
        call. = FALSE
      )
    }
    Y <- as.integer(round(Y))
    num_categories <- max(Y)
  }

  has_intercept <- "(Intercept)" %in% colnames(X_full)
  if (family %in% c("ordered", "sequential")) has_intercept <- FALSE
  fixed_colnames <- colnames(X_full)
  X <- X_full
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  bars <- if (has_random) .expand_glmer_bars(findbars(formula)) else list()
  num_bars <- length(bars)
  suffix <- function(b) if (num_bars > 1) paste0("_", b) else ""

  random_terms <- list()
  random_data <- list()
  num_ranef_list <- list()
  num_groups_list <- list()
  ranef_names_list <- list()
  group_labels_list <- list()

  if (num_bars > 0) {
    for (b in seq_len(num_bars)) {
      bar <- bars[[b]]
      re_form <- as.formula(paste("~", .deparse_glmer_term(bar[[2]])))
      Z_mat <- model.matrix(re_form, mf)

      group_var_label <- .deparse_glmer_term(bar[[3]])
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
      num_ranef <- ncol(Z_mat)
      Zt <- .make_glmer_Zt(Z_mat, group_idx, num_groups)

      random_data[[paste0("Zt", suffix(b))]] <- Zt
      random_data[[paste0("group_idx", suffix(b))]] <- group_idx
      num_ranef_list[[b]] <- num_ranef
      num_groups_list[[b]] <- num_groups

      r_names <- colnames(Z_mat)
      r_names[r_names == "(Intercept)"] <- "Int"
      ranef_names_list[[b]] <- r_names
      group_labels_list[[b]] <- group_var_label

      random_terms[[b]] <- list(
        Z = Z_mat,
        Zt = Zt,
        group_idx = group_idx,
        group = flist,
        group_label = group_var_label,
        num_groups = num_groups,
        num_ranef = num_ranef,
        ranef_names = r_names
      )
    }
  }

  list(
    Y = Y,
    X = X,
    trials = trials,
    offset = offset,
    N = N,
    mf = mf,
    X_assign = X_assign,
    X_terms = X_terms,
    fixed_colnames = fixed_colnames,
    fixed_names = colnames(X),
    has_intercept = has_intercept,
    num_categories = num_categories,
    random = list(
      bars = bars,
      num_bars = num_bars,
      terms = random_terms,
      data = random_data,
      num_ranef_list = num_ranef_list,
      num_groups_list = num_groups_list,
      ranef_names_list = ranef_names_list,
      group_labels_list = group_labels_list
    )
  )
}

#' RTMB-based GLMM wrapper function
#'
#' @param formula lme4-style formula (e.g., Y ~ X + (1 | GID))
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian",
#'   "binomial", "poisson", "ordered", "sequential")
#' @param laplace Logical; whether to marginalize random effects using Laplace approximation
#' @param prior An object of class "rtmb_prior" specifying the prior
#'   distribution. Use `prior_flat()`, `prior_normal()`, `prior_weak()`,
#'   `prior_rhs()`, or `prior_ssp()`. `prior_jzs()` is available for
#'   continuous families. Default is `prior_flat()`.
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Required when using weakly informative or regularized priors with continuous models.
#' @param init List of initial values (generated automatically based on glm if omitted)

#' @param gmc Character vector of variable names for Grand Mean Centering (GMC). If "all", all numeric variables are centered.
#' @param centering Alias for `gmc`.
#' @param cwc List for Centering Within Cluster (CWC). Should contain \code{cluster} (group variable) and \code{pars} (variable names to center).
#'   You can also use \code{cwc = list(ID, "x")} or \code{cwc = list(ID, "all")};
#'   \code{"all"} centers all numeric fixed-effect variables within the cluster.
#' @param view Character vector of parameter names to prioritize in summary.
#' @param factors Character vector of variable names to be treated as factors.
#' @param contrasts Character string specifying the contrast type ("treatment" or "sum").
#' @param sigma_by Character vector specifying variables to group residual variance by (heteroscedasticity).
#' @param resid_corr Residual correlation structure: "ar1" (Autoregressive), "cs" (Compound Symmetry), "toep" (Toeplitz), or "un" (Unstructured).
#' @param resid_time Variable name for time points in residual correlation.
#' @param resid_group Variable name for grouping in residual correlation.
#' @param within Optional list for wide-to-long conversion.
#' @param generate Optional expression for generated quantities.
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param .force_sum Logical; internal use only.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param missing Missing value handling strategy: "listwise".
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glmer <- function(formula, data, family = "gaussian", laplace = FALSE,
                       prior = prior_flat(),
                       y_range = NULL,
                       init = NULL,
                       fixed = NULL,
                       gmc = NULL,
                       centering = NULL,
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
                       missing = c("listwise", "fiml"),
                       WAIC = FALSE,
                       .force_sum = FALSE) {

  cwc <- if (base::missing(cwc)) NULL else .resolve_glmer_cwc_spec(substitute(cwc), parent.frame(), names(data))

  if (!is.null(centering)) {
    if (!is.null(gmc) && !identical(gmc, centering)) {
      stop("Specify only one of 'gmc' or 'centering', or use identical values.", call. = FALSE)
    }
    gmc <- centering
  }

  valid_families <- c(
    "gaussian", "lognormal", "student_t", "bernoulli", "binomial",
    "poisson", "neg_binomial", "gamma", "ordered", "sequential"
  )
  if (!is.character(family) || length(family) != 1L || is.na(family) ||
      !(family %in% valid_families)) {
    stop(
      sprintf(
        "Invalid 'family' value. Valid options are: %s",
        paste(valid_families, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!is.null(contrasts) &&
      (!is.character(contrasts) || length(contrasts) != 1L ||
       !(contrasts %in% c("treatment", "sum")))) {
    stop("'contrasts' must be either 'treatment', 'sum', or NULL.", call. = FALSE)
  }

  # Expand dot (.) in formula if present
  if (!is.null(data) && inherits(formula, "formula")) {
    terms_expanded <- try(terms(formula, data = data), silent = TRUE)
    if (!inherits(terms_expanded, "try-error")) {
      formula <- formula(terms_expanded)
    }
  }

  missing <- match.arg(missing)
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
        .glmer_fixed_numeric_vars(formula, data_centered)
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
      group_var_name <- if (is.character(cluster_var) && length(cluster_var) == 1L &&
                            cluster_var %in% names(data_centered)) {
        cluster_var
      } else {
        NULL
      }

      group_id <- if (!is.null(group_var_name)) {
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

      if (is.character(target_pars) && length(target_pars) == 1L && identical(target_pars, "all")) {
        target_pars <- .glmer_fixed_numeric_vars(
          formula,
          data_centered,
          exclude = group_var_name %||% character()
        )
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

  prior <- .validate_prior_type(
    prior,
    allowed = c("flat", "normal", "weak", "rhs", "ssp", "jzs"),
    context = "rtmb_glmer()"
  )

  prior_type <- prior$type
  if (identical(prior_type, "jzs") &&
      !(family %in% c("gaussian", "lognormal", "student_t"))) {
    stop(
      "prior_jzs() is currently supported only for continuous rtmb_glmer() families ",
      "('gaussian', 'lognormal', and 'student_t').",
      call. = FALSE
    )
  }

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

  glmer_terms <- make_glmer_re_terms(
    formula = formula,
    data = data,
    family = family,
    resid_group = resid_group,
    resid_time = resid_time,
    within = within,
    factors = factors
  )
  Y <- glmer_terms$Y
  X <- glmer_terms$X
  trials <- glmer_terms$trials
  X_assign <- glmer_terms$X_assign
  X_terms <- glmer_terms$X_terms
  offset <- glmer_terms$offset
  N <- glmer_terms$N
  mf <- glmer_terms$mf
  fixed_colnames <- glmer_terms$fixed_colnames
  fixed_names <- glmer_terms$fixed_names
  has_intercept <- glmer_terms$has_intercept
  ordered_num_categories <- glmer_terms$num_categories
  bars <- glmer_terms$random$bars
  num_bars <- glmer_terms$random$num_bars
  suffix <- function(b) if (num_bars > 1) paste0("_", b) else ""
  num_ranef_list <- glmer_terms$random$num_ranef_list
  num_groups_list <- glmer_terms$random$num_groups_list
  ranef_names_list <- glmer_terms$random$ranef_names_list
  group_labels_list <- glmer_terms$random$group_labels_list
  
  jzs_V_val <- NULL
  K_tmp <- ncol(X)
  sequential_transitions_tmp <- if (identical(family, "sequential")) ordered_num_categories - 1 else 1

  # --- Robust SE Weights ---
  dat_weights <- list(robust_obs_weight = rep(1, N))
  if (has_random && num_bars > 0) {
    for (b in 1:num_bars) {
      weight_name <- paste0("robust_re_weight", suffix(b))
      dat_weights[[weight_name]] <- rep(1, num_groups_list[[b]])
    }
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

    b_est <- con_est_list$b
    if (is.null(b_est) || length(b_est) == 0 || any(is.na(as.numeric(b_est)))) {
      b_est <- if (identical(family, "sequential")) matrix(0, K_tmp, sequential_transitions_tmp) else rep(0, K_tmp)
    } else if (identical(family, "sequential")) {
      b_est <- matrix(as.numeric(b_est), K_tmp, sequential_transitions_tmp)
    } else {
      b_est <- as.numeric(b_est)
    }

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
        b_init <- if (has_intercept) init_coef[-1] else init_coef
        if (family == "sequential") {
          init$b <- matrix(b_init, K_tmp, sequential_transitions_tmp)
        } else {
          init$b <- b_init
        }
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

    }, error = function(e) init <<- NULL)
  }

  setup_df <- .select_glmer_setup_data(formula, data, resid_group = resid_group, resid_time = resid_time)
  class(setup_df) <- c("rtmb_setup_df", class(setup_df))
  dat <- list(df = setup_df, formula = formula,
              sigma_idx = sigma_idx, num_sigma_groups = num_sigma_groups)
  if (has_random && !identical(family, "gaussian")) dat$family_name <- family
  if (has_random && !is.null(resid_group)) dat$resid_group_name <- resid_group
  if (has_random && !is.null(resid_time)) dat$resid_time_name <- resid_time
  if (has_random && !is.null(within)) dat$within_info <- within
  if (has_random && !is.null(factors)) dat$factors_info <- factors

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

  if (!is.null(sigma_idx)) {
    dat$sigma_idx <- sigma_idx
    dat$num_sigma_groups <- num_sigma_groups
  }

  # --- Setup AST ---
  setup_exprs <- list()
  Y_setup_raw <- model.response(mf)
  
  if (is.character(Y_setup_raw)) {
    stop("The response variable must be numeric (or a factor for classification models). Character variables are not supported.", call. = FALSE)
  }
  
  if (is.factor(Y_setup_raw) && !(family %in% c("binomial", "bernoulli", "ordered", "sequential"))) {
    stop(sprintf("The response variable for family '%s' must be numeric. Factor variables are not supported.", family), call. = FALSE)
  }
  if (has_random) {
    glmer_re_call_args <- list(formula = as.name("formula"), data = as.name("df"))
    if (!identical(family, "gaussian")) glmer_re_call_args$family <- as.name("family_name")
    if (!is.null(resid_group)) glmer_re_call_args$resid_group <- as.name("resid_group_name")
    if (!is.null(resid_time)) glmer_re_call_args$resid_time <- as.name("resid_time_name")
    if (!is.null(within)) glmer_re_call_args$within <- as.name("within_info")
    if (!is.null(factors)) glmer_re_call_args$factors <- as.name("factors_info")
    glmer_re_call_args$missing <- missing
    setup_exprs[[length(setup_exprs) + 1]] <- as.call(list(as.name("<-"), as.name("res"), as.call(c(list(as.name("make_glmer_re_terms")), glmer_re_call_args))))
    setup_exprs[[length(setup_exprs) + 1]] <- quote(Y <- res$Y)
    if (family == "binomial") {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(trials <- res$trials)
    }
    setup_exprs[[length(setup_exprs) + 1]] <- quote(X <- res$X)
    if (!is.null(offset)) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(offset <- res$offset)
    }
    setup_exprs[[length(setup_exprs) + 1]] <- quote(N <- res$N)
    setup_exprs[[length(setup_exprs) + 1]] <- quote(K <- ncol(X))
    if (family %in% c("ordered", "sequential")) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(num_categories <- res$num_categories)
    }
    setup_exprs[[length(setup_exprs) + 1]] <- "# Random-effect design"
    for (b in 1:num_bars) {
      group_idx_name <- as.name(paste0("group_idx", suffix(b)))
      Z_mat_name <- as.name(paste0("Z_mat", suffix(b)))
      num_groups_name <- as.name(paste0("num_groups", suffix(b)))
      num_ranef_name <- as.name(paste0("num_ranef", suffix(b)))

      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(group_idx_name) <- res$random$terms[[.(b)]]$group_idx)
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(Z_mat_name) <- res$random$terms[[.(b)]]$Z)
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(num_groups_name) <- res$random$terms[[.(b)]]$num_groups)
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(num_ranef_name) <- res$random$terms[[.(b)]]$num_ranef)
    }
  } else {
    na_action_expr <- if (missing == "listwise") quote(na.omit) else quote(na.pass)
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(na_act <- .(na_action_expr))
    setup_exprs[[length(setup_exprs) + 1]] <- quote(mf <- model.frame(formula, df, na.action = na_act))
    setup_exprs[[length(setup_exprs) + 1]] <- quote(Y <- model.response(mf))
    if (is.matrix(Y_setup_raw)) {
      if (ncol(Y_setup_raw) == 2 && family == "binomial") {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(trials <- as.numeric(Y[, 1] + Y[, 2]))
        setup_exprs[[length(setup_exprs) + 1]] <- quote(Y <- as.numeric(Y[, 1]))
      } else if (ncol(Y_setup_raw) == 1) {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(Y <- as.numeric(Y[, 1]))
        if (family == "binomial") {
          setup_exprs[[length(setup_exprs) + 1]] <- quote(trials <- rep(1, length(Y)))
        }
      }
    } else {
      if (is.factor(Y_setup_raw)) {
        if (family %in% c("bernoulli", "binomial")) {
          setup_exprs[[length(setup_exprs) + 1]] <- quote(Y <- as.numeric(Y) - 1)
        } else {
          setup_exprs[[length(setup_exprs) + 1]] <- quote(Y <- as.numeric(Y))
        }
      } else if (!is.numeric(Y_setup_raw)) {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(Y <- as.numeric(Y))
      }
      if (family == "binomial") {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(trials <- rep(1, length(Y)))
      }
    }
    if (family %in% c("ordered", "sequential")) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(Y <- as.integer(round(Y)))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(num_categories <- max(Y))
    }
    if ("(Intercept)" %in% fixed_colnames) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_full <- model.matrix(formula, mf))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X <- X_full[, colnames(X_full) != "(Intercept)", drop = FALSE])
    } else {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X <- model.matrix(formula, mf))
    }
    if (!is.null(offset)) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(offset <- model.offset(mf))
    }
    setup_exprs[[length(setup_exprs) + 1]] <- quote(N <- length(Y))
    setup_exprs[[length(setup_exprs) + 1]] <- quote(K <- ncol(X))
  }

  if (use_weak_info) {
    setup_exprs[[length(setup_exprs) + 1]] <- "# Prior setting"
    if (family %in% c("bernoulli", "binomial", "ordered", "sequential")) {
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
      if (family == "sequential") {
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(beta_prior_sd <- matrix(.(prior$max_beta) * base_scale / X_sd, K, num_categories - 1))
      } else {
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(beta_prior_sd <- .(prior$max_beta) * base_scale / X_sd)
      }
      
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
        if (family == "sequential") {
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(p0 <- min(.(prior$expected_vars), K * (num_categories - 1) - 1))
          setup_exprs[[length(setup_exprs) + 1]] <- quote(K_eff <- K * (num_categories - 1))
        } else {
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(p0 <- min(.(prior$expected_vars), K - 1))
          setup_exprs[[length(setup_exprs) + 1]] <- quote(K_eff <- K)
        }
        setup_exprs[[length(setup_exprs) + 1]] <- quote(if (p0 < 1) p0 <- 1)
        setup_exprs[[length(setup_exprs) + 1]] <- quote(tau0 <- (p0 / (K_eff - p0)) * (base_scale / sqrt(N)))
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(half_slab_df <- .(prior$slab_df) / 2.0)
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(half_slab_scale2 <- (.(prior$slab_df) * .(prior$slab_scale)^2) / 2.0)
      } else if (regularization == "ssp") {
        if (family == "sequential") {
          setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_scale <- matrix(base_scale / X_sd, K, num_categories - 1))
        } else {
          setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_scale <- base_scale / X_sd)
        }
      }
    }
  }

  if (use_centering && K_tmp > 0) {
    setup_exprs[[length(setup_exprs) + 1]] <- "# Centering"
    setup_exprs[[length(setup_exprs) + 1]] <- quote(X_mean <- apply(X, 2, mean))
    if (has_intercept) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_c <- X - rep(1, N) %*% t(X_mean))
    }
  }
  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(dat)
  eval(setup_ast, tmp_env)
  N <- tmp_env$N; K <- tmp_env$K
  num_categories <- tmp_env$num_categories
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
      if (family == "sequential") param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(c(K, num_categories - 1)))
      else param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(K))
    } else if (regularization == "none") {
      if (family == "sequential") param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(c(K, num_categories - 1)))
      else param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(K))
    } else if (regularization == "rhs") {
      if (family == "sequential") {
        param_exprs[[length(param_exprs) + 1]] <- quote(z <- Dim(c(K, num_categories - 1)))
        param_exprs[[length(param_exprs) + 1]] <- quote(lambda <- Dim(c(K, num_categories - 1), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- quote(w_lambda <- Dim(c(K, num_categories - 1), lower = 0))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- quote(z <- Dim(K))
        param_exprs[[length(param_exprs) + 1]] <- quote(lambda <- Dim(K, lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- quote(w_lambda <- Dim(K, lower = 0))
      }
      param_exprs[[length(param_exprs) + 1]] <- quote(tau_hs <- Dim(1, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(w_tau <- Dim(1, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(c2 <- Dim(1, lower = 0))
    } else if (regularization == "ssp") {
      if (family == "sequential") {
        param_exprs[[length(param_exprs) + 1]] <- quote(beta_raw <- Dim(c(K, num_categories - 1)))
        param_exprs[[length(param_exprs) + 1]] <- quote(r <- Dim(c(K, num_categories - 1), lower = 0.001, upper = 0.999))
        param_exprs[[length(param_exprs) + 1]] <- quote(tau <- Dim(c(K, num_categories - 1), lower = 0))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- quote(beta_raw <- Dim(K))
        param_exprs[[length(param_exprs) + 1]] <- quote(r <- Dim(K, lower = 0.001, upper = 0.999))
        param_exprs[[length(param_exprs) + 1]] <- quote(tau <- Dim(K, lower = 0))
      }
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
  if (family == "ordered") {
    param_exprs[[length(param_exprs) + 1]] <- quote(cutpoints <- Dim(num_categories - 1, type = "ordered"))
  }
  if (family == "sequential") {
    param_exprs[[length(param_exprs) + 1]] <- quote(cutpoints <- Dim(num_categories - 1))
  }
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
  if (family == "sequential") {
    if (K > 0) transform_exprs[[1]] <- quote(eta <- X %*% b)
    else transform_exprs[[1]] <- quote(eta <- rtmb_array(0, dim = c(N, num_categories - 1)))
  } else if (has_intercept) {
    if (use_centering) {
      if (K > 0) transform_exprs[[1]] <- quote(eta <- Intercept_c + X_c %*% b)
      else transform_exprs[[1]] <- quote(eta <- rep(Intercept_c, N))
    } else {
      if (K > 0) transform_exprs[[1]] <- quote(eta <- Intercept + X %*% b)
      else transform_exprs[[1]] <- quote(eta <- rep(Intercept, N))
    }
  } else {
    if (K > 0) transform_exprs[[1]] <- quote(eta <- X %*% b)
    else transform_exprs[[1]] <- quote(eta <- rtmb_vector(0, N))
  }

  if (has_random) {
    for (b in 1:num_bars) {
      Z_mat_name <- as.name(paste0("Z_mat", suffix(b)))
      r_re_name <- as.name(paste0("r_re", suffix(b)))
      group_idx_name <- as.name(paste0("group_idx", suffix(b)))
      sd_name <- as.name(paste0("sd", suffix(b)))

      if (num_ranef_list[[b]] > 1) {
        if (family == "sequential") {
          transform_exprs[[length(transform_exprs) + 1]] <- bquote(
            for (i in 1:N) {
              for (k_seq in 1:(num_categories - 1)) {
                eta[i, k_seq] <- eta[i, k_seq] + sum(.(Z_mat_name)[i, ] * .(r_re_name)[.(group_idx_name)[i], ] * .(sd_name))
              }
            }
          )
        } else {
          transform_exprs[[length(transform_exprs) + 1]] <- bquote(
            for (i in 1:N) {
              eta[i] <- eta[i] + sum(.(Z_mat_name)[i, ] * .(r_re_name)[.(group_idx_name)[i], ] * .(sd_name))
            }
          )
        }
      } else {
        if (family == "sequential") {
          transform_exprs[[length(transform_exprs) + 1]] <- bquote(
            for (k_seq in 1:(num_categories - 1)) {
              eta[, k_seq] <- eta[, k_seq] + .(Z_mat_name)[, 1] * .(r_re_name)[.(group_idx_name)] * .(sd_name)
            }
          )
        } else {
          transform_exprs[[length(transform_exprs) + 1]] <- bquote(eta <- eta + .(Z_mat_name)[, 1] * .(r_re_name)[.(group_idx_name)] * .(sd_name))
        }
      }
    }
  }
  if (!is.null(offset)) {
    if (family == "sequential") {
      transform_exprs[[length(transform_exprs) + 1]] <- quote(
        for (k_seq in 1:(num_categories - 1)) {
          eta[, k_seq] <- eta[, k_seq] + offset
        }
      )
    } else {
      transform_exprs[[length(transform_exprs) + 1]] <- quote(eta <- eta + offset)
    }
  }

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

        Vi <- rtmb_array(0, dim = c(ni, ni))
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
                                 "ordered" = quote(Y ~ ordered_logistic(eta, cutpoints)),
                                 "sequential" = quote(Y ~ sequential_logistic(eta, cutpoints))
    )
  }

  ll_random_exprs <- list()
  if (has_random) {
    for (b in 1:num_bars) {
        num_groups_val <- num_groups_list[[b]]
        num_ranef_val <- num_ranef_list[[b]]
        if (num_ranef_val == 1) {
          r_re_name <- as.name(paste0("r_re", suffix(b)))
          ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(.(r_re_name) ~ normal(0, 1))
        } else {
          r_re_name <- as.name(paste0("r_re", suffix(b)))
          CF_corr_name <- as.name(paste0("CF_corr", suffix(b)))
          if (!is_flat && !is.null(prior$lkj_eta)) {
            ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(.(CF_corr_name) ~ lkj_CF_corr(.(prior$lkj_eta)))
          }
          ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(for (j in 1:.(num_groups_val)) .(r_re_name)[j, ] ~ multi_normal_CF(rep(0, .(num_ranef_val)), rep(1, .(num_ranef_val)), .(CF_corr_name)))
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
    if (family %in% c("ordered", "sequential") && !is.null(prior$cutpoint_sd)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))

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
    if (family %in% c("ordered", "sequential") && !is.null(prior$cutpoint_sd)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))
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

  waic_ast <- NULL
  if (isTRUE(WAIC)) {
    if (!is.null(resid_corr)) {
      waic_exprs <- c(transform_exprs, list(bquote({
        log_lik <- rtmb_vector(0, num_groups_resid)
        has_sig_idx <- !is.null(sigma_idx)
        .(if (resid_corr == "un") quote(R_mat_resid <- L_resid %*% t(L_resid)) else quote(NULL))
        for (i in 1:num_groups_resid) {
          idx <- which(group_resid == i)
          ni <- length(idx)
          ti <- time_resid[idx]
          mu_i <- eta[idx]
          Vi <- rtmb_array(0, dim = c(ni, ni))
          for (r in 1:ni) {
            s_r <- if (has_sig_idx) sigma[sigma_idx[idx[r]]] else sigma
            for (c in 1:ni) {
              s_c <- if (has_sig_idx) sigma[sigma_idx[idx[c]]] else sigma
              .(vi_inner_expr)
            }
          }
          diag(Vi) <- diag(Vi) + 1e-3
          log_lik[i] <- multi_normal_lpdf(Y[idx], mu_i, Vi)
        }
      })))
      waic_ast <- as.call(c(list(as.name("{")), waic_exprs))
    } else {
      waic_ll <- switch(family,
        "gaussian" = if (!is.null(sigma_idx)) quote(log_lik <- normal_lpdf(Y, eta, sigma[sigma_idx], sum = FALSE)) else quote(log_lik <- normal_lpdf(Y, eta, sigma, sum = FALSE)),
        "lognormal" = if (!is.null(sigma_idx)) quote(log_lik <- lognormal_lpdf(Y, eta, sigma[sigma_idx], sum = FALSE)) else quote(log_lik <- lognormal_lpdf(Y, eta, sigma, sum = FALSE)),
        "student_t" = if (!is.null(sigma_idx)) quote(log_lik <- student_t_lpdf(Y, nu, eta, sigma[sigma_idx], sum = FALSE)) else quote(log_lik <- student_t_lpdf(Y, nu, eta, sigma, sum = FALSE)),
        "gamma" = quote(log_lik <- gamma_lpdf(Y, shape, shape / exp(eta), sum = FALSE)),
        "bernoulli" = quote(log_lik <- bernoulli_logit_lpmf(Y, eta, sum = FALSE)),
        "binomial" = quote(log_lik <- binomial_logit_lpmf(Y, trials, eta, sum = FALSE)),
        "poisson" = quote(log_lik <- poisson_lpmf(Y, exp(eta), sum = FALSE)),
        "neg_binomial" = quote(log_lik <- neg_binomial_2_lpmf(Y, exp(eta), phi, sum = FALSE)),
        "ordered" = quote(log_lik <- ordered_logistic_lpmf(Y, eta, cutpoints, sum = FALSE)),
        "sequential" = quote(log_lik <- sequential_logistic_lpmf(Y, eta, cutpoints, sum = FALSE))
      )
      waic_ast <- as.call(c(list(as.name("{")), c(transform_exprs, list(waic_ll))))
    }
  }

  code_obj <- list(setup = setup_ast, parameters = param_ast)
  code_obj$setup_env <- .rtmb_setup_env(environment(), setup_ast, exclude = names(dat))
  if (!is.null(tran_ast)) code_obj$transform <- tran_ast
  code_obj$model <- model_ast
  generate <- .rtmb_waic_generate_ast(generate, waic_ast)
  if (!is.null(generate)) code_obj$generate <- generate

  par_names_list <- list()
  if (num_sigma_groups > 1) {
    par_names_list$sigma <- sigma_group_names
  }
  if (K > 0) {
    if (family == "sequential") {
      transition_names <- paste0(seq_len(num_categories - 1), "->", 2:num_categories)
      par_names_list$b <- list(fixed_names, transition_names)
    } else {
      par_names_list$b <- fixed_names
    }
    if (regularization == "rhs") {
      if (family == "sequential") {
        par_names_list$z <- list(fixed_names, transition_names)
        par_names_list$lambda <- list(fixed_names, transition_names)
        par_names_list$w_lambda <- list(fixed_names, transition_names)
      } else {
        par_names_list$z <- fixed_names; par_names_list$lambda <- fixed_names; par_names_list$w_lambda <- fixed_names
      }
    } else if (regularization == "ssp") {
      if (family == "sequential") {
        par_names_list$beta_raw <- list(fixed_names, transition_names)
        par_names_list$r <- list(fixed_names, transition_names)
        par_names_list$tau <- list(fixed_names, transition_names)
      } else {
        par_names_list$beta_raw <- fixed_names; par_names_list$r <- fixed_names; par_names_list$tau <- fixed_names
      }
    }
  }
  if (family %in% c("ordered", "sequential")) {
    par_names_list$cutpoints <- paste0("C", seq_len(num_categories - 1))
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
  if (has_intercept) {
    view_vars <- c("Intercept")
  }
  if (K > 0) view_vars <- c(view_vars, "b")
  if (family %in% c("ordered", "sequential")) view_vars <- c(view_vars, "cutpoints")
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
  ordered_data <- utils::modifyList(ordered_data, dat_weights)
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

  # --- mixed-model metadata for robust standard errors ---
  if (obj$type %in% c("lmer", "glmer")) {
    obj$extra$glmer_info <- list(
      family = family,
      num_bars = num_bars,
      group_labels = group_labels_list,
      weight_names = if (num_bars > 0) paste0("robust_re_weight", sapply(1:num_bars, suffix)) else character(0),
      obs_weight_name = "robust_obs_weight",
      group_idx_names = if (num_bars > 0) paste0("group_idx", sapply(1:num_bars, suffix)) else character(0),
      num_groups_list = num_groups_list
    )
  }

  return(obj)
}
