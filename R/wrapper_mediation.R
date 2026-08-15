#' RTMB-based Mediation Analysis Wrapper
#'
#' @description
#' `rtmb_mediation` performs mediation analysis by simultaneously estimating multiple
#' GLM regression equations. It automatically identifies mediation paths and calculates
#' indirect, direct, and total effects.
#'
#' @param formula A list of formulas defining the regression paths (e.g.,
#'   `list(M ~ X, Y ~ X + M)`). Each equation may optionally include one
#'   random intercept, such as `(1 | ID)`.
#' @param data A data frame containing the variables.
#' @param family A single character string or a list of character strings specifying the error distribution
#'   for each equation (e.g., `family = list("gaussian", "binomial")`). Default is "gaussian".
#' @param prior An object of class "rtmb_prior" specifying the prior distribution.
#' @param y_range Theoretical minimum and maximum values of the response variable.
#' @param fixed A named list of parameter values to fix (optional).
#' @param view Character vector of parameter names to prioritize in summary.
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param gmc Character vector naming predictors to grand-mean center, or
#'   `"all"` to center all numeric predictors used by the mediation equations.
#' @param centering Alias for `gmc`.
#' @param cwc Centering-within-cluster specification. Use, for example,
#'   `list(cluster = ID, pars = c("X", "M"))` or `list(ID, "X")`.
#'   Cluster means are not added automatically.
#' @param ... Reserved; unused arguments are rejected.
#'
#' @details
#' The function identifies mediation paths by looking for variables that are
#' responses in one equation and predictors in another. Indirect effects are
#' calculated as the product of coefficients along these paths (\eqn{a * b}).
#'
#' Random intercepts may be included in all equations or in only a subset of
#' equations. When more than one equation contains a random intercept, all
#' random intercepts must use the same grouping variable and their correlations
#' are estimated jointly. Random slopes, multiple random-effect terms within an
#' equation, and different grouping variables across equations are not yet
#' supported.
#' The mediation-specific classical bootstrap is currently unavailable for
#' random-intercept and CWC models.
#'
#' `gmc` (or its alias `centering`) and `cwc` are applied to predictor uses of
#' the selected variables before each equation's model matrix is constructed.
#' Their response uses remain on the original scale. Thus, if `M` is the
#' response in one equation and a predictor in another, `cwc = list(ID, "M")`
#' centers `M` only in the latter role. When both GMC and CWC target the same
#' variable, GMC is applied first. Between-cluster means must be created by the
#' user and included explicitly in the formulas.
#'
#' \strong{Uncertainty Estimation}:
#' When using `$optimize(ci_method = "sampling")`, the function provides asymmetric
#' confidence intervals for indirect effects based on the distribution of products,
#' which is more accurate than the standard Sobel test (Delta Method).
#'
#' @return An `RTMB_Model` object.
#' @example inst/examples/ex_mediation.R
#' @export
rtmb_mediation <- function(formula, data, family = "gaussian", prior = prior_flat(),
                           y_range = NULL, fixed = NULL, view = NULL,
                           WAIC = FALSE, gmc = NULL, centering = NULL,
                           cwc = NULL, ...) {

  .check_unused_dots(..., .fn = "rtmb_mediation()")
  cwc_expr <- if (base::missing(cwc)) quote(NULL) else substitute(cwc)

  if (!is.list(formula)) stop("formula must be a list of formulas (e.g., list(M ~ X, Y ~ X + M)).")
  n_eq <- length(formula)
  if (n_eq < 1L || !all(vapply(formula, inherits, logical(1), what = "formula"))) {
    stop("'formula' must be a non-empty list of formula objects.", call. = FALSE)
  }

  random_bars <- lapply(formula, findbars)
  n_random_terms <- lengths(random_bars)
  if (any(n_random_terms > 1L)) {
    bad <- which(n_random_terms > 1L)
    stop(
      "Each mediation equation may contain at most one random-effect term. ",
      "Equation(s) ", paste(bad, collapse = ", "), " contain more than one.",
      call. = FALSE
    )
  }

  random_eq <- unname(which(n_random_terms == 1L))
  has_random <- length(random_eq) > 0L
  group_var <- NULL

  if (has_random) {
    group_vars <- character(length(random_eq))

    for (j in seq_along(random_eq)) {
      i <- random_eq[[j]]
      bar <- random_bars[[i]][[1L]]
      re_formula <- stats::as.formula(
        as.call(list(as.name("~"), bar[[2L]])),
        env = environment(formula[[i]])
      )
      re_terms <- stats::terms(re_formula)
      is_intercept_only <-
        identical(attr(re_terms, "intercept"), 1L) &&
        length(attr(re_terms, "term.labels")) == 0L

      if (!is_intercept_only) {
        stop(
          "Only random intercepts are currently supported in 'rtmb_mediation'. ",
          "Use '(1 | group)' in equation ", i, ".",
          call. = FALSE
        )
      }

      group_expr <- bar[[3L]]
      if (!is.symbol(group_expr)) {
        stop(
          "The grouping variable in equation ", i,
          " must be a single column name, such as '(1 | ID)'.",
          call. = FALSE
        )
      }
      group_vars[[j]] <- as.character(group_expr)
    }

    if (length(unique(group_vars)) > 1L) {
      stop(
        "Random intercepts across mediation equations must use the same grouping variable. ",
        "Found: ", paste(unique(group_vars), collapse = ", "), ".",
        call. = FALSE
      )
    }
    group_var <- group_vars[[1L]]
  }

  fixed_formulas <- lapply(formula, function(f) {
    fixed_f <- if (is.null(findbars(f))) f else nobars(f)
    if (!inherits(fixed_f, "formula")) {
      fixed_f <- stats::as.formula(fixed_f, env = environment(f))
    }
    fixed_f
  })

  # Validate: No duplicate response variables
  resp_check <- vapply(formula, function(f) as.character(f[[2]]), character(1))
  if (anyDuplicated(resp_check)) {
    dup_vars <- resp_check[duplicated(resp_check)]
    stop(
      sprintf("Duplicate response variable(s) detected: %s. Each equation must have a unique response.",
              paste(unique(dup_vars), collapse = ", ")),
      call. = FALSE
    )
  }

  if (is.null(prior)) prior <- prior_flat()

  # Automatically switch to prior_weak() if y_range is provided and prior is default flat
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && identical(prior$type, "flat")) {
    prior <- prior_weak()
  }

  prior <- .validate_prior_type(
    prior,
    allowed = c("flat", "normal", "weak"),
    context = "rtmb_mediation()"
  )

  # Prepare family list
  if (!is.list(family)) {
    family_list <- rep(list(family), n_eq)
  } else {
    if (length(family) != n_eq) stop("Length of family list must match length of formula list.")
    family_list <- family
  }
  valid_families <- c("gaussian", "bernoulli", "binomial", "poisson")
  bad_families <- unlist(family_list, use.names = FALSE)
  bad_families <- bad_families[
    !vapply(bad_families, function(x) is.character(x) && length(x) == 1L &&
              !is.na(x) && x %in% valid_families, logical(1))
  ]
  if (length(bad_families) > 0L) {
    stop(
      sprintf(
        "Invalid 'family' value in mediation model: %s. Valid options are: %s",
        paste(unique(bad_families), collapse = ", "),
        paste(valid_families, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(centering)) {
    if (!is.null(gmc) && !identical(gmc, centering)) {
      stop(
        "Specify only one of 'gmc' or 'centering', or use identical values.",
        call. = FALSE
      )
    }
    gmc <- centering
  }
  if (!is.null(gmc) &&
      (!is.character(gmc) || anyNA(gmc) || any(!nzchar(gmc)))) {
    stop("'gmc'/'centering' must be a character vector, 'all', or NULL.", call. = FALSE)
  }

  setup_input_data <- as.data.frame(data)
  cwc <- .normalize_glmer_cwc_spec(.resolve_glmer_cwc_spec(
    cwc_expr,
    parent.frame(),
    names(setup_input_data)
  ))
  cwc <- .match_glmer_cwc_cluster_column(cwc, setup_input_data)

  predictor_vars <- unique(unlist(lapply(fixed_formulas, function(f) {
    all.vars(f[[3L]], functions = FALSE)
  }), use.names = FALSE))
  available_predictors <- intersect(predictor_vars, names(setup_input_data))
  numeric_predictors <- available_predictors[
    vapply(setup_input_data[available_predictors], is.numeric, logical(1))
  ]

  validate_center_targets <- function(targets, option) {
    targets <- unique(as.character(targets))
    missing_targets <- setdiff(targets, names(setup_input_data))
    if (length(missing_targets) > 0L) {
      stop(
        "Variable(s) specified in '", option, "' were not found in data: ",
        paste(missing_targets, collapse = ", "), ".",
        call. = FALSE
      )
    }
    non_predictors <- setdiff(targets, predictor_vars)
    if (length(non_predictors) > 0L) {
      stop(
        "Variable(s) specified in '", option,
        "' are not predictors in any mediation equation: ",
        paste(non_predictors, collapse = ", "), ".",
        call. = FALSE
      )
    }
    non_numeric <- targets[!vapply(setup_input_data[targets], is.numeric, logical(1))]
    if (length(non_numeric) > 0L) {
      stop(
        "Only numeric predictors can be centered. Non-numeric variable(s) in '",
        option, "': ", paste(non_numeric, collapse = ", "), ".",
        call. = FALSE
      )
    }
    targets
  }

  target_gmc <- character(0)
  if (!is.null(gmc)) {
    target_gmc <- if (identical(gmc, "all")) {
      numeric_predictors
    } else {
      validate_center_targets(gmc, "gmc/centering")
    }
  }

  target_cwc <- character(0)
  cwc_cluster_setup_name <- NULL
  if (!is.null(cwc)) {
    cluster_var <- cwc$cluster
    if (is.character(cluster_var) && length(cluster_var) == 1L) {
      if (!(cluster_var %in% names(setup_input_data))) {
        stop(
          "Cluster variable '", cluster_var, "' for CWC was not found in data.",
          call. = FALSE
        )
      }
      cwc_cluster_setup_name <- cluster_var
    } else {
      if (length(cluster_var) != nrow(setup_input_data)) {
        stop(
          "A CWC cluster vector must have the same length as the model data (",
          nrow(setup_input_data), ").",
          call. = FALSE
        )
      }
      cwc_cluster_setup_name <- ".mediation_cwc_cluster"
      while (cwc_cluster_setup_name %in% names(setup_input_data)) {
        cwc_cluster_setup_name <- paste0(cwc_cluster_setup_name, "_")
      }
      setup_input_data[[cwc_cluster_setup_name]] <- cluster_var
    }

    target_cwc <- if (identical(cwc$pars, "all")) {
      setdiff(numeric_predictors, cwc_cluster_setup_name)
    } else {
      validate_center_targets(cwc$pars, "cwc")
    }
    cwc <- list(cluster = cwc_cluster_setup_name, pars = target_cwc)
  }

  setup_vars <- unique(c(
    unlist(lapply(formula, all.vars), use.names = FALSE),
    target_gmc,
    target_cwc,
    cwc_cluster_setup_name
  ))
  missing_vars <- setdiff(setup_vars, names(setup_input_data))
  if (length(missing_vars) > 0) {
    stop(
      "The following variables in formula are not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  setup_df <- stats::na.omit(setup_input_data[, setup_vars, drop = FALSE])
  class(setup_df) <- c("rtmb_setup_df", class(setup_df))

  predictor_df <- setup_df
  for (name in target_gmc) {
    predictor_df[[name]] <- center_grand_mean(predictor_df[[name]])
  }
  if (length(target_cwc) > 0L) {
    cluster_id <- setup_df[[cwc_cluster_setup_name]]
    for (name in target_cwc) {
      predictor_df[[name]] <- center_within_cluster(predictor_df[[name]], cluster_id)
    }
  }

  N <- nrow(setup_df)
  num_groups <- 0L
  group_levels <- character(0)
  if (has_random) {
    group_factor <- droplevels(as.factor(setup_df[[group_var]]))
    num_groups <- nlevels(group_factor)
    group_levels <- levels(group_factor)
    if (num_groups < 2L) {
      stop(
        "Random-intercept mediation requires at least two levels in grouping variable '",
        group_var, "'.",
        call. = FALSE
      )
    }
  }

  model_data <- setup_input_data
  resp_names <- character(n_eq)
  X_list <- list()
  X_colnames <- list()
  half_d_y_values <- vector("list", n_eq)
  mid_y_values <- vector("list", n_eq)
  has_predictor_centering <- length(target_gmc) > 0L || length(target_cwc) > 0L

  # 1. Parse Formulas and Prepare Data
  for (i in 1:n_eq) {
    f <- fixed_formulas[[i]]
    mf <- model.frame(f, data = setup_df)
    y_name <- as.character(formula[[i]][[2]])
    resp_names[i] <- y_name

    X_mat <- if (has_predictor_centering) {
      model.matrix(stats::delete.response(stats::terms(f)), data = predictor_df)
    } else {
      model.matrix(f, data = mf)
    }
    cols <- colnames(X_mat)
    cols[cols == "(Intercept)"] <- "Intercept"
    X_colnames[[i]] <- cols
    X_list[[i]] <- X_mat

    # Add range variables for weak priors dynamically
    if (inherits(prior, "rtmb_prior") && prior$type == "weak") {
       f_type <- family_list[[i]]
       if (!(f_type %in% c("bernoulli", "binomial", "poisson"))) {
         range_i <- if (is.list(y_range)) y_range[[y_name]] else y_range
         if (is.null(range_i)) {
           stop(paste0("y_range is required for response variable '", y_name, "' when using weakly informative priors. ",
                       "Please provide y_range as a vector or a named list (e.g., y_range = list(", y_name, " = c(1, 5)))."))
         }
         half_d_y_values[[i]] <- diff(range_i) / 2
         mid_y_values[[i]] <- mean(range_i)
       }
    }
  }

  prior_type <- prior$type
  if (prior_type == "normal") {
    if (is.null(prior$mu_sd) && !is.null(prior$Intercept_sd)) {
      prior$mu_sd <- prior$Intercept_sd
    }
  }
  if (prior_type == "weak") {
    if (is.null(prior$max_beta)) prior$max_beta <- 1.0
    if (is.null(prior$sd_ratio)) prior$sd_ratio <- 0.5
  }

  # 2. Setup AST Block
  setup_exprs <- list(
    bquote(df <- stats::na.omit(as.data.frame(.data)[, .(setup_vars), drop = FALSE])),
    quote(N <- nrow(df))
  )

  if (has_predictor_centering) {
    setup_exprs[[length(setup_exprs) + 1L]] <- quote(predictor_df <- df)
  }
  if (length(target_gmc) > 0L) {
    setup_exprs[[length(setup_exprs) + 1L]] <- "# Grand-mean centering of predictors"
    for (name in target_gmc) {
      target <- bquote(predictor_df[[.(name)]])
      center_call <- as.call(list(as.name("center_grand_mean"), target))
      setup_exprs[[length(setup_exprs) + 1L]] <- as.call(list(
        as.name("<-"), target, center_call
      ))
    }
  }
  if (length(target_cwc) > 0L) {
    setup_exprs[[length(setup_exprs) + 1L]] <- "# Centering predictors within cluster"
    cluster_expr <- bquote(df[[.(cwc_cluster_setup_name)]])
    for (name in target_cwc) {
      target <- bquote(predictor_df[[.(name)]])
      center_call <- as.call(list(
        as.name("center_within_cluster"), target, cluster_expr
      ))
      setup_exprs[[length(setup_exprs) + 1L]] <- as.call(list(
        as.name("<-"), target, center_call
      ))
    }
  }

  if (has_random) {
    setup_exprs[[length(setup_exprs) + 1L]] <-
      bquote(mediation_group <- droplevels(as.factor(df[[.(group_var)]])))
    setup_exprs[[length(setup_exprs) + 1L]] <-
      quote(mediation_group_idx <- as.integer(mediation_group))
    setup_exprs[[length(setup_exprs) + 1L]] <-
      quote(mediation_num_groups <- nlevels(mediation_group))
  }

  for (i in 1:n_eq) {
    mf_name <- as.name(paste0("mf_", i))
    Y_name <- as.name(paste0("Y_", i))
    X_name <- as.name(paste0("X_", i))
    formula_i <- fixed_formulas[[i]]

    setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(mf_name) <- model.frame(.(formula_i), df))
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(Y_name) <- as.numeric(model.response(.(mf_name))))
    if (has_predictor_centering) {
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(
        .(X_name) <- model.matrix(
          stats::delete.response(stats::terms(.(formula_i))),
          predictor_df
        )
      )
    } else {
      setup_exprs[[length(setup_exprs) + 1]] <-
        bquote(.(X_name) <- model.matrix(.(formula_i), .(mf_name)))
    }
  }

  for (i in 1:n_eq) {
    f_type <- family_list[[i]]
    p_name <- paste0("b", i)
    X_name <- as.name(paste0("X_", i))

    if (prior_type == "weak" && !(f_type %in% c("bernoulli", "binomial", "poisson"))) {
      half_d_y_name <- as.name(paste0("half_d_y_", i))
      mid_y_name <- as.name(paste0("mid_y_", i))
      setup_exprs[[length(setup_exprs) + 1]] <-
        bquote(.(half_d_y_name) <- .(half_d_y_values[[i]]))
      setup_exprs[[length(setup_exprs) + 1]] <-
        bquote(.(mid_y_name) <- .(mid_y_values[[i]]))
    }

    has_b_prior <- prior_type %in% c("weak", "normal") || !is.null(prior$b_sd) || !is.null(prior$Intercept_sd)
    has_sigma_prior <- f_type == "gaussian" && (prior_type %in% c("weak", "normal") || !is.null(prior$sigma_rate))

    if (has_b_prior) {
      X_sd_name <- as.name(paste0("X_sd_", i))
      X_mean_name <- as.name(paste0("X_mean_", i))
      b_prior_sd_name <- as.name(paste0(p_name, "_prior_sd"))
      b_prior_mean_name <- as.name(paste0(p_name, "_prior_mean"))
      intercept_prior_sd_name <- as.name(paste0("intercept_prior_sd_", i))

      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_sd_name) <- apply(.(X_name), 2, sd))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_sd_name)[.(X_sd_name) == 0] <- 1)
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_mean_name) <- rep(0, ncol(.(X_name))))

      if (prior_type == "weak") {
        X_c_name <- as.name(paste0("X_c_", i))
        mid_y_name <- as.name(paste0("mid_y_", i))

        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_mean_name) <- apply(.(X_name), 2, mean))
        has_intercept <- "Intercept" %in% X_colnames[[i]]
        if (has_intercept) {
           idx <- which(X_colnames[[i]] == "Intercept")
           setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_mean_name)[.(idx)] <- 0)
        }
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_c_name) <- .(X_name) - rep(1, N) %*% t(.(X_mean_name)))

        base_scale_name <- as.name(paste0("base_scale_", i))
        if (f_type %in% c("bernoulli", "binomial")) {
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(base_scale_name) <- pi / sqrt(3))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(intercept_prior_sd_name) <- .(base_scale_name) * .(prior$max_beta))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name) <- (.(prior$max_beta) * .(base_scale_name)) / .(X_sd_name))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(mid_y_name) <- 0)
        } else if (f_type == "poisson") {
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(base_scale_name) <- 1.0)
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(intercept_prior_sd_name) <- 1.0)
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name) <- 1.0 / .(X_sd_name))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(mid_y_name) <- 0)
        } else {
          half_d_y_name <- as.name(paste0("half_d_y_", i))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(base_scale_name) <- .(half_d_y_name) * .(prior$sd_ratio))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(intercept_prior_sd_name) <- .(half_d_y_name))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name) <- (.(prior$max_beta) * .(base_scale_name)) / .(X_sd_name))
          if (has_sigma_prior) {
            sigma_rate_name <- as.name(paste0("sigma", i, "_rate"))
            setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(sigma_rate_name) <- 1.0 / .(base_scale_name))
          }
        }
        if (has_intercept) {
           idx <- which(X_colnames[[i]] == "Intercept")
           setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_mean_name)[.(idx)] <- .(mid_y_name))
        }
      } else if (prior_type == "normal") {
        # Manual prior with explicit SDs
        b_sd_val <- if (!is.null(prior$b_sd)) prior$b_sd else 10
        int_sd_val <- if (!is.null(prior$mu_sd)) prior$mu_sd else 10
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name) <- rep(.(b_sd_val), ncol(.(X_name))))
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(intercept_prior_sd_name) <- .(int_sd_val))
        
        X_mean_name <- as.name(paste0("X_mean_", i))
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_mean_name) <- apply(.(X_name), 2, mean))
        has_intercept <- "Intercept" %in% X_colnames[[i]]
        if (has_intercept) {
           idx <- which(X_colnames[[i]] == "Intercept")
           setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_mean_name)[.(idx)] <- 0)
        }
        X_c_name <- as.name(paste0("X_c_", i))
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_c_name) <- .(X_name) - rep(1, N) %*% t(.(X_mean_name)))
      } else {
        # flat: No prior added by setup_exprs (handled by prior_exprs check later if any)
        # However, in current logic, if has_b_prior is FALSE, we don't enter here.
        # If it's flat, has_b_prior will be FALSE.
        NULL
      }
      # Overwrite intercept SD if needed
      has_intercept <- "Intercept" %in% X_colnames[[i]]
      if (has_intercept) {
         idx <- which(X_colnames[[i]] == "Intercept")
         setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name)[.(idx)] <- .(intercept_prior_sd_name))
      }
    }

    if (has_sigma_prior && prior_type == "normal") {
       sigma_rate_val <- if (!is.null(prior$sigma_rate)) prior$sigma_rate else 1
       sigma_rate_name <- as.name(paste0("sigma", i, "_rate"))
       setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(sigma_rate_name) <- .(sigma_rate_val))
    }
  }

  has_sd_re_prior <- FALSE
  if (has_random && prior_type != "flat") {
    if (!is.null(prior$tau_rate)) {
      setup_exprs[[length(setup_exprs) + 1L]] <-
        bquote(sd_re_rate <- .(prior$tau_rate))
      has_sd_re_prior <- TRUE
    } else if (prior_type == "weak") {
      rate_exprs <- lapply(random_eq, function(i) {
        bquote(1 / .(as.name(paste0("base_scale_", i))))
      })
      rate_vector <- as.call(c(list(as.name("c")), rate_exprs))
      setup_exprs[[length(setup_exprs) + 1L]] <- as.call(list(
        as.name("<-"), as.name("sd_re_rate"), rate_vector
      ))
      has_sd_re_prior <- TRUE
    }
  }

  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  # 3. Parameters, Transform and Model Block AST
  param_exprs <- list()
  tran_exprs <- list()
  model_exprs <- list()
  generate_exprs <- list()

  effect_names <- character(0)

  v_names <- list()
  init_list <- list()
  b_vars <- c()
  s_vars <- c()

  n_random_eq <- as.numeric(length(random_eq))
  random_labels <- character(0)
  if (has_random) {
    random_labels <- paste0(resp_names[random_eq], ":Intercept|", group_var)
    param_exprs[[length(param_exprs) + 1L]] <-
      bquote(sd_re <- Dim(.(n_random_eq), lower = 0))

    if (n_random_eq == 1L) {
      param_exprs[[length(param_exprs) + 1L]] <-
        quote(r_re <- Dim(mediation_num_groups, random = TRUE))
      init_list$r_re <- rep(0, num_groups)
      v_names$r_re <- group_levels
    } else {
      param_exprs[[length(param_exprs) + 1L]] <-
        bquote(CF_corr_re <- Dim(c(.(n_random_eq), .(n_random_eq)), type = "CF_corr"))
      param_exprs[[length(param_exprs) + 1L]] <-
        bquote(r_re <- Dim(c(mediation_num_groups, .(n_random_eq)), random = TRUE))
      init_list$r_re <- matrix(0, nrow = num_groups, ncol = n_random_eq)
      v_names$r_re <- list(group_levels, random_labels)
      v_names$corr_re <- random_labels
    }

    init_list$sd_re <- rep(1, n_random_eq)
    v_names$sd_re <- random_labels
  }

  for (i in 1:n_eq) {
    y_name <- resp_names[i]
    p_name <- paste0("b", i)
    p_c_name <- paste0("b_c", i)
    s_name <- paste0("sigma", i)
    f_type <- family_list[[i]]
    P_dim <- ncol(X_list[[i]])

    has_intercept <- "Intercept" %in% X_colnames[[i]]
    is_centered <- (prior_type %in% c("weak", "normal") && has_intercept)
    target_p_name <- if (is_centered) p_c_name else p_name

    param_exprs[[length(param_exprs) + 1]] <- bquote(.(as.name(target_p_name)) <- Dim(.(P_dim)))

    if (is_centered) {
       v_names[[p_c_name]] <- X_colnames[[i]]
       v_names[[p_name]] <- X_colnames[[i]]
       init_list[[p_c_name]] <- rep(0, P_dim)

       idx <- which(X_colnames[[i]] == "Intercept")
       X_mean_name <- as.name(paste0("X_mean_", i))

       tran_exprs[[length(tran_exprs) + 1]] <- bquote(.(as.name(p_name)) <- .(as.name(p_c_name)))
       tran_exprs[[length(tran_exprs) + 1]] <- bquote(.(as.name(p_name))[.(idx)] <- .(as.name(p_c_name))[.(idx)] - sum(.(X_mean_name) * .(as.name(p_c_name))))

       lin_pred_expr <- bquote(.(as.name(paste0("X_c_", i))) %*% .(as.name(p_c_name)))
    } else {
       v_names[[p_name]] <- X_colnames[[i]]
       init_list[[p_name]] <- rep(0, P_dim)
       lin_pred_expr <- bquote(.(as.name(paste0("X_", i))) %*% .(as.name(p_name)))
    }

    random_pos <- match(i, random_eq)
    if (!is.na(random_pos)) {
      if (n_random_eq == 1L) {
        lin_pred_expr <- bquote(
          .(lin_pred_expr) + sd_re[1] * r_re[mediation_group_idx]
        )
      } else {
        random_pos <- as.numeric(random_pos)
        lin_pred_expr <- bquote(
          .(lin_pred_expr) +
            sd_re[.(random_pos)] * r_re[mediation_group_idx, .(random_pos)]
        )
      }
    }

    if (f_type == "gaussian") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(.(as.name(s_name)) <- Dim(lower = 0))
      init_list[[s_name]] <- 1.0
      model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(paste0("Y_", i))) ~ normal(.(lin_pred_expr), .(as.name(s_name))))
      if (isTRUE(WAIC)) {
        generate_exprs[[length(generate_exprs) + 1]] <- as.call(list(
          as.name("<-"),
          as.name(paste0("log_lik_", i)),
          bquote(normal_lpdf(.(as.name(paste0("Y_", i))), .(lin_pred_expr), .(as.name(s_name)), sum = FALSE))
        ))
      }
    } else if (f_type %in% c("binomial", "bernoulli")) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(paste0("Y_", i))) ~ bernoulli_logit(.(lin_pred_expr)))
      if (isTRUE(WAIC)) {
        generate_exprs[[length(generate_exprs) + 1]] <- as.call(list(
          as.name("<-"),
          as.name(paste0("log_lik_", i)),
          bquote(bernoulli_logit_lpmf(.(as.name(paste0("Y_", i))), .(lin_pred_expr), sum = FALSE))
        ))
      }
    } else if (f_type == "poisson") {
      model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(paste0("Y_", i))) ~ poisson_log(.(lin_pred_expr)))
      if (isTRUE(WAIC)) {
        generate_exprs[[length(generate_exprs) + 1]] <- as.call(list(
          as.name("<-"),
          as.name(paste0("log_lik_", i)),
          bquote(poisson_lpmf(.(as.name(paste0("Y_", i))), exp(.(lin_pred_expr)), sum = FALSE))
        ))
      }
    }

    has_b_prior <- prior_type %in% c("weak", "normal")
    if (has_b_prior) {
      b_prior_sd_name <- as.name(paste0(p_name, "_prior_sd"))
      b_prior_mean_name <- as.name(paste0(p_name, "_prior_mean"))
      model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(target_p_name)) ~ normal(.(b_prior_mean_name), .(b_prior_sd_name)))
    }

    if (f_type == "gaussian") {
      has_sigma_prior <- prior_type %in% c("weak", "normal")
      if (has_sigma_prior) {
        rate_name <- as.name(paste0(s_name, "_rate"))
        model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(s_name)) ~ exponential(.(rate_name)))
      }
    }

    b_vars <- c(b_vars, p_name)
    if (f_type == "gaussian") s_vars <- c(s_vars, s_name)
  }

  if (has_random) {
    if (n_random_eq == 1L) {
      model_exprs[[length(model_exprs) + 1L]] <- quote(r_re ~ normal(0, 1))
    } else {
      tran_exprs[[length(tran_exprs) + 1L]] <-
        quote(corr_re <- CF_corr_re %*% t(CF_corr_re))

      lkj_eta <- if (prior_type == "flat") NULL else prior$lkj_eta %||% 1
      if (!is.null(lkj_eta)) {
        model_exprs[[length(model_exprs) + 1L]] <-
          bquote(CF_corr_re ~ lkj_CF_corr(.(lkj_eta)))
      }
      model_exprs[[length(model_exprs) + 1L]] <- bquote(
        for (g in 1:mediation_num_groups) {
          r_re[g, ] ~ multi_normal_CF(
            rep(0, .(n_random_eq)),
            rep(1, .(n_random_eq)),
            CF_corr_re
          )
        }
      )
    }

    if (has_sd_re_prior) {
      model_exprs[[length(model_exprs) + 1L]] <-
        quote(sd_re ~ exponential(sd_re_rate))
    }
  }

  # 4. Path Identification and Indirect Effects
  all_preds <- unique(unlist(X_colnames))
  all_resps <- unique(resp_names)
  mediators <- intersect(all_resps, all_preds)
  indeps <- setdiff(all_preds, c(all_resps, "Intercept"))
  if (length(indeps) == 0 && length(mediators) > 0) indeps <- setdiff(X_colnames[[1]], "Intercept")

  for (iv in indeps) {
    for (m in mediators) {
      idx_m_resp <- which(resp_names == m)
      idx_m_pred <- which(sapply(X_colnames, function(x) m %in% x))

      if (length(idx_m_resp) > 0 && length(idx_m_pred) > 0) {
        pos_iv <- which(X_colnames[[idx_m_resp]] == iv)
        if (length(pos_iv) > 0) {
          a_val <- bquote(.(as.name(paste0("b", idx_m_resp)))[.(pos_iv)])
          for (dv_idx in idx_m_pred) {
            if (dv_idx == idx_m_resp) next
            dv_name <- resp_names[dv_idx]
            pos_m <- which(X_colnames[[dv_idx]] == m)
            b_val <- bquote(.(as.name(paste0("b", dv_idx)))[.(pos_m)])
            ie_name <- paste0("IE_", iv, "_", m, "_", dv_name)
            tran_exprs[[ie_name]] <- bquote(.(as.name(ie_name)) <- .(a_val) * .(b_val))
            effect_names <- unique(c(effect_names, ie_name))
            pos_iv_direct <- which(X_colnames[[dv_idx]] == iv)
            if (length(pos_iv_direct) > 0) {
              de_name <- paste0("DE_", iv, "_", dv_name)
              de_val <- bquote(.(as.name(paste0("b", dv_idx)))[.(pos_iv_direct)])
              tran_exprs[[de_name]] <- bquote(.(as.name(de_name)) <- .(de_val))
              effect_names <- unique(c(effect_names, de_name))
              te_name <- paste0("TE_", iv, "_", m, "_", dv_name)
              tran_exprs[[te_name]] <- bquote(.(as.name(te_name)) <- .(as.name(ie_name)) + .(de_val))
              effect_names <- unique(c(effect_names, te_name))
            }
          }
        }
      }
    }
  }

  equation_df <- vapply(seq_len(n_eq), function(i) {
    if (family_list[[i]] == "gaussian") {
      as.numeric(N - qr(X_list[[i]])$rank)
    } else {
      Inf
    }
  }, numeric(1))
  names(equation_df) <- paste0("eq", seq_len(n_eq))

  # Fixed Gaussian equations use their own residual df. Derived effects inherit
  # df from the coefficients that contribute to their delta-method gradient.
  df_map <- NULL
  if (!has_random) {
    df_map <- list()
    for (i in seq_len(n_eq)) {
      df_map[[paste0("b", i)]] <- equation_df[[i]]
      if (family_list[[i]] == "gaussian") {
        df_map[[paste0("sigma", i)]] <- equation_df[[i]]
      }
    }
  }

  mdl_code <- list(
    setup = setup_ast,
    parameters = as.call(c(list(as.name("{")), param_exprs))
  )
  if (length(tran_exprs) > 0) {
    mdl_code$transform <- as.call(c(list(as.name("{")), tran_exprs))
  }
  mdl_code$model <- as.call(c(list(as.name("{")), model_exprs))

  if (length(generate_exprs) > 0) {
    if (isTRUE(WAIC)) {
      log_lik_names <- lapply(seq_len(n_eq), function(i) as.name(paste0("log_lik_", i)))
      joint_log_lik <- if (length(log_lik_names) == 1L) {
        log_lik_names[[1L]]
      } else {
        Reduce(function(a, b) as.call(list(as.name("+"), a, b)), log_lik_names)
      }
      generate_exprs[[length(generate_exprs) + 1]] <- joint_log_lik
      generate_exprs[[length(generate_exprs)]] <- as.call(list(as.name("<-"), as.name("log_lik"), generate_exprs[[length(generate_exprs)]]))
    }
    gen_ast <- as.call(c(list(as.name("{")), generate_exprs))
    mdl_code$generate <- if (isTRUE(WAIC)) .rtmb_waic_generate_ast(NULL, gen_ast) else gen_ast
  }
  mdl_code$env <- parent.frame()
  mdl_code$setup_env <- .rtmb_setup_env(environment(), setup_ast, exclude = names(model_data))

  random_view <- if (has_random) {
    c("sd_re", if (n_random_eq > 1L) "corr_re" else character(0))
  } else {
    character(0)
  }
  view_order <- c(b_vars, effect_names, s_vars, random_view)
  if (!is.null(view)) {
    view_order <- unique(c(view, view_order))
  }

  mdl <- rtmb_model(data = model_data, code = mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order, silent = FALSE)
  mdl$formula <- formula
  mdl$raw_data <- setup_df
  mdl$family <- family_list

  mdl$type <- "mediation"
  mdl$extra <- list(
    source = "wrapper",
    prior_type = prior$type,
    marginal = if (prior_type %in% c("weak", "normal")) paste0("b_c", 1:n_eq) else paste0("b", 1:n_eq),
    mediation = list(
      formula = formula,
      family = family_list,
      view = view,
      n_eq = n_eq,
      responses = resp_names,
      has_random = has_random,
      random_equations = random_eq,
      group = group_var,
      random_structure = if (has_random) "intercept" else NULL,
      gmc = target_gmc,
      centering = target_gmc,
      cwc = cwc
    )
  )

  if (has_random) {
    mdl$extra$mediation$num_groups <- num_groups
    mdl$extra$mediation$group_levels <- group_levels
    mdl$extra$mediation$random_labels <- random_labels
  }

  mdl$extra$df_map <- df_map
  mdl$extra$effect_names <- effect_names
  mdl$extra$equation_df <- equation_df

  return(mdl)
}
