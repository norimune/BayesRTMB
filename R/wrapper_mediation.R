#' RTMB-based Mediation Analysis Wrapper
#'
#' @description
#' `rtmb_mediation` performs mediation analysis by simultaneously estimating multiple
#' GLM regression equations. It automatically identifies mediation paths and calculates
#' indirect, direct, and total effects.
#'
#' @param formula A list of formulas defining the regression paths (e.g., `list(M ~ X, Y ~ X + M)`).
#' @param data A data frame containing the variables.
#' @param family A single character string or a list of character strings specifying the error distribution
#'   for each equation (e.g., `family = list("gaussian", "binomial")`). Default is "gaussian".
#' @param prior An object of class "rtmb_prior" specifying the prior distribution.
#' @param y_range Theoretical minimum and maximum values of the response variable.
#' @param fixed A named list of parameter values to fix (optional).
#' @param view Character vector of parameter names to prioritize in summary.
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Reserved; unused arguments are rejected.
#'
#' @details
#' The function identifies mediation paths by looking for variables that are
#' responses in one equation and predictors in another. Indirect effects are
#' calculated as the product of coefficients along these paths (\eqn{a * b}).
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
                           WAIC = FALSE, ...) {

  .check_unused_dots(..., .fn = "rtmb_mediation()")

  if (!is.list(formula)) stop("formula must be a list of formulas (e.g., list(M ~ X, Y ~ X + M)).")
  n_eq <- length(formula)

  # Validate: No random effects allowed in this version
  for (i in 1:n_eq) {
    if (any(grepl("\\|", as.character(formula[[i]])))) {
      stop("Multilevel mediation (random effects) is not yet supported in 'rtmb_mediation'. Please use fixed-effects formulas.")
    }
  }

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

  setup_vars <- unique(unlist(lapply(formula, all.vars), use.names = FALSE))
  missing_vars <- setdiff(setup_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "The following variables in formula are not found in data: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  setup_df <- stats::na.omit(as.data.frame(data)[, setup_vars, drop = FALSE])
  class(setup_df) <- c("rtmb_setup_df", class(setup_df))

  N <- nrow(setup_df)
  data_list <- list(df = setup_df)
  resp_names <- character(n_eq)
  X_list <- list()
  X_colnames <- list()

  # 1. Parse Formulas and Prepare Data
  for (i in 1:n_eq) {
    f <- formula[[i]]
    mf <- model.frame(f, data = data)
    y_name <- as.character(f[[2]])
    resp_names[i] <- y_name

    X_mat <- model.matrix(f, data = data)
    data_list[[paste0("formula_", i)]] <- f
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
         data_list[[paste0("half_d_y_", i)]] <- diff(range_i) / 2
         data_list[[paste0("mid_y_", i)]] <- mean(range_i)
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
  setup_exprs <- list()
  setup_exprs[[1]] <- quote(N <- nrow(df))

  for (i in 1:n_eq) {
    mf_name <- as.name(paste0("mf_", i))
    formula_name <- as.name(paste0("formula_", i))
    Y_name <- as.name(paste0("Y_", i))
    X_name <- as.name(paste0("X_", i))

    setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(mf_name) <- model.frame(.(formula_name), df))
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(Y_name) <- as.numeric(model.response(.(mf_name))))
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_name) <- model.matrix(.(formula_name), .(mf_name)))
  }

  for (i in 1:n_eq) {
    f_type <- family_list[[i]]
    p_name <- paste0("b", i)
    X_name <- as.name(paste0("X_", i))

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

  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(data_list)
  eval(setup_ast, tmp_env)

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

  # Calculate specific degrees of freedom for each equation
  df_map <- list()
  for (i in 1:n_eq) {
    df_val <- N - ncol(X_list[[i]])
    df_map[[paste0("b", i)]] <- df_val
    df_map[[paste0("sigma", i)]] <- df_val
  }
  # For IE/DE/TE, use the DF of the outcome equation
  for (iv in indeps) {
    for (m in mediators) {
      idx_m_resp <- which(resp_names == m)
      idx_m_pred <- which(sapply(X_colnames, function(x) m %in% x))
      if (length(idx_m_resp) > 0 && length(idx_m_pred) > 0) {
        for (dv_idx in idx_m_pred) {
          if (dv_idx == idx_m_resp) next
          dv_name <- resp_names[dv_idx]
          df_dv <- N - ncol(X_list[[dv_idx]])
          df_map[[paste0("IE_", iv, "_", m, "_", dv_name)]] <- df_dv
          df_map[[paste0("DE_", iv, "_", dv_name)]] <- df_dv
          df_map[[paste0("TE_", iv, "_", m, "_", dv_name)]] <- df_dv
        }
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
  mdl_code$env <- tmp_env
  mdl_code$setup_env <- .rtmb_setup_env(environment(), setup_ast, exclude = names(data_list))

  view_order <- c(b_vars, effect_names, s_vars)
  if (!is.null(view)) {
    view_order <- unique(c(view, view_order))
  }

  mdl <- rtmb_model(data = data_list, code = mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order, silent = FALSE)
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
      responses = resp_names
    )
  )

  mdl$extra$df_map <- df_map
  mdl$extra$effect_names <- effect_names
  mdl$extra$equation_df <- setNames(
    vapply(seq_len(n_eq), function(i) N - ncol(X_list[[i]]), numeric(1)),
    paste0("eq", seq_len(n_eq))
  )

  return(mdl)
}
