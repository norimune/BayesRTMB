#' Fit a Correlation Model using RTMB
#'
#' @description
#' `rtmb_corr` fits a correlation model to estimate means, standard deviations, and correlation structures.
#' It supports simple correlation, multilevel correlation, and classical frequentist estimation.
#'
#' @param x A matrix, data frame, formula, or expression (e.g., \code{cbind(V1, V2)}) of response variables.
#' @param data An optional data frame containing the variables.
#' @param ID A character string or expression specifying the group ID variable for multilevel models.
#' @param covariates Optional numeric matrix or data frame of covariates to be included in the joint MVN model.
#' @param method Correlation method for \code{classic()}: \code{"pearson"}, \code{"spearman"}, or \code{"reml"}.
#' @param prior Prior configuration object: `prior_flat()`, `prior_normal()`, or `prior_weak()`. Default is `prior_flat()`.
#' @param y_range Optional numeric vector or matrix defining the theoretical range (min, max) of response variables.
#' Required when using \code{prior_weak()}. Can be a vector of length 2 (applies to all variables) or a matrix/list of length P.
#' @param init Optional list of initial values.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param missing Missing value handling strategy: "listwise" (default), "pairwise", or "fiml" (Full Information Maximum Likelihood).
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#'
#' @example inst/examples/ex_corr.R
#' @export
rtmb_corr <- function(x = NULL, data = NULL, ID = NULL,
                      covariates = NULL,
                       method = c("pearson", "spearman", "reml"),
                       prior = prior_flat(), y_range = NULL,
                       init = NULL, fixed = NULL,
                       missing = c("listwise", "fiml", "pairwise"),
                       WAIC = FALSE) {

  missing <- match.arg(missing)
  method <- match.arg(method)
  if (missing == "pairwise" && method == "reml") {
    stop("missing = 'pairwise' is only available when method is 'pearson' or 'spearman'.")
  }
  if (isTRUE(WAIC) && missing != "listwise") {
    stop("WAIC = TRUE is currently supported for rtmb_corr() only with missing = 'listwise'.", call. = FALSE)
  }
  x_expr <- substitute(x)
  id_expr <- substitute(ID)
  data_eval <- data
  if (!is.null(data_eval) && is.matrix(data_eval)) {
    data_eval <- as.data.frame(data_eval)
  }

  eval_col_range <- function(expr, data_env) {
    if (is.null(data_env) || !is.call(expr)) return(NULL)
    data_names <- names(data_env)
    if (is.null(data_names)) return(NULL)

    range_expr <- expr
    if (identical(expr[[1]], as.name("cbind")) && length(expr) == 2L) {
      range_expr <- expr[[2]]
    }
    if (!is.call(range_expr) || !identical(range_expr[[1]], as.name(":"))) {
      return(NULL)
    }

    left <- as.character(range_expr[[2]])
    right <- as.character(range_expr[[3]])
    if (length(left) != 1L || length(right) != 1L) return(NULL)
    if (!left %in% data_names || !right %in% data_names) return(NULL)

    idx <- match(c(left, right), data_names)
    data_env[, seq.int(idx[1L], idx[2L]), drop = FALSE]
  }

  # Evaluation logic for response variables (Y_mat)
  if (!is.null(data_eval)) {
    # Evaluate x in the context of data (NSE support for cbind(a, b))
    Y_mat <- eval_col_range(x_expr, data_eval)
    if (is.null(Y_mat)) {
      Y_mat <- eval(x_expr, data_eval, parent.frame())
    }
    if (is.null(Y_mat)) {
      Y_mat <- data_eval
    }
  } else {
    Y_mat <- x
  }

  # Handle if Y_mat is a formula
  if (inherits(Y_mat, "formula")) {
    formula <- Y_mat
    lhs_expr <- formula[[2]]

    # 1. Try to evaluate LHS directly (handles cbind(y1, y2) or Y_df)
    Y_mat <- try(eval(lhs_expr, data_eval, parent.frame()), silent = TRUE)

    # 2. If eval(lhs) failed, try a safer model.frame call
    if (inherits(Y_mat, "try-error") || is.null(Y_mat)) {
      resp_formula <- bquote(.(lhs_expr) ~ 1)
      mf_y <- try(model.frame(resp_formula, data = data_eval, na.action = na.pass), silent = TRUE)
      if (!inherits(mf_y, "try-error")) {
         Y_mat <- model.response(mf_y)
         if (is.null(Y_mat)) Y_mat <- mf_y[, 1, drop = FALSE]
      } else {
         # Last resort: check if it's a name in data
         id_name <- as.character(lhs_expr)
         if (!is.null(data_eval) && id_name %in% names(data_eval)) {
            Y_mat <- data_eval[[id_name]]
         }
      }
    }

    # Final check if it's a list (not df)
    if (is.list(Y_mat) && !is.data.frame(Y_mat)) {
       Y_mat <- do.call(cbind, Y_mat)
    }

    # If formula has RHS, extract covariates
    if (length(formula) == 3 && is.null(covariates)) {
       rhs_expr <- formula[[3]]
       cov_formula <- bquote(~ .(rhs_expr))
       mf_x <- try(model.frame(cov_formula, data = data_eval, na.action = na.pass), silent = TRUE)
       if (inherits(mf_x, "try-error")) {
          # Try to subset data to only RHS variables to avoid errors from other complex columns
          rhs_vars <- all.vars(rhs_expr)
          data_sub <- if (!is.null(data_eval)) data_eval[, intersect(rhs_vars, names(data_eval)), drop = FALSE] else data_eval
          mf_x <- model.frame(cov_formula, data = data_sub, na.action = na.pass)
       }
       covariates <- model.matrix(attr(mf_x, "terms"), mf_x)
       # Remove intercept
       if ("(Intercept)" %in% colnames(covariates)) {
         covariates <- covariates[, colnames(covariates) != "(Intercept)", drop = FALSE]
       }
    }
  } else {
    formula <- NULL
    # If Y_mat is a character vector of names, subset from data
    if (is.character(Y_mat) && length(Y_mat) > 1 && !is.null(data_eval)) {
      Y_mat <- data_eval[, Y_mat, drop = FALSE]
    }
    # Handle list of vectors
    if (is.list(Y_mat) && !is.data.frame(Y_mat)) {
       Y_mat <- do.call(cbind, Y_mat)
    }
  }

  # Handle covariates (if provided separately or extracted from formula)
  X_mat <- NULL
  if (!is.null(covariates)) {
     if (is.matrix(covariates) || is.data.frame(covariates)) {
        X_mat <- as.matrix(covariates)
     } else if (inherits(covariates, "formula")) {
        mf_x <- model.frame(covariates, data = data_eval, na.action = na.pass)
        X_mat <- model.matrix(covariates, mf_x)
        if ("(Intercept)" %in% colnames(X_mat)) {
          X_mat <- X_mat[, colnames(X_mat) != "(Intercept)", drop = FALSE]
        }
     } else if (is.character(covariates) && !is.null(data_eval)) {
        X_mat <- as.matrix(data_eval[, covariates, drop = FALSE])
     }
  }

  if (!is.null(X_mat)) {
    # Ensure numeric
    X_mat <- X_mat[, sapply(as.data.frame(X_mat), is.numeric), drop = FALSE]
  }

  # Parse ID (NSE support for ID = group)
  id_val <- NULL
  if (!is.null(id_expr)) {
    if (!is.null(data_eval)) {
      # Try to evaluate ID in data
      id_val <- try(eval(id_expr, data_eval, parent.frame()), silent = TRUE)
      if (inherits(id_val, "try-error") || is.null(id_val)) {
        # Fallback: check if ID name is in data
        id_name <- as.character(id_expr)
        if (id_name %in% names(data_eval)) {
          id_val <- data_eval[[id_name]]
        } else {
          id_val <- ID # Literal value
        }
      }
    } else if (is.data.frame(Y_mat)) {
      # If data is NULL but x is a dataframe, check if ID is a column name
      id_name <- as.character(id_expr)
      if (id_name %in% names(Y_mat)) {
        id_val <- Y_mat[[id_name]]
        # We will remove this column from Y_mat later
      } else {
        id_val <- ID
      }
    } else {
      id_val <- ID
    }

    # If the response matrix still contains the ID column, remove it
    id_name_str <- as.character(id_expr)
    if (is.data.frame(Y_mat) && id_name_str %in% colnames(Y_mat)) {
       Y_mat[[id_name_str]] <- NULL
    }
  }

  if (is.data.frame(Y_mat)) {
    if (!all(sapply(Y_mat, is.numeric))) {
      stop("All variables in the response data must be numeric. Character or factor variables are not supported.", call. = FALSE)
    }
  } else if (!is.numeric(Y_mat) && !is.logical(Y_mat)) {
    stop("The response data matrix must be numeric.", call. = FALSE)
  }

  if (missing == "listwise") {
    if (!is.null(X_mat)) {
      valid_idx <- stats::complete.cases(Y_mat, X_mat)
      Y_mat <- Y_mat[valid_idx, , drop = FALSE]
      X_mat <- X_mat[valid_idx, , drop = FALSE]
      if (length(id_val) == length(valid_idx)) id_val <- id_val[valid_idx]
    } else {
      valid_idx <- stats::complete.cases(Y_mat)
      Y_mat <- Y_mat[valid_idx, , drop = FALSE]
      if (length(id_val) == length(valid_idx)) id_val <- id_val[valid_idx]
    }
  } else if (missing == "pairwise" && anyNA(Y_mat)) {
     # No row deletion for pairwise
  }

  Y_mat <- as.matrix(Y_mat)

  P_y <- ncol(Y_mat)
  target_names <- colnames(Y_mat)
  if (is.null(target_names)) target_names <- paste0("Y", 1:P_y)

  # Combine with covariates for Joint MVN
  if (!is.null(X_mat)) {
     P_x <- ncol(X_mat)
     control_names <- colnames(X_mat)
     if (is.null(control_names)) control_names <- paste0("X", 1:P_x)
     Y_mat <- cbind(Y_mat, X_mat)
     colnames(Y_mat) <- c(target_names, control_names)
  } else {
     P_x <- 0
     control_names <- NULL
  }

  N <- nrow(Y_mat)
  P <- ncol(Y_mat)
  if (P < 1) stop("No numeric columns found for correlation analysis.")
  if (N < 2) {
    stop("Correlation analysis requires at least two complete observations.", call. = FALSE)
  }
  if (P == 1 && is.null(id_val)) {
    stop("Correlation analysis requires at least two variables.", call. = FALSE)
  }

  var_names <- colnames(Y_mat)
  if (is.null(var_names)) var_names <- paste0("V", 1:P)
  colnames(Y_mat) <- var_names



  # --- RTMB Implementation ---
  if (is.null(prior)) {
    prior <- prior_flat()
  }

  if (!is.null(id_val)) {
     # Multilevel correlation mode
     if (length(id_val) != N) {
       stop("'ID' must have the same length as the response data after missing-value handling.", call. = FALSE)
     }
     if (anyNA(id_val)) {
       stop("'ID' must not contain missing values after missing-value handling.", call. = FALSE)
     }
     group_factor <- as.factor(id_val)
     group_id <- as.integer(group_factor)
     group_names <- levels(group_factor)
     J <- length(group_names)
     group_counts <- tabulate(group_id, nbins = J)
     if (J < 2L) {
       stop("Multilevel correlation requires at least two ID groups.", call. = FALSE)
     }
     if (all(group_counts < 2L)) {
       stop(
         "Multilevel correlation requires at least one ID group with two or more observations. ",
         "All supplied ID groups contain only one observation, so within-group variation cannot be estimated.",
         call. = FALSE
       )
     }

     # Automatically switch to prior_weak() if y_range is provided and prior is default flat
     if (!is.null(y_range) && inherits(prior, "rtmb_prior") && identical(prior$type, "flat")) {
       prior <- prior_weak()
     }

     prior <- .validate_prior_type(
       prior,
       allowed = c("flat", "normal", "weak"),
       context = "rtmb_corr()"
     )

     prior_type <- prior$type

     is_flat <- identical(prior_type, "flat")
     is_normal <- identical(prior_type, "normal")
     is_weak <- identical(prior_type, "weak")

     if (is_normal) {
       if (is.null(prior$mu_sd) && !is.null(prior$Intercept_sd)) {
         prior$mu_sd <- prior$Intercept_sd
       }
     }

     use_weak_info <- prior_type %in% c("weak")
     multivariate <- P > 1

     setup_exprs <- list()
     if (use_weak_info) {
       if (is.null(y_range)) {
         stop("Specifying 'y_range' is required when using prior_weak().")
       }
       if (is.list(y_range)) {
         Y_range_mat <- do.call(rbind, lapply(y_range, as.numeric))
       } else if (is.vector(y_range) && length(y_range) == 2) {
         Y_range_mat <- matrix(rep(as.numeric(y_range), each = P), P, 2)
       } else {
         Y_range_mat <- as.matrix(y_range)
       }

       half_d_y <- (Y_range_mat[, 2] - Y_range_mat[, 1]) / 2
       mid_y_val <- (Y_range_mat[, 2] + Y_range_mat[, 1]) / 2
       base_scale <- half_d_y * prior$sd_ratio

       setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mid_y)
       setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- alpha_prior_sd)
       setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate_vec <- sigma_rate_vec)
     }
     setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

     param_exprs <- list(as.name("{"))
     param_exprs[[length(param_exprs) + 1]] <- bquote(mu <- Dim(.(P), random = TRUE))
     param_exprs[[length(param_exprs) + 1]] <- bquote(u <- Dim(c(J, .(P)), random = TRUE))

     param_exprs[[length(param_exprs) + 1]] <- bquote(sigma_between <- Dim(.(P), lower = 0))
     param_exprs[[length(param_exprs) + 1]] <- bquote(sigma_within <- Dim(.(P), lower = 0))

     if (multivariate) {
       param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr_between <- Dim(c(.(P), .(P)), type = "CF_corr"))
       param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr_within <- Dim(c(.(P), .(P)), type = "CF_corr"))
     }
     param_ast <- as.call(param_exprs)

     model_exprs <- list(as.name("{"))
      if (multivariate) {
        if (!is.null(prior$lkj_eta)) {
          model_exprs[[length(model_exprs) + 1]] <- bquote(L_corr_between ~ lkj_CF_corr(.(prior$lkj_eta)))
          model_exprs[[length(model_exprs) + 1]] <- bquote(L_corr_within ~ lkj_CF_corr(.(prior$lkj_eta)))
        }
        model_exprs[[length(model_exprs) + 1]] <- quote(u ~ multi_normal_CF(mean = rep(0, P), sd = sigma_between, CF_Omega = L_corr_between))
      } else {
       model_exprs[[length(model_exprs) + 1]] <- quote(u ~ normal(0, sigma_between))
     }

     model_exprs[[length(model_exprs) + 1]] <- quote(Y_pred <- u[group_id, ] + rep(mu, each = N))
     if (multivariate) {
       model_exprs[[length(model_exprs) + 1]] <- quote(Y ~ multi_normal_CF(mean = Y_pred, sd = sigma_within, CF_Omega = L_corr_within))
     } else {
       model_exprs[[length(model_exprs) + 1]] <- quote(Y ~ normal(Y_pred, sigma_within))
     }

     if (prior_type == "weak") {
       model_exprs[[length(model_exprs) + 1]] <- quote(sigma_between ~ exponential(sigma_rate_vec))
       model_exprs[[length(model_exprs) + 1]] <- quote(sigma_within ~ exponential(sigma_rate_vec))
       model_exprs[[length(model_exprs) + 1]] <- quote(mu ~ normal(mid_y, alpha_prior_sd))
      } else if (prior_type == "normal") {
        if (!is.null(prior$sigma_rate)) {
          model_exprs[[length(model_exprs) + 1]] <- bquote(sigma_between ~ exponential(.(prior$sigma_rate)))
          model_exprs[[length(model_exprs) + 1]] <- bquote(sigma_within ~ exponential(.(prior$sigma_rate)))
        }
        if (!is.null(prior$mu_sd) || !is.null(prior$Intercept_sd)) {
          mu_sd_val <- if (!is.null(prior$mu_sd)) prior$mu_sd else prior$Intercept_sd
          model_exprs[[length(model_exprs) + 1]] <- bquote(mu ~ normal(0, .(mu_sd_val)))
        }
      }
     model_ast <- as.call(model_exprs)

      transform_exprs <- list(as.name("{"))
      transform_exprs[[length(transform_exprs) + 1]] <- quote(ICC <- (sigma_between^2) / (sigma_between^2 + sigma_within^2))
      if (multivariate) {
        transform_exprs[[length(transform_exprs) + 1]] <- quote(B_corr <- L_corr_between %*% t(L_corr_between))
        transform_exprs[[length(transform_exprs) + 1]] <- quote(W_corr <- L_corr_within %*% t(L_corr_within))

        if (P_x > 0) {
           transform_exprs[[length(transform_exprs) + 1]] <- quote(R_yy <- W_corr[1:P_y, 1:P_y])
           transform_exprs[[length(transform_exprs) + 1]] <- quote(R_yx <- W_corr[1:P_y, (P_y+1):P])
           transform_exprs[[length(transform_exprs) + 1]] <- quote(R_xx <- W_corr[(P_y+1):P, (P_y+1):P])
           transform_exprs[[length(transform_exprs) + 1]] <- quote(P_cov_w <- R_yy - R_yx %*% solve(R_xx) %*% t(R_yx))
           transform_exprs[[length(transform_exprs) + 1]] <- quote(D_w <- diag(1 / sqrt(diag(P_cov_w))))
           transform_exprs[[length(transform_exprs) + 1]] <- quote(W_pcorr <- D_w %*% P_cov_w %*% D_w)

           transform_exprs[[length(transform_exprs) + 1]] <- quote(R_yy_b <- B_corr[1:P_y, 1:P_y])
           transform_exprs[[length(transform_exprs) + 1]] <- quote(R_yx_b <- B_corr[1:P_y, (P_y+1):P])
           transform_exprs[[length(transform_exprs) + 1]] <- quote(R_xx_b <- B_corr[(P_y+1):P, (P_y+1):P])
           transform_exprs[[length(transform_exprs) + 1]] <- quote(P_cov_b <- R_yy_b - R_yx_b %*% solve(R_xx_b) %*% t(R_yx_b))
           transform_exprs[[length(transform_exprs) + 1]] <- quote(D_b <- diag(1 / sqrt(diag(P_cov_b))))
           transform_exprs[[length(transform_exprs) + 1]] <- quote(B_pcorr <- D_b %*% P_cov_b %*% D_b)
           transform_exprs[[length(transform_exprs) + 1]] <- quote(report(ICC, B_corr, W_corr, W_pcorr, B_pcorr))
        } else {
           transform_exprs[[length(transform_exprs) + 1]] <- quote(report(ICC, B_corr, W_corr))
        }
      } else {
        transform_exprs[[length(transform_exprs) + 1]] <- quote(report(ICC))
      }
      transform_ast <- as.call(transform_exprs)

     generate_ast <- NULL
     if (isTRUE(WAIC)) {
       generate_ast <- if (multivariate) {
         quote({
           Y_pred <- u[group_id, ] + rep(mu, each = N)
           log_lik <- multi_normal_CF_lpdf(Y, Y_pred, sigma_within, L_corr_within, sum = FALSE)
         })
       } else {
         quote({
           Y_pred <- u[group_id, ] + rep(mu, each = N)
           log_lik <- normal_lpdf(Y, Y_pred, sigma_within, sum = FALSE)
         })
       }
     }

     data_list <- list(Y = Y_mat, group_id = group_id, N = N, P = P, J = J, P_y = P_y, P_x = P_x)
     if (use_weak_info) {
       data_list$mid_y <- mid_y_val
       data_list$alpha_prior_sd <- half_d_y
       data_list$sigma_rate_vec <- 1.0 / base_scale
     }

     mdl_code <- list(setup = setup_ast, parameters = param_ast, transform = transform_ast, model = model_ast, env = parent.frame())
     mdl_code$setup_env <- .rtmb_setup_env(environment(), setup_ast, exclude = names(data_list))
     if (!is.null(generate_ast)) mdl_code$generate <- .rtmb_waic_generate_ast(NULL, generate_ast)
     class(mdl_code) <- "rtmb_code"

     v_names <- list(mu = var_names, sigma_between = var_names, sigma_within = var_names)
     v_names$u <- list(group_names, var_names)
     v_names$ICC <- var_names

     if (multivariate) {
       v_names$B_corr <- list(var_names, var_names)
       v_names$W_corr <- list(var_names, var_names)
       if (P_x > 0) {
           v_names$B_pcorr <- list(target_names, target_names)
           v_names$W_pcorr <- list(target_names, target_names)
       }
     }

     init_list <- list(mu = colMeans(Y_mat))
     init_list$sigma_between <- apply(Y_mat, 2, sd) * 0.5
     init_list$sigma_within <- apply(Y_mat, 2, sd) * 0.5
     init_list$u <- matrix(0, J, P)

     view_order <- c("pcorr", "B_corr", "W_corr", "corr", "ICC", "mu", "sigma_between", "sigma_within", "sigma")

     obj <- rtmb_model(data_list, mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order)
     obj$raw_data <- data

     obj$type <- "corr"
     obj$extra <- list(
       source = "wrapper",
       prior_type = prior$type,
       marginal = "mu",
       corr_method = "reml",
       missing = missing
     )

     # Set degrees of freedom map for BW method
     # Between: J, Within: N - J - P
     df_map <- list()
     # Level-2 (Between)
     df_between <- pmax(J, 1)
     df_map$mu <- df_between
     df_map$sigma_between <- df_between
     df_map$B_corr <- df_between
     df_map$B_pcorr <- df_between
     df_map$ICC <- df_between
     
     # Level-1 (Within)
     df_within <- pmax(N - J - P, 1)
     df_map$sigma_within <- df_within
     df_map$W_corr <- df_within
     df_map$W_pcorr <- df_within
     
     obj$extra$df_map <- df_map

     return(obj)
  } else {
     # Simple correlation mode
     # Automatically switch to prior_weak() if y_range is provided and prior is default flat
     if (!is.null(y_range) && inherits(prior, "rtmb_prior") && identical(prior$type, "flat")) {
       prior <- prior_weak()
     }

     prior <- .validate_prior_type(
       prior,
       allowed = c("flat", "normal", "weak"),
       context = "rtmb_corr()"
     )

     prior_type <- prior$type

     is_flat <- identical(prior_type, "flat")
     is_normal <- identical(prior_type, "normal")
     is_weak <- identical(prior_type, "weak")

     if (is_normal) {
       if (is.null(prior$mu_sd) && !is.null(prior$Intercept_sd)) {
         prior$mu_sd <- prior$Intercept_sd
       }
     }

     setup_exprs <- list(as.name("{"))
     setup_exprs[[length(setup_exprs) + 1]] <- quote(N <- nrow(Y))
     setup_exprs[[length(setup_exprs) + 1]] <- quote(P <- ncol(Y))
     if (missing == "listwise") {
       setup_exprs[[length(setup_exprs) + 1]] <- quote(Y_bar <- colMeans(Y))
       setup_exprs[[length(setup_exprs) + 1]] <- quote(S_Y <- cov(Y) * (N - 1))
     } else if (missing == "pairwise") {
       setup_exprs[[length(setup_exprs) + 1]] <- quote(Y_bar <- colMeans(Y, na.rm = TRUE))
       setup_exprs[[length(setup_exprs) + 1]] <- quote(S_Y <- cov(Y, use = "pairwise.complete.obs") * (N - 1))
     }

     if (prior_type == "weak") {
       if (is.null(y_range)) {
         stop("y_range is required when using prior_weak().")
       }
       if (is.list(y_range)) {
         y_range_mat <- do.call(rbind, y_range)
       } else if (is.vector(y_range) && length(y_range) == 2) {
         y_range_mat <- matrix(y_range, P, 2, byrow = TRUE)
       } else {
         y_range_mat <- as.matrix(y_range)
       }

       mid_y_val <- (y_range_mat[, 1] + y_range_mat[, 2]) / 2
       alpha_prior_sd_val <- (y_range_mat[, 2] - y_range_mat[, 1]) / 2
       sigma_rate_val <- 1 / (alpha_prior_sd_val * 0.5)

       setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mid_y)
       setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- alpha_prior_sd)
       setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate_vec <- sigma_rate_vec)
     }
     setup_ast <- as.call(setup_exprs)

     param_exprs <- list(as.name("{"))
     param_exprs[[length(param_exprs) + 1]] <- quote(mean <- Dim(P, random = TRUE))
     param_exprs[[length(param_exprs) + 1]] <- quote(sd   <- Dim(P, lower = 0))

     if (P == 2) {
       param_exprs[[length(param_exprs) + 1]] <- bquote(corr <- Dim(lower = -1, upper = 1))
     } else {
       param_exprs[[length(param_exprs) + 1]] <- quote(CF_corr <- Dim(c(P, P), type = "CF_corr"))
     }
     param_ast <- as.call(param_exprs)

     model_exprs <- list(as.name("{"))
     if (prior_type == "weak") {
       model_exprs[[length(model_exprs) + 1]] <- quote(mean ~ normal(mid_y, alpha_prior_sd))
       model_exprs[[length(model_exprs) + 1]] <- quote(sd ~ exponential(sigma_rate_vec))
      } else if (prior_type == "normal") {
        if (!is.null(prior$mu_sd) || !is.null(prior$Intercept_sd)) {
          mu_sd_val <- if (!is.null(prior$mu_sd)) prior$mu_sd else prior$Intercept_sd
          model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(0, .(mu_sd_val)))
        }
        if (!is.null(prior$sigma_rate)) {
          model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(.(prior$sigma_rate)))
        }
      }

      if (P == 2) {
        if (!is.null(prior$lkj_eta)) {
          model_exprs[[length(model_exprs) + 1]] <- bquote(corr ~ lkj_corr(.(prior$lkj_eta)))
        }
        model_exprs[[length(model_exprs) + 1]] <- quote(CF_corr <- matrix(corr * 0, nrow = 2, ncol = 2))
        model_exprs[[length(model_exprs) + 1]] <- quote(CF_corr[1, 1] <- 1)
        model_exprs[[length(model_exprs) + 1]] <- quote(CF_corr[2, 1] <- corr)
        model_exprs[[length(model_exprs) + 1]] <- quote(CF_corr[2, 2] <- sqrt(1 - corr^2))
      } else {
        if (!is.null(prior$lkj_eta)) {
          model_exprs[[length(model_exprs) + 1]] <- bquote(CF_corr ~ lkj_CF_corr(.(prior$lkj_eta)))
        }
      }

      if (missing %in% c("listwise", "pairwise")) {
        model_exprs[[length(model_exprs) + 1]] <- quote(S_Y ~ sufficient_multi_normal_CF(N, Y_bar, mean, sd, CF_corr))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- quote(Y ~ fiml_multi_normal_CF(mean, sd, CF_corr))
      }
     model_ast <- as.call(model_exprs)

     transform_exprs <- list(as.name("{"))
     if (P > 2) {
       transform_exprs[[length(transform_exprs) + 1]] <- quote(corr <- CF_corr %*% t(CF_corr))

       if (P_x > 0) {
          transform_exprs[[length(transform_exprs) + 1]] <- quote(R_yy <- corr[1:P_y, 1:P_y])
          transform_exprs[[length(transform_exprs) + 1]] <- quote(R_yx <- corr[1:P_y, (P_y+1):P])
          transform_exprs[[length(transform_exprs) + 1]] <- quote(R_xx <- corr[(P_y+1):P, (P_y+1):P])
          transform_exprs[[length(transform_exprs) + 1]] <- quote(P_cov <- R_yy - R_yx %*% solve(R_xx) %*% t(R_yx))
          transform_exprs[[length(transform_exprs) + 1]] <- quote(D <- diag(1 / sqrt(diag(P_cov))))
          transform_exprs[[length(transform_exprs) + 1]] <- quote(pcorr <- D %*% P_cov %*% D)
          transform_exprs[[length(transform_exprs) + 1]] <- quote(report(corr, pcorr))
       } else {
          transform_exprs[[length(transform_exprs) + 1]] <- quote(report(corr))
       }
     }
     transform_ast <- as.call(transform_exprs)

     generate_ast <- NULL
     if (isTRUE(WAIC)) {
       generate_ast <- if (P == 1) {
         quote({
           log_lik <- normal_lpdf(Y, mean, sd, sum = FALSE)
         })
       } else {
         quote({
           log_lik <- multi_normal_CF_lpdf(Y, mean, sd, CF_corr, sum = FALSE)
         })
       }
     }

     mdl_code <- list(setup = setup_ast, parameters = param_ast, transform = transform_ast, model = model_ast, env = parent.frame())
     if (!is.null(generate_ast)) mdl_code$generate <- .rtmb_waic_generate_ast(NULL, generate_ast)
     class(mdl_code) <- "rtmb_code"

     dat_list <- list(Y = Y_mat, P_y = P_y, P_x = P_x)
     if (prior_type == "weak") {
       dat_list$mid_y <- mid_y_val
       dat_list$alpha_prior_sd <- alpha_prior_sd_val
       dat_list$sigma_rate_vec <- sigma_rate_val
     }
     mdl_code$setup_env <- .rtmb_setup_env(environment(), setup_ast, exclude = names(dat_list))

     v_names <- list(mean = var_names, sd = var_names)
     if (P == 2) {
       v_names$corr <- "rho"
     } else {
       v_names$corr <- list(var_names, var_names)
       if (P_x > 0) v_names$pcorr <- list(target_names, target_names)
     }

     init_list <- if (is.null(init)) {
       list(mean = colMeans(Y_mat), sd = apply(Y_mat, 2, sd))
     } else init

     view_vars <- if (P_x > 0) c("pcorr", "B_corr", "W_corr", "corr", "mean", "sd") else c("B_corr", "W_corr", "corr", "mean", "sd")
     obj <- rtmb_model(data = dat_list, code = mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_vars)

    obj$raw_data <- data
    obj$type <- "corr"
    obj$extra <- list(
      source = "wrapper",
      prior_type = prior$type,
      marginal = "mean",
      corr_method = method,
      missing = missing
    )

    return(obj)
  }
}
