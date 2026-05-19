#' MAP fit object
#'
#' An R6 class storing optimization results from maximum a posteriori
#' (MAP) estimation.
#'
#' @field model The `RTMB_Model` object used for estimation.
#' @field par_vec Parameter vector on the unconstrained scale (constrained values unlisted).
#' @field par Parameter list on the constrained scale.
#' @field par_unc Parameter vector on the unconstrained scale (raw unconstrained values).
#' @field ci_method Method used for CI estimation ("wald", "sampling", or "none").
#' @field objective RTMB objective function object.
#' @field log_ml Log marginal likelihood or related model criterion.
#' @field convergence Optimizer convergence code.
#' @field sd_rep Standard deviation report object.
#' @field df_fixed Summary table for fixed-effect parameters.
#' @field random_effects Random effect estimates.
#' @field df_transform Summary table for transformed parameter estimates.
#' @field df_generate Summary table for generated quantity estimates.
#' @field opt_history A vector of optimize objective history.
#' @field transform List of transformed parameters maintaining their original dimensions.
#' @field generate List of generated quantities maintaining their original dimensions.
#' @field se_samples List of simulated samples for standard error estimation.
#' @field laplace Logical; whether Laplace approximation was used.
#' @field vcov_unc Variance-covariance matrix of parameters in unconstrained space.
#' @field map List; the parameter mapping used.
#' @field marginal_vars Character vector of parameter names requested through `optimize(marginal = ...)`.
#' @field laplace_random_vars Character vector of all parameter names passed to `MakeADFun(random = ...)` during Laplace approximation.
#' @field idx_fix_active Numeric vector; mapping between active parameters and full unconstrained vector.
#' @field show_df Logical; whether to display degrees of freedom in the summary output.
#' @field view Character vector of parameter names to prioritize in summary.
#' @field fallback_needed Logical; whether Hessian/SE fallback was used during optimization.
#'
#' @export
MAP_Fit <- R6::R6Class(
  classname = "map_fit",
  inherit = RTMB_Fit_Base,

  public = list(
    # --- Fields ---
    model          = NULL,
    par_vec        = NULL,
    par            = NULL,
    objective      = NULL,
    log_ml         = NULL,
    convergence    = NULL,
    sd_rep         = NULL,
    df_fixed       = NULL,
    random_effects = NULL,
    df_transform   = NULL,
    df_generate    = NULL,
    opt_history    = NULL,
    transform      = NULL,
    generate       = NULL,
    se_samples     = NULL,
    par_unc        = NULL,
    ci_method      = NULL,
    laplace        = NULL,
    map            = NULL,
    vcov_unc       = NULL,
    marginal_vars  = NULL,
    laplace_random_vars = NULL,
    idx_fix_active = NULL,
    show_df        = TRUE,
    view           = NULL,
    fallback_needed = NULL,

    #' @description Get point estimate for a target parameter.
    #' @param target Target parameter name.
    #' @param ... Additional arguments, ignored for MAP fits.
    #' @return Matrix, array, vector, or scalar point estimate.
    get_point_estimate = function(target, ...) {
      if (!is.null(self$par[[target]])) return(self$par[[target]])
      if (!is.null(self$transform[[target]])) return(self$transform[[target]])
      if (!is.null(self$generate[[target]])) return(self$generate[[target]])
      stop("Parameter not found: ", target)
    },



    #' @description Create a new `MAP_Fit` object.
    #' @param model The `RTMB_Model` object used for estimation.
    #' @param par_vec Parameter vector on the unconstrained scale (constrained values unlisted).
    #' @param par Parameter list on the constrained scale.
    #' @param objective The objective function value at the optimum.
    #' @param log_ml Log marginal likelihood.
    #' @param convergence Optimizer convergence code.
    #' @param sd_rep The `sdreport` object from TMB.
    #' @param df_fixed Data frame of fixed effects estimates and CIs.
    #' @param random_effects Data frame of random effects estimates and CIs.
    #' @param df_transform Data frame of transformed parameters.
    #' @param df_generate Data frame of generated quantities.
    #' @param opt_history Data frame of optimization history.
    #' @param transform List of transformed parameters maintaining their original dimensions.
    #' @param generate List of generated quantities maintaining their original dimensions.
    #' @param se_samples List of simulated samples for standard error estimation.
    #' @param par_unc Parameter vector on the unconstrained scale (raw values).
    #' @param ci_method Method used for CI estimation ("wald", "sampling", or "none").
    #' @param laplace Logical; whether Laplace approximation was used.
    #' @param vcov_unc Variance-covariance matrix of parameters in unconstrained space.
    #' @param map List; the parameter mapping used.
    #' @param marginal_vars Character vector of parameter names requested through `optimize(marginal = ...)`.
    #' @param laplace_random_vars Character vector of all parameter names passed to `MakeADFun(random = ...)` during Laplace approximation.
    #' @param idx_fix_active Numeric vector; mapping between active parameters and full unconstrained vector.
    #' @param show_df Logical; whether to display degrees of freedom in the summary output. Default is TRUE.
    #' @param view Character vector of parameter names to prioritize in summary.
    #' @param fallback_needed Logical; whether Hessian/SE fallback was used during optimization.
    initialize = function(model, par_vec = NULL, par = NULL, objective = NULL, log_ml = NULL,
                          convergence = NULL, sd_rep = NULL, df_fixed = NULL, random_effects = NULL,
                          df_transform = NULL, df_generate = NULL, opt_history = NULL,
                          transform = NULL, generate = NULL, se_samples = NULL, par_unc = NULL,
                          ci_method = "wald", laplace = TRUE, map = NULL, vcov_unc = NULL,
                          marginal_vars = NULL, laplace_random_vars = NULL, idx_fix_active = NULL,
                          show_df = TRUE, view = NULL, fallback_needed = NULL) {
      self$model <- model
      self$par_vec <- par_vec
      self$par <- par
      self$objective <- objective
      self$log_ml <- log_ml
      self$convergence <- convergence
      self$sd_rep <- sd_rep
      self$df_fixed <- df_fixed
      self$random_effects <- random_effects
      self$df_transform <- df_transform
      self$df_generate <- df_generate
      self$opt_history <- opt_history
      self$transform <- transform
      self$generate <- generate
      self$se_samples <- se_samples
      self$par_unc <- par_unc
      self$ci_method <- ci_method
      self$vcov_unc <- vcov_unc
      self$laplace <- laplace
      self$map <- map
      self$marginal_vars <- marginal_vars
      self$laplace_random_vars <- laplace_random_vars
      self$idx_fix_active <- idx_fix_active
      self$show_df <- show_df
      self$view <- view
      self$fallback_needed <- fallback_needed

      if (!is.null(self$par_vec)) self$par_vec <- Re(self$par_vec)
      if (!is.null(self$par)) self$par <- lapply(self$par, Re)
      if (!is.null(self$random_effects) && !is.data.frame(self$random_effects)) {
        self$random_effects <- lapply(self$random_effects, Re)
      }
      if (!is.null(self$transform)) self$transform <- lapply(self$transform, Re)
      if (!is.null(self$generate)) self$generate <- lapply(self$generate, Re)
    },

    #' @description Return random effect estimates as a named list.
    #' @return A named list of random effect estimates.
    ranef = function() {
      if (is.null(self$random_effects)) {
        message("No random effects found in this model.")
        return(NULL)
      }
      
      # Extract only the parameters that are marked as random
      random_names <- names(self$model$par_list)[sapply(self$model$par_list, function(x) isTRUE(x$random))]
      return(self$par[random_names])
    },

    #' @description Extract samples from the asymptotic posterior distribution.
    #' @param pars Character or numeric vector specifying the names or indices of parameters to extract.
    #' @param inc_random Logical; whether to include random effects.
    #' @param inc_transform Logical; whether to include transformed parameters.
    #' @param inc_generate Logical; whether to include generated quantities.
    #' @param ... Ignored.
    #' @return An array of samples [iterations, 1, parameters].
    draws = function(pars = NULL, inc_random = FALSE, inc_transform = TRUE, inc_generate = TRUE, ...) {
      # 1. Determine number of samples
      num_samples <- 1
      if (!is.null(self$se_samples)) {
        num_samples <- nrow(self$se_samples$con[[1]])
      }

      # 2. Identify candidate base variables
      all_base_vars <- names(self$model$par_list)
      if (!inc_random) {
        random_flags <- sapply(self$model$par_list, function(x) isTRUE(x$random))
        all_base_vars <- all_base_vars[!random_flags]
      }
      
      if (inc_transform && !is.null(self$transform)) all_base_vars <- c(all_base_vars, names(self$transform))
      if (inc_generate && !is.null(self$generate)) all_base_vars <- c(all_base_vars, names(self$generate))

      # 3. Narrow down base variables if pars is provided as character
      target_base_vars <- all_base_vars
      if (!is.null(pars) && is.character(pars)) {
        pars_base <- unique(gsub("\\[.*\\]$", "", pars))
        # We need to include base variables that are directly in pars_base,
        # but also we don't know if a flat name matches until we expand.
        # To be safe and efficient, we take the union of all_base_vars and pars_base.
        target_base_vars <- intersect(all_base_vars, pars_base)
        if (length(target_base_vars) == 0) target_base_vars <- all_base_vars
      }

      # 4. Collect and flatten samples for target base variables
      flat_list <- list()
      flat_names_vec <- c()

      for (v in target_base_vars) {
        # Determine group
        group <- if (v %in% names(self$model$par_list)) "con" 
                 else if (v %in% names(self$transform)) "tran"
                 else if (v %in% names(self$generate)) "gq"
                 else NULL
        
        if (is.null(group)) next

        # Get samples matrix [num_samples, length(v)]
        samps <- if (!is.null(self$se_samples[[group]][[v]])) {
          self$se_samples[[group]][[v]]
        } else {
          val <- if (group == "con") self$par[[v]] else self[[if(group=="gq")"generate" else "transform"]][[v]]
          matrix(rep(as.numeric(val), each = num_samples), nrow = num_samples)
        }

        flat_list[[v]] <- samps
        
        # Get flat names
        v_info <- self$model$par_list[[v]]
        v_dim <- if (!is.null(v_info)) v_info$dim else dim(self[[if(group=="gq")"generate" else "transform"]][[v]])
        if (is.null(v_dim)) v_dim <- length(as.numeric(self[[if(group=="gq")"generate" else "transform"]][[v]]))
        v_names_def <- self$model$par_names[[v]]
        flat_names_vec <- c(flat_names_vec, generate_flat_names(v, v_dim, v_names_def))
      }

      if (length(flat_list) == 0) stop("No matching parameters found.")

      res_mat <- do.call(cbind, unname(flat_list))
      colnames(res_mat) <- flat_names_vec

      # 5. Final filtering by pars (handles both base names and flattened names)
      if (!is.null(pars)) {
        if (is.numeric(pars)) {
          target_idx <- pars[pars >= 1 & pars <= ncol(res_mat)]
        } else if (is.character(pars)) {
          base_names <- gsub("\\[.*\\]$", "", flat_names_vec)
          target_idx <- integer(0)
          for (p in pars) {
            match_idx <- which(flat_names_vec == p | base_names == p)
            target_idx <- c(target_idx, match_idx)
          }
          target_idx <- unique(target_idx)
        } else {
          stop("'pars' must be numeric or character.")
        }
        
        if (length(target_idx) == 0) stop("No matching parameters found in 'pars'.")
        res_mat <- res_mat[, target_idx, drop = FALSE]
        flat_names_vec <- flat_names_vec[target_idx]
      }

      # 6. Convert to array [iterations, chains=1, parameters]
      res_array <- array(res_mat, dim = c(num_samples, 1, ncol(res_mat)))
      dimnames(res_array) <- list(iteration = 1:num_samples, chain = 1, variable = flat_names_vec)

      return(res_array)
    },

    #' @description Summarize MAP estimates.
    #' @param pars Character vector specifying the names of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param max_rows Maximum number of rows to print in summaries. Default is 10.
    #' @param digits Number of digits to print.
    #' @param ranef Logical; whether to also display random effect estimates. Default is FALSE.
    #' @param view Character vector of parameter names to prioritize or filter by.
    #' @return A summary object, typically a data frame.
    summary = function(pars = NULL, max_rows = 10, digits = 5, ranef = FALSE, view = NULL) {
      # Defensive check for pars
      vars <- if (is.null(pars)) names(self$par) else pars
      if (is.numeric(vars) && !is.null(self$par)) vars <- names(self$par)[vars]
      vars <- as.character(vars) # Fix for environment subsetting error

      cat("\nCall:\nMAP Estimation via RTMB\n")
      cat(sprintf("\nNegative Log-Posterior: %.2f\n", self$objective))

      if (!is.null(self$log_ml) && !is.na(self$log_ml)) {
        cat(sprintf("Approx. Log Marginal Likelihood (Laplace): %.2f\n", self$log_ml))
      } else {
        cat("Approx. Log Marginal Likelihood (Laplace): NA\n")
      }

      if (!is.null(self$random_effects) && !isTRUE(ranef)) {
        cat("Note: Random effects are stored in $random_effects (use ranef = TRUE to show them)\n")
      }

      target_views <- c()
      if (!is.null(view)) {
        target_views <- c(target_views, view)
      } else if (!is.null(self$view)) {
        target_views <- c(target_views, self$view)
      } else if (!is.null(self$model$view)) {
        target_views <- c(target_views, self$model$view)
      }

      view_matches_ranef <- FALSE
      if (length(target_views) > 0 && is.data.frame(self$random_effects) && nrow(self$random_effects) > 0) {
        re_names <- rownames(self$random_effects)
        re_base_names <- gsub("\\[.*\\]$", "", re_names)
        for (v in target_views) {
          if (any(re_names == v | re_base_names == v)) {
            view_matches_ranef <- TRUE
            break
          }
        }
      }

      all_dfs <- list()
      if (!is.null(self$df_fixed) && nrow(self$df_fixed) > 0) all_dfs$fixed <- self$df_fixed
      
      # Automatically include random effects if requested by pars, ranef=TRUE, or view.
      include_ranef <- isTRUE(ranef) || !is.null(pars) || view_matches_ranef
      if (include_ranef && is.data.frame(self$random_effects) && nrow(self$random_effects) > 0) {
        all_dfs$random <- self$random_effects
      }
      if (!is.null(self$df_transform) && nrow(self$df_transform) > 0) all_dfs$transform <- self$df_transform
      if (!is.null(self$df_generate) && nrow(self$df_generate) > 0) all_dfs$generate <- self$df_generate

      if (length(all_dfs) == 0) {
        cat("\nNo parameters to display.\n")
        return(invisible(NULL))
      }

      df_combined <- do.call(rbind, unname(all_dfs))

      if (!is.null(pars)) {
        var_names <- rownames(df_combined)
        base_names <- gsub("\\[.*\\]$", "", var_names)
        
        ordered_idx <- integer(0)
        for (p in vars) {
          match_idx <- which(var_names == p | base_names == p)
          ordered_idx <- c(ordered_idx, match_idx)
        }
        df_combined <- df_combined[unique(ordered_idx), , drop = FALSE]
      }

      if (nrow(df_combined) == 0) {
        cat("\nNo matching parameters found.\n")
        return(invisible(df_combined))
      }

      if (nrow(df_combined) > 0 && is.null(pars)) {
        var_names <- rownames(df_combined)
        base_names <- gsub("\\[.*\\]$", "", var_names)

        priority_idx <- integer(0)
        for (v in target_views) {
          match_idx <- which(var_names == v | base_names == v)
          priority_idx <- c(priority_idx, match_idx)
        }
        priority_idx <- unique(priority_idx)
        other_idx <- setdiff(seq_along(var_names), priority_idx)
        df_combined <- df_combined[c(priority_idx, other_idx), , drop = FALSE]
      }

      ci_label <- "95% Wald CI"
      if (!is.null(self$ci_method)) {
        if (self$ci_method == "sampling") ci_label <- "95% Sampling-based CI"
      }
      cat(sprintf("\nPoint Estimates and %s:\n", ci_label))

      num_cols <- sapply(df_combined, is.numeric)
      df_combined[num_cols] <- lapply(df_combined[num_cols], function(x) {
        x[abs(x) < 1e-12 & !is.na(x)] <- 0
        return(x)
      })

      # Display all estimated quantities including SE and CI (excluding only metadata if needed)
      df_display <- df_combined
      
      # Ensure consistent column ordering for readability (support both DF and df)
      # Order: Estimate, Std. Error, [CIs], DF/df, [Other]
      desired_order <- c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "Lower 2.5%", "Upper 97.5%", "df")
      current_names <- colnames(df_display)
      order_to_use <- intersect(desired_order, current_names)
      other_names <- setdiff(current_names, desired_order)
      df_display <- df_display[, c(order_to_use, other_names), drop = FALSE]

      out_df <- data.frame(variable = rownames(df_display), df_display, check.names = FALSE, stringsAsFactors = FALSE)
      
      # --- Hide DF column if requested ---
      if (!isTRUE(self$show_df) && "df" %in% names(out_df)) {
         out_df$df <- NULL
      }
      
      rownames(out_df) <- NULL

      if (!is.null(max_rows) && nrow(out_df) > max_rows) {
        out_df <- head(out_df, max_rows)
      }

      class(out_df) <- c("summary_BayesRTMB", "data.frame")
      print(out_df, digits = digits)

      cat("\n")
      invisible(df_combined)
    },


    #' @description Print a brief summary of the fitted object.
    #' @param pars Character vector specifying the names of parameters to summarize.
    #' @param max_rows Maximum number of rows to print in summaries.
    #' @param digits Number of digits to print.
    #' @param ... Additional arguments passed to the `summary` method.
    #' @return The object itself, invisibly.
    print = function(pars = NULL, max_rows = 10, digits = 5, ...) {
      self$summary(pars = pars, max_rows = max_rows, digits = digits, ...)
      invisible(self)
    },

    #' @description Compute generated quantities from the MAP estimate.
    #' @param code An `rtmb_code(\{ ... \})` or `\{ ... \}` block containing the logic to be calculated using the MAP estimate.
    #' @return The `MAP_Fit` object itself (invisibly). Results are added or updated in the `generate` list and `df_generate`.
    generated_quantities = function(code) {
      raw_code <- substitute(code)
      if (is.name(raw_code)) {
        evaluated <- tryCatch(eval(raw_code, envir = parent.frame()), error = function(e) NULL)
        if (is.language(evaluated) || is.call(evaluated)) code <- evaluated
      }

      gen_ast <- if (is.call(code) && identical(code[[1]], as.name("rtmb_code"))) code$generate else code

      gen_fn <- eval(bquote(transform_code(.(gen_ast))))
      environment(gen_fn) <- parent.env(globalenv())

      p_list <- self$par
      if (!is.null(self$transform)) p_list <- c(p_list, self$transform)
      if (!is.null(self$generate)) p_list <- c(p_list, self$generate)

      res <- gen_fn(self$model$data, p_list)

      if (is.null(self$generate)) self$generate <- list()
      for (n in names(res)) {
        self$generate[[n]] <- res[[n]]
      }

      flat_base <- unlist(res, use.names = FALSE)
      names_vec <- c()
      for (name in names(res)) {
        val <- res[[name]]
        dim_val <- dim(val)
        if (is.null(dim_val)) dim_val <- length(val)
        names_def <- self$model$par_names[[name]]
        names_vec <- c(names_vec, generate_flat_names(name, dim_val, names_def))
      }

      se_out <- rep(NA, length(flat_base))
      low_out <- rep(NA, length(flat_base))
      up_out <- rep(NA, length(flat_base))

      # Apply se_sampling if available
      if (!is.null(self$se_samples)) {
        num_samples <- nrow(self$se_samples$con[[1]])
        samps_res <- list()

        for (s in 1:num_samples) {
          tmp_p_list <- list()
          for (name in names(self$se_samples$con)) {
            val <- self$se_samples$con[[name]][s, ]
            if (length(dim(self$par[[name]])) > 1) dim(val) <- dim(self$par[[name]])
            tmp_p_list[[name]] <- val
          }
          if (!is.null(self$se_samples$tran)) {
            for (name in names(self$se_samples$tran)) {
              val <- self$se_samples$tran[[name]][s, ]
              if (length(dim(self$transform[[name]])) > 1) dim(val) <- dim(self$transform[[name]])
              tmp_p_list[[name]] <- val
            }
          }
          if (!is.null(self$se_samples$gq)) {
            for (name in names(self$se_samples$gq)) {
              val <- self$se_samples$gq[[name]][s, ]
              if (length(dim(self$generate[[name]])) > 1) dim(val) <- dim(self$generate[[name]])
              tmp_p_list[[name]] <- val
            }
          }

          tmp_res <- tryCatch(gen_fn(self$model$data, tmp_p_list), error = function(e) NULL)
          if (!is.null(tmp_res)) {
            for (name in names(tmp_res)) {
              val <- as.numeric(tmp_res[[name]])
              if (is.null(samps_res[[name]])) samps_res[[name]] <- matrix(NA, num_samples, length(val))
              samps_res[[name]][s, ] <- val
            }
          }
        }

        if (length(samps_res) > 0) {
          if (is.null(self$se_samples$gq)) self$se_samples$gq <- list()
          for (name in names(samps_res)) {
            self$se_samples$gq[[name]] <- samps_res[[name]]
          }

          mat_list <- list()
          for (name in names(res)) {
            if (!is.null(samps_res[[name]])) {
              mat_list[[name]] <- samps_res[[name]]
            } else {
              mat_list[[name]] <- matrix(NA, num_samples, length(res[[name]]))
            }
          }
          mat_all <- do.call(cbind, unname(mat_list))

          se_out <- apply(mat_all, 2, sd, na.rm = TRUE)
          low_out <- apply(mat_all, 2, quantile, probs = 0.025, na.rm = TRUE)
          up_out <- apply(mat_all, 2, quantile, probs = 0.975, na.rm = TRUE)
        }
      }

      new_df <- data.frame(
        Estimate     = flat_base,
        `Std. Error` = se_out,
        `Lower 95%`  = low_out,
        `Upper 95%`  = up_out,
        row.names    = names_vec,
        check.names  = FALSE
      )

      if (is.null(self$df_generate)) {
        self$df_generate <- new_df
      } else {
        overlap <- intersect(rownames(self$df_generate), rownames(new_df))
        if (length(overlap) > 0) {
          self$df_generate[overlap, ] <- new_df[overlap, ]
          new_df <- new_df[!(rownames(new_df) %in% overlap), , drop = FALSE]
        }
        if (nrow(new_df) > 0) {
          self$df_generate <- rbind(self$df_generate, new_df)
        }
      }

      cat("Generated quantities updated.\n")
      invisible(self)
    },

    #' @description Compute an approximate WAIC from sampling-based uncertainty propagation.
    #' @return A `waic_BayesRTMB` object.
    WAIC = function() {
      if (is.null(self$se_samples) || is.null(self$se_samples$gq) ||
          is.null(self$se_samples$gq$log_lik)) {
        stop(
          "WAIC for MAP_Fit requires pointwise `log_lik` samples. ",
          "Use optimize(se_method = 'sampling') with WAIC = TRUE, or use sample()/variational().",
          call. = FALSE
        )
      }
      .compute_waic_from_log_lik(self$se_samples$gq$log_lik)
    },

    #' @description Run basic diagnostics for the MAP fit.
    #' @param ... Additional arguments passed to `diagnose_map_fit()`.
    #' @return A `diagnose_BayesRTMB` object.
    diagnose = function(...) {
      diagnose_map_fit(self, ...)
    },

    #' @description Calculate Profile Likelihood confidence intervals for specific parameters.
    #' @param pars Character vector of parameter names to profile. If NULL, all fixed parameters are profiled.
    #' @param level Confidence level (default is 0.95).
    #' @param trace Logical; whether to print profiling progress. Default is FALSE.
    #' @param digits Integer; number of decimal places to print. Default is 5.
    #' @param show_plot Logical; whether to plot the profile likelihood curves. Default is FALSE.
    #' @param quiet Logical; whether to suppress text output. Default is FALSE.
    #' @param jacobian Character; "none" (default), "random", or "all". Whether to include Jacobian adjustments for transformations.
    #' @param ... Additional arguments passed to TMB::tmbprofile (e.g., ytol).
    #' @return A data frame containing the profile-based confidence intervals, with the raw profile objects stored in the "profiles" attribute.
    profile = function(pars = NULL, level = 0.95, trace = FALSE, digits = 5, show_plot = FALSE, quiet = FALSE, jacobian = "none", ...) {
      if (!quiet) cat("Estimating confidence intervals via Profile Likelihood...\n")
      
      # 1. Re-build ad_obj using optimized values and specific jacobian adjustment
      # Use stored marginal_vars to ensure marginal likelihood is profiled if applicable
      if (jacobian == "none" && isTRUE(self$laplace)) jacobian <- "random"
      
      # Check if any target parameters are currently marginalized
      if (!is.null(pars) && !is.null(self$marginal_vars)) {
          target_base_names <- unique(gsub("\\[.*\\]$", "", pars))
          conflicting_vars <- intersect(target_base_names, self$marginal_vars)
          if (length(conflicting_vars) > 0) {
              stop(sprintf("Parameter(s) '%s' are marginalized (integrated out) and cannot be profiled. 
To profile these parameters, please re-run optimize() without including them in 'marginal'.", 
                           paste(conflicting_vars, collapse = ", ")), call. = FALSE)
          }
      }
      
      # Use stored laplace_random_vars to ensure all random components are integrated in the re-built object
      profile_random_vars <- self$laplace_random_vars %||% self$marginal_vars
      
      ad_setup <- self$model$build_ad_obj(init = self$par, laplace = self$laplace, 
                                          map = self$map, jacobian_target = jacobian,
                                          .marginal_vars = profile_random_vars)
      ad_obj <- ad_setup$ad_obj
      
      # Sync internal state
      ad_obj$fn(ad_obj$par)
      
      # 2. Identify indices
      all_fixed_names <- rownames(self$df_fixed)
      if (is.null(all_fixed_names)) {
        stop("No fixed parameters found to profile.")
      }
      
      if (is.null(pars)) {
        target_indices <- seq_along(all_fixed_names)
        target_names <- all_fixed_names
      } else {
        target_indices <- integer(0)
        target_names <- character(0)
        
        for (p in pars) {
          idx <- which(all_fixed_names == p)
          if (length(idx) == 0) {
            idx <- grep(paste0("^", p, "(\\[|$)"), all_fixed_names)
          }
          if (length(idx) > 0) {
            target_indices <- c(target_indices, idx)
            target_names <- c(target_names, all_fixed_names[idx])
          } else {
            warning(sprintf("Parameter '%s' not found in estimated fixed effects.", p), call. = FALSE)
          }
        }
        target_indices <- unique(target_indices)
        target_names <- unique(target_names)
      }
      
      if (length(target_indices) == 0) {
        stop("No valid parameters selected for profiling.", call. = FALSE)
      }
      
      res_mat <- matrix(NA, nrow = length(target_indices), ncol = 2)
      low_label <- sprintf("Lower %g%%", (1 - level)/2 * 100)
      up_label <- sprintf("Upper %g%%", (1 - (1 - level)/2) * 100)
      colnames(res_mat) <- c(low_label, up_label)
      rownames(res_mat) <- target_names

      # Use the stored idx_fix_active for reliable index mapping
      idx_fix_active <- self$idx_fix_active
      
      if (is.null(idx_fix_active)) {
        stop("idx_fix_active not found in MAP_Fit object. Profile calculation unavailable.")
      }

      prof_list <- list()
      
      # 3. Profile
      for (i in seq_along(target_indices)) {
        idx <- target_indices[i]
        p_name <- target_names[i]
        if (!quiet) {
          cat(sprintf("  Profiling %d/%d (%s)...                                                                \r", i, length(target_indices), p_name))
          utils::flush.console()
        }
        
        args_list <- list(...)
        if (!"ytol" %in% names(args_list)) {
          args_list$ytol <- stats::qchisq(level, df = 1) / 2 + 2
        }
        args_list$obj <- ad_obj
        args_list$name <- idx
        args_list$trace <- trace
        
        prof <- tryCatch(suppressWarnings(do.call(TMB::tmbprofile, args_list)), error = function(e) NULL)
        if (!is.null(prof)) {
          prof_list[[p_name]] <- prof
          if (show_plot) {
            plot(prof, main = paste("Profile Likelihood:", p_name))
          }
          ci <- tryCatch(suppressWarnings(confint(prof, level = level)), error = function(e) NULL)
          if (!is.null(ci)) {
            # Transform CIs back to constrained scale
            transform_val <- function(u) {
              if (is.na(u)) return(NA)
              u_full <- self$par_unc
              u_full[idx_fix_active[idx]] <- u

              c_list <- to_constrained(unconstrained_vector_to_list(u_full, self$model$par_list), self$model$par_list)

              # Extract the specific value for p_name
              p_base <- gsub("\\[.*\\]$", "", p_name)
              p_val <- as.numeric(c_list[[p_base]])
              
              # Find the index of p_name within the flattened names of p_base
              p_info <- self$model$par_list[[p_base]]
              p_names_def <- self$model$par_names[[p_base]]
              all_p_names <- generate_flat_names(p_base, p_info$dim, p_names_def)
              
              p_idx <- which(all_p_names == p_name)
              if (length(p_idx) == 0) p_idx <- 1
              
              return(p_val[p_idx])
            }
            low_c <- transform_val(ci[1])
            up_c <- transform_val(ci[2])
            res_mat[i, ] <- c(min(low_c, up_c, na.rm = TRUE), max(low_c, up_c, na.rm = TRUE))
          }
        }
      }
      if (!quiet) cat("\nDone.\n")
      
      # 4. Summary
      res_df <- data.frame(
        variable = target_names,
        Estimate = format(round(self$df_fixed[target_names, "Estimate"], digits), nsmall = digits),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      res_df[[low_label]] <- format(round(res_mat[, 1], digits), nsmall = digits)
      res_df[[up_label]] <- format(round(res_mat[, 2], digits), nsmall = digits)
      
      rownames(res_df) <- NULL
      
      if (!quiet) {
        cat("\nProfile Likelihood Confidence Intervals:\n")
        # Format for pretty printing: left align the variable names
        res_df_print <- res_df
        res_df_print$variable <- format(res_df_print$variable, justify = "left")
        print(res_df_print, row.names = FALSE, right = FALSE)
        cat("\n")
      }
      
      attr(res_df, "profiles") <- prof_list
      return(invisible(res_df))
    }
  )
)

#' @export
summary.map_fit <- function(object, ...) {
  object$summary(...)
}

#' @export
print.map_fit <- function(x, ...) {
  x$summary()
  invisible(x)
}
