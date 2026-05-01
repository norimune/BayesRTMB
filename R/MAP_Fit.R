#' MAP fit object
#'
#' An R6 class storing optimization results from maximum a posteriori
#' (MAP) estimation.
#
#' @field model The `RTMB_Model` object used for estimation.
#' @field par_vec Parameter vector on the unconstrained scale (constrained values unlisted).
#' @field par Parameter list on the constrained scale.
#' @field par_unc Parameter vector on the unconstrained scale (raw unconstrained values).
#' @field ci_method Method used for CI estimation ("wald", "profile", or "sampling").
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
#' @field map List; the parameter mapping used.
#'
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

    #' @description Get point estimate for a target parameter (internal use).
    #' @param target Target parameter name.
    #' @return Matrix or array of point estimate.
    get_point_estimate = function(target) {
      if (!is.null(self$par[[target]])) return(self$par[[target]])
      if (!is.null(self$transform[[target]])) return(self$transform[[target]])
      if (!is.null(self$generate[[target]])) return(self$generate[[target]])
      stop("Parameter not found: ", target)
    },

    #' @description Return point estimates (EAP is not applicable).
    #' @param ... Ignored.
    #' @return A named list of point estimates.
    EAP = function(...) {
      warning("EAP is not applicable for MAP_Fit. Returning point estimates instead.")
      return(c(self$par, self$transform, self$generate))
    },

    #' @description Return point estimates (MAP sampling method is not applicable).
    #' @param ... Ignored.
    #' @return A named list of point estimates.
    MAP = function(...) {
      warning("Sampling-based MAP is not applicable for MAP_Fit. Returning point estimates instead.")
      return(c(self$par, self$transform, self$generate))
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
    #' @param ci_method Method used for CI estimation ("wald" or "sampling").
    #' @param laplace Logical; whether Laplace approximation was used.
    #' @param map List; the parameter mapping used.
    initialize = function(model, par_vec, par, objective, log_ml, convergence, sd_rep, df_fixed,
                          random_effects, df_transform = NULL, df_generate = NULL, opt_history = NULL,
                          transform = NULL, generate = NULL, se_samples = NULL, par_unc = NULL, 
                          ci_method = "wald", laplace = TRUE, map = NULL) {
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
      self$laplace <- laplace
      self$map <- map

      if (!is.null(self$par_vec)) self$par_vec <- Re(self$par_vec)
      if (!is.null(self$par)) self$par <- lapply(self$par, Re)
      if (!is.null(self$random_effects) && !is.data.frame(self$random_effects)) {
        self$random_effects <- lapply(self$random_effects, Re)
      }
      if (!is.null(self$transform)) self$transform <- lapply(self$transform, Re)
      if (!is.null(self$generate)) self$generate <- lapply(self$generate, Re)
      class(self) <- c(class(self), "RTMB_Fit_Base")
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

    #' @description Summarize MAP estimates.
    #' @param pars Character vector specifying the names of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param max_rows Maximum number of rows to print in summaries. Default is 10.
    #' @param digits Number of digits to print.
    #' @param ranef Logical; whether to also display random effect estimates. Default is FALSE.
    #' @return A summary object, typically a data frame.
    summary = function(pars = NULL, max_rows = 10, digits = 5, ranef = FALSE) {
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

      all_dfs <- list()
      if (!is.null(self$df_fixed) && nrow(self$df_fixed) > 0) all_dfs$fixed <- self$df_fixed
      if (isTRUE(ranef) && is.data.frame(self$random_effects) && nrow(self$random_effects) > 0) {
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
        base_names <- gsub("\\[.*\\]$", "", rownames(df_combined))
        keep_idx <- rownames(df_combined) %in% pars | base_names %in% pars
        df_combined <- df_combined[keep_idx, , drop = FALSE]
      }

      if (nrow(df_combined) == 0) {
        cat("\nNo matching parameters found.\n")
        return(invisible(df_combined))
      }

      if (nrow(df_combined) > 0) {
        var_names <- rownames(df_combined)
        base_names <- gsub("\\[.*\\]$", "", var_names)

        target_views <- c()
        if (!is.null(self$model$view)) {
          target_views <- c(target_views, self$model$view)
        }

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
        if (self$ci_method == "profile") ci_label <- "95% Profile Likelihood CI"
        if (self$ci_method == "sampling") ci_label <- "95% Sampling-based CI"
      }
      cat(sprintf("\nPoint Estimates and %s:\n", ci_label))

      num_cols <- sapply(df_combined, is.numeric)
      df_combined[num_cols] <- lapply(df_combined[num_cols], function(x) {
        x[abs(x) < 1e-12 & !is.na(x)] <- 0
        return(x)
      })

      if ("DF" %in% colnames(df_combined)) {
        df_combined$DF <- ifelse(is.infinite(df_combined$DF), "Inf", as.character(round(df_combined$DF)))
      }

      out_df <- data.frame(variable = rownames(df_combined), df_combined, check.names = FALSE, stringsAsFactors = FALSE)
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

    #' @description Calculate Profile Likelihood confidence intervals for specific parameters.
    #' @param pars Character vector of parameter names to profile. If NULL, all fixed parameters are profiled.
    #' @param level Confidence level (default is 0.95).
    #' @param trace Logical; whether to print profiling progress. Default is FALSE.
    #' @param digits Integer; number of decimal places to print. Default is 5.
    #' @param show_plot Logical; whether to plot the profile likelihood curves. Default is FALSE.
    #' @param quiet Logical; whether to suppress text output. Default is FALSE.
    #' @param ... Additional arguments passed to TMB::tmbprofile (e.g., ytol).
    #' @return A data frame containing the profile-based confidence intervals, with the raw profile objects stored in the "profiles" attribute.
    profile = function(pars = NULL, level = 0.95, trace = FALSE, digits = 5, show_plot = FALSE, quiet = FALSE, ...) {
      if (!quiet) cat("Estimating confidence intervals via Profile Likelihood...\n")
      
      # 1. Re-build ad_obj using optimized values and no jacobian adjustment
      ad_setup <- self$model$build_ad_obj(init = self$par, laplace = self$laplace, 
                                          map = self$map, jacobian_target = "none")
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
            res_mat[i, ] <- ci
          }
        }
      }
      if (!quiet) cat("\nDone.\n")
      
      # 4. Summary
      res_df <- data.frame(
        Parameter = target_names,
        Estimate = format(round(self$df_fixed[target_names, "Estimate"], digits), nsmall = digits),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      res_df[[low_label]] <- format(round(res_mat[, 1], digits), nsmall = digits)
      res_df[[up_label]] <- format(round(res_mat[, 2], digits), nsmall = digits)
      
      rownames(res_df) <- NULL
      
      if (!quiet) {
        cat("\nProfile Likelihood Confidence Intervals:\n")
        print(res_df, digits = digits)
        cat("\n")
      }
      
      attr(res_df, "profiles") <- prof_list
      return(invisible(res_df))
    },

    #' @description Plot the posterior density for specific parameters based on Profile Likelihood.
    #' @param pars Character vector of parameter names to plot.
    #' @param level Confidence level for the range to plot (default is 0.999).
    #' @param ... Additional arguments passed to profile() (e.g. ytol).
    plot_density = function(pars = NULL, level = 0.999, ...) {
      # 1. Get profile data
      # Use show_plot = FALSE to avoid deviance plots
      res_ci <- self$profile(pars = pars, level = level, show_plot = FALSE, quiet = TRUE, ...)
      profiles <- attr(res_ci, "profiles")
      
      if (length(profiles) == 0) return(invisible(NULL))
      
      ad_setup <- self$model$build_ad_obj(init = self$par, laplace = self$laplace, 
                                          map = self$map, jacobian_target = "none")
      ad_obj <- ad_setup$ad_obj
      
      n_pars <- length(profiles)
      if (n_pars > 1) {
        old_par <- graphics::par(mfrow = grDevices::n2mfrow(n_pars))
        on.exit(graphics::par(old_par))
      }
      
      for (p_name in names(profiles)) {
        prof <- profiles[[p_name]]
        u_vals <- prof[[1]]
        deviance <- prof$value
        lik <- exp(-0.5 * deviance)
        
        c_vals <- numeric(length(u_vals))
        dens <- numeric(length(u_vals))
        
        all_fixed_names <- rownames(self$df_fixed)
        idx_in_active <- which(all_fixed_names == p_name)
        
        all_indices <- seq_along(ad_obj$env$last.par)
        fixed_indices <- if (is.null(ad_obj$env$random)) all_indices else all_indices[-ad_obj$env$random]
        idx_full <- fixed_indices[idx_in_active]
        
        opt_unc_list <- to_unconstrained(self$par, self$model$par_list)
        
        for (j in seq_along(u_vals)) {
          eps <- 1e-5
          calc_c <- function(u) {
            tmp_unc_list <- opt_unc_list
            curr_idx <- 0
            for (name in names(tmp_unc_list)) {
              L <- length(tmp_unc_list[[name]])
              if (idx_full > curr_idx && idx_full <= curr_idx + L) {
                tmp_unc_list[[name]][idx_full - curr_idx] <- u
                break
              }
              curr_idx <- curr_idx + L
            }
            tmp_con <- to_constrained(tmp_unc_list, self$model$par_list)
            curr_idx_c <- 0
            for (name in names(tmp_con)) {
              L_c <- length(tmp_con[[name]])
              if (idx_full > curr_idx_c && idx_full <= curr_idx_c + L_c) {
                return(as.numeric(tmp_con[[name]])[idx_full - curr_idx_c])
              }
              curr_idx_c <- curr_idx_c + L_c
            }
            return(NA)
          }
          
          phi <- calc_c(u_vals[j])
          phi_plus <- calc_c(u_vals[j] + eps)
          dphi_du <- (phi_plus - phi) / eps
          
          c_vals[j] <- phi
          dens[j] <- lik[j] / abs(dphi_du)
        }
        
        valid <- !is.na(c_vals) & !is.na(dens) & !is.infinite(dens)
        c_vals <- c_vals[valid]; dens <- dens[valid]
        
        if (length(c_vals) < 2) {
          warning(sprintf("Not enough valid points for '%s'.", p_name))
          next
        }
        
        dx <- diff(c_vals); mid_dens <- (dens[-1] + dens[-length(dens)]) / 2
        area <- sum(dx * mid_dens)
        if (!is.na(area) && area > 0) dens <- dens / area
        
        graphics::plot(c_vals, dens, type = "l", lwd = 2, col = "#2c3e50",
                       xlab = p_name, ylab = "Posterior Density",
                       main = paste("Posterior Distribution:", p_name))
        graphics::grid()
        graphics::polygon(c(c_vals, rev(c_vals)), c(dens, rep(0, length(dens))),
                          col = "#3498db33", border = NA)
        graphics::abline(v = self$df_fixed[p_name, "Estimate"], col = "#e74c3c", lty = 2, lwd = 1.5)
      }
      return(invisible(NULL))
    }
  )
)
