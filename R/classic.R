#' Classic Model Class for Frequentist Estimation
#'
#' @description
#' An R6 class representing a classical (frequentist) statistical model.
#' This class is used when `classic = TRUE` is specified in wrapper functions
#' like `rtmb_lm`. It bypasses the RTMB automatic differentiation engine
#' and uses standard R functions (e.g., `stats::lm`) for estimation.
#'
#' @export
Classic_Model <- R6::R6Class(
  classname = "Classic_Model",
  public = list(
    #' @field type Character string specifying the model type (e.g., "lm").
    type = NULL,
    #' @field formula The formula used for the model.
    formula = NULL,
    #' @field data The data used for estimation.
    data = NULL,
    #' @field family The distribution family.
    family = NULL,
    #' @field view Parameter names to prioritize in summary display.
    view = NULL,
    #' @field obj Optional underlying model object (e.g., RTMB_Model for lmer).
    obj = NULL,

    #' @description Create a new `Classic_Model` object.
    #' @param type Character string specifying the model type (e.g., "lm").
    #' @param formula The formula used for the model.
    #' @param data The data used for estimation.
    #' @param family The distribution family.
    #' @param view Parameter names to prioritize in summary display.
    #' @param obj Optional underlying model object.
    initialize = function(type, formula, data, family = "gaussian", view = NULL, obj = NULL) {
      self$type <- type
      self$formula <- formula
      self$data <- data
      self$family <- family
      self$view <- view
      self$obj <- obj
    },

    #' @description Perform estimation using classical methods.
    #' @return A `Classic_Fit` object containing the results.
    estimate = function() {
      if (self$type == "lm") {
        if (!identical(self$family, "gaussian")) {
          stop("Classic mode ('classic = TRUE') for linear models only supports the Gaussian family.")
        }

        cat("Estimating via Ordinary Least Squares (OLS)...\n")
        # Ensure stats is used
        fit <- stats::lm(self$formula, data = self$data)
        return(Classic_Fit$new(self, fit))

      } else if (self$type == "lmer") {
        # REML estimation using RTMB
        obj <- self$obj
        if (is.null(obj)) stop("RTMB model object is missing for lmer estimation.")

        # Identify fixed effects
        # We assume the wrapper has already prepared the obj with appropriate random flags
        fe_params <- c()
        if ("Intercept" %in% names(obj$par_list)) fe_params <- c(fe_params, "Intercept")
        if ("Intercept_c" %in% names(obj$par_list)) fe_params <- c(fe_params, "Intercept_c")
        if ("b" %in% names(obj$par_list)) fe_params <- c(fe_params, "b")
        if ("beta_raw" %in% names(obj$par_list)) fe_params <- c(fe_params, "beta_raw")

        # Build and Optimize
        # Note: In the wrapper, we already set the random flags for fe_params.
        ad_setup <- obj$build_ad_obj(laplace = TRUE, jacobian_target = "none")
        ad_obj <- ad_setup$ad_obj
        
        cat("Estimating via Restricted Maximum Likelihood (REML) using RTMB...\n")
        opt <- nlminb(ad_obj$par, ad_obj$fn, ad_obj$gr)
        
        # Extract Estimates and SEs
        sd_rep <- RTMB::sdreport(ad_obj)
        smry_ran <- summary(sd_rep, select = "random")
        smry_fix <- summary(sd_rep, select = "fixed")
        
        ran_names <- rownames(smry_ran)
        fe_idx_in_ran <- which(gsub("\\[.*\\]$", "", ran_names) %in% fe_params)
        
        # Satterthwaite Degrees of Freedom and Conditional SEs
        reml_res <- obj$calculate_reml_satterthwaite_df(ad_obj, opt$par, fe_idx_in_ran)
        fe_dfs <- reml_res$df
        se_fe <- reml_res$se
        
        # Construct Result Data Frame
        est_fe <- smry_ran[fe_idx_in_ran, "Estimate"]
        t_fe <- est_fe / se_fe
        p_fe <- 2 * pt(-abs(t_fe), df = fe_dfs)
        
        df_combined <- data.frame(
          Estimate = est_fe,
          `Std. Error` = se_fe,
          `Lower 95%` = est_fe + qt(0.025, df = fe_dfs) * se_fe,
          `Upper 95%` = est_fe + qt(0.975, df = fe_dfs) * se_fe,
          `t value` = t_fe,
          df = fe_dfs,
          Pr = p_fe,
          check.names = FALSE
        )
        
        # Rename fixed effects
        fe_row_names <- ran_names[fe_idx_in_ran]
        b_counter <- 1
        # We need access to fixed_names. It should be in obj$par_names$b or similar.
        fixed_names <- obj$par_names$b %||% obj$par_names$beta_raw
        
        for (i in seq_along(fe_row_names)) {
          nm <- fe_row_names[i]
          if (nm %in% c("Intercept", "Intercept_c")) {
            fe_row_names[i] <- "Intercept"
            next
          }
          if (nm == "b" || nm == "beta_raw" || grepl("^b\\[", nm)) {
            if (!is.null(fixed_names) && b_counter <= length(fixed_names)) {
              fe_row_names[i] <- paste0("b[", fixed_names[b_counter], "]")
              b_counter <- b_counter + 1
            }
          }
        }
        rownames(df_combined) <- fe_row_names
        
        # Add Variance Components
        vc_idx <- which(grepl("sigma|sd|corr", rownames(smry_fix)))
        if (length(vc_idx) > 0) {
          vc_names <- rownames(smry_fix)[vc_idx]
          vc_est_raw <- smry_fix[vc_idx, "Estimate"]
          vc_est <- ifelse(grepl("sigma|sd", vc_names), exp(vc_est_raw), vc_est_raw)
          
          # Improve names using par_names metadata if available
          new_vc_names <- vc_names
          if (!is.null(obj$par_names)) {
            for (i in seq_along(new_vc_names)) {
              vnm <- new_vc_names[i]
              if (grepl("^sd(\\d+)?$", vnm)) {
                idx_str <- gsub("sd", "", vnm)
                idx <- if (idx_str == "") 1 else as.integer(idx_str)
                if (!is.null(obj$par_names$sd) && !is.na(idx) && idx <= length(obj$par_names$sd)) {
                  new_vc_names[i] <- paste0("sd[", obj$par_names$sd[idx], "]")
                }
              } else if (grepl("^corr(\\d+)?$", vnm)) {
                idx_str <- gsub("corr", "", vnm)
                idx <- if (idx_str == "") 1 else as.integer(idx_str)
                if (!is.null(obj$par_names$corr) && !is.na(idx) && idx <= length(obj$par_names$corr)) {
                  new_vc_names[i] <- paste0("corr[", obj$par_names$corr[idx], "]")
                }
              }
            }
          }
          
          df_vc <- data.frame(
            Estimate = vc_est,
            `Std. Error` = NA,
            `Lower 95%` = NA,
            `Upper 95%` = NA,
            `t value` = NA,
            df = NA,
            Pr = NA,
            check.names = FALSE
          )
          rownames(df_vc) <- new_vc_names
          df_combined <- rbind(df_combined, df_vc)
        }
        
        return(Classic_Fit$new(self, df_combined))

      } else if (self$type == "corr") {
        # For corr, the calculation is already performed and stored in self$data by the wrapper
        return(Classic_Fit$new(self, self$data))
      } else {
        stop("Unknown classic model type: ", self$type)
      }
    }
  )
)

#' Classic Fit Class for Frequentist Estimation Results
#'
#' @description
#' An R6 class for storing and summarizing results from a `Classic_Model`.
#'
#' @export
Classic_Fit <- R6::R6Class(
  classname = "Classic_Fit",
  public = list(
    #' @field model The `Classic_Model` object.
    model = NULL,
    #' @field fit The raw fit object (e.g., an `lm` object).
    fit = NULL,

    #' @description Create a new `Classic_Fit` object.
    #' @param model The `Classic_Model` object.
    #' @param fit The raw fit object (e.g., an `lm` object).
    initialize = function(model, fit) {
      self$model <- model
      self$fit <- fit
    },

    #' @description Display a summary of the estimation results.
    #' @param pars Character vector specifying the names of parameters to summarize.
    #' @param max_rows Maximum number of rows to print.
    #' @param digits Number of digits to print.
    #' @param ... Additional arguments.
    summary = function(pars = NULL, max_rows = 10, digits = 5, ...) {
      df_combined <- NULL
      if (self$model$type == "lm") {
        s_lm <- summary(self$fit)
        coef_mat <- s_lm$coefficients
        ci_mat <- stats::confint(self$fit)

        df_combined <- data.frame(
          Estimate = coef_mat[, 1],
          `Std. Error` = coef_mat[, 2],
          `Lower 95%` = ci_mat[, 1],
          `Upper 95%` = ci_mat[, 2],
          `t value` = coef_mat[, 3],
          df = s_lm$df[2],
          Pr = coef_mat[, 4],
          check.names = FALSE
        )

        # Consistent naming with BayesRTMB
        raw_names <- rownames(coef_mat)
        new_names <- character(length(raw_names))
        for (i in seq_along(raw_names)) {
          if (raw_names[i] == "(Intercept)") {
            new_names[i] <- "Intercept"
          } else {
            new_names[i] <- paste0("b[", raw_names[i], "]")
          }
        }
        rownames(df_combined) <- new_names

        # Add sigma
        sigma_val <- s_lm$sigma
        df_combined <- rbind(df_combined, data.frame(
           Estimate = sigma_val,
           `Std. Error` = NA,
           `Lower 95%` = NA,
           `Upper 95%` = NA,
           `t value` = NA,
           df = NA,
           Pr = NA,
           row.names = "sigma",
           check.names = FALSE
        ))
        
        cat("\nCall:\nClassical Estimation (Frequentist)\n")
        
      } else if (self$model$type == "lmer") {
        df_combined <- self$fit
        cat("\nCall:\nClassical Mixed Model (Frequentist REML via RTMB)\n")

      } else if (self$model$type == "corr") {
        df_combined <- self$fit
        cat("\nCall:\nClassical Correlation (Frequentist via Fisher's Z)\n")
      }
      
      if (!is.null(df_combined)) {
        # Reorder columns: Estimate, SE, Lower 95%, Upper 95%, t value, df, Pr
        target_cols <- c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "t value", "df", "Pr")
        existing_cols <- intersect(target_cols, colnames(df_combined))
        other_cols <- setdiff(colnames(df_combined), target_cols)
        df_combined <- df_combined[, c(existing_cols, other_cols), drop = FALSE]

        # Add significance markers
        if ("Pr" %in% colnames(df_combined)) {
          p_vals <- df_combined[["Pr"]]
          sig <- ifelse(is.na(p_vals), "", 
                        ifelse(p_vals < .001, "***",
                               ifelse(p_vals < .01, "**",
                                      ifelse(p_vals < .05, "*",
                                             ifelse(p_vals < .1, ".", "")))))
          df_combined$sig <- sig
          # Use a single space for column name to look like it's attached to p-value
          colnames(df_combined)[colnames(df_combined) == "sig"] <- " "
        }
        if (!is.null(pars)) {
          var_names <- rownames(df_combined)
          base_names <- gsub("\\[.*\\]$", "", var_names)
          match_idx <- which(var_names %in% pars | base_names %in% pars)
          if (length(match_idx) > 0) df_combined <- df_combined[match_idx, , drop = FALSE]
        }

        cat(sprintf("\nPoint Estimates and 95%% Classic CI:\n"))

        # Apply view order if available
        if (!is.null(self$model$view) && is.null(pars)) {
          var_names <- rownames(df_combined)
          base_names <- gsub("\\[.*\\]$", "", var_names)

          priority_idx <- c()
          for (v in self$model$view) {
            priority_idx <- c(priority_idx, which(base_names == v))
          }
          priority_idx <- unique(priority_idx)
          other_idx <- setdiff(seq_len(nrow(df_combined)), priority_idx)
          df_combined <- df_combined[c(priority_idx, other_idx), , drop = FALSE]
        }

        num_cols <- sapply(df_combined, is.numeric)
        df_combined[num_cols] <- lapply(df_combined[num_cols], function(x) {
          x[abs(x) < 1e-12 & !is.na(x)] <- 0
          return(x)
        })

        out_df <- cbind(data.frame(variable = rownames(df_combined), stringsAsFactors = FALSE), df_combined)
        rownames(out_df) <- NULL

        if (!is.null(max_rows) && nrow(out_df) > max_rows) {
          out_df <- head(out_df, max_rows)
        }

        class(out_df) <- c("summary_BayesRTMB", "data.frame")
        print(out_df, digits = digits)

        cat("\n")
        invisible(df_combined)
      }
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
    }
  )
)
