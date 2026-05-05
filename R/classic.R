#' Classic Model Class for Frequentist Estimation
#'
#' @description
#' An R6 class representing a classical (frequentist) statistical model.
#' This class is used when `classic = TRUE` is specified in wrapper functions
#' like `rtmb_lm`.
#'
#' @field type Character string specifying the model type (e.g., "lm").
#' @field formula The formula used for the model.
#' @field data The data used for estimation.
#' @field family The distribution family.
#' @field view Parameter names to prioritize in summary display.
#' @field obj Optional underlying model object (e.g., RTMB_Model for lmer).
#' @field refit_fn Optional function to re-fit the model on new data.
#' @field extra Optional list of additional parameters for specific model types.
#'
#' @export
Classic_Model <- R6::R6Class(
  classname = "Classic_Model",
  public = list(
    type = NULL,
    formula = NULL,
    data = NULL,
    family = NULL,
    view = NULL,
    obj = NULL,
    refit_fn = NULL,
    extra = list(),

    #' @description Create a new `Classic_Model` object.
    #' @param type Character string specifying the model type (e.g., "lm", "lmer").
    #' @param formula The formula used for the model.
    #' @param data The data used for estimation.
    #' @param family The distribution family.
    #' @param view Parameter names to prioritize in summary display.
    #' @param obj Optional underlying model object (e.g., RTMB_Model for lmer).
    #' @param refit_fn Optional function to re-fit the model on new data.
    #' @param extra Optional list of additional parameters for specific model types.
    initialize = function(type, formula, data, family = "gaussian", view = NULL, obj = NULL, refit_fn = NULL, extra = list()) {
      self$type <- type
      self$formula <- formula
      self$data <- data
      self$family <- family
      self$view <- view
      self$obj <- obj
      self$refit_fn <- refit_fn
      self$extra <- extra
    },

    #' @description Perform estimation using classical methods.
    #' @param bootstrap Logical; whether to perform non-parametric bootstrap.
    #' @param n_boot Integer; number of bootstrap samples.
    #' @return A `Classic_Fit` object containing the results.
    estimate = function(bootstrap = FALSE, n_boot = 1000) {
      # 1. Perform original fit
      fit_res <- self$.perform_fit(self$data)
      
      if (!bootstrap) {
        return(Classic_Fit$new(self, fit_res$df_combined, vcov = fit_res$V_beta))
      }

      # 2. Perform Bootstrap
      cat(sprintf("Performing non-parametric bootstrap (%d samples)...\n", n_boot))
      
      boot_results <- list()
      pb <- txtProgressBar(min = 0, max = n_boot, style = 3)
      
      # Suppress messages during bootstrap
      old_opt <- options(BayesRTMB.silent = TRUE)
      on.exit(options(old_opt), add = TRUE)
      
      for (i in 1:n_boot) {
        resampled_data <- self$.resample_data()
        
        tryCatch({
          # Capture output to keep it clean
          capture.output({
            boot_fit <- if (!is.null(self$refit_fn)) {
               self$refit_fn(resampled_data)
            } else {
               self$.perform_fit(resampled_data)
            }
          })
          
          if (is.data.frame(boot_fit)) {
            boot_results[[length(boot_results) + 1]] <- boot_fit$Estimate
          } else if (inherits(boot_fit, "Classic_Fit")) {
            boot_results[[length(boot_results) + 1]] <- boot_fit$fit$Estimate
          } else if (inherits(boot_fit, "lm")) {
            boot_results[[length(boot_results) + 1]] <- stats::coef(boot_fit)
          }
        }, error = function(e) {
          # Skip failed fits
        })
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      if (length(boot_results) == 0) stop("All bootstrap fits failed.")
      
      # Aggregate results
      boot_mat <- do.call(rbind, boot_results)
      
      boot_se <- apply(boot_mat, 2, sd, na.rm = TRUE)
      boot_lower <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
      boot_upper <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
      
      # Prepare result dataframe
      if (is.data.frame(fit_res)) {
        res_df <- fit_res
      } else if (inherits(fit_res, "lm")) {
        s_lm <- summary(fit_res)
        res_df <- data.frame(
          Estimate = stats::coef(fit_res),
          `Std. Error` = s_lm$coefficients[, 2],
          check.names = FALSE
        )
      } else {
        res_df <- data.frame(Estimate = boot_results[[1]], check.names = FALSE)
      }

      # Update with bootstrap stats
      n_rows <- min(nrow(res_df), length(boot_se))
      res_df$`Std. Error`[1:n_rows] <- boot_se[1:n_rows]
      res_df$`Lower 95%`[1:n_rows] <- boot_lower[1:n_rows]
      res_df$`Upper 95%`[1:n_rows] <- boot_upper[1:n_rows]
      res_df$`t value`[1:n_rows] <- res_df$Estimate[1:n_rows] / res_df$`Std. Error`[1:n_rows]
      if (!is.null(res_df$df)) {
        res_df$Pr[1:n_rows] <- 2 * pt(-abs(res_df$`t value`[1:n_rows]), df = res_df$df[1:n_rows])
      }
      
      res <- Classic_Fit$new(self, res_df, vcov = fit_res$V_beta)
      res$bootstrap_results <- boot_mat
      return(res)
    },

    #' @description (Internal) Resample data for bootstrap.
    .resample_data = function() {
      N <- nrow(self$data)
      
      if (self$type %in% c("lmer", "lmm", "glmm")) {
        bars <- findbars(self$formula)
        if (is.null(bars)) {
          idx <- sample(1:N, N, replace = TRUE)
          return(self$data[idx, , drop = FALSE])
        }
        
        # Cluster bootstrap: sample by the highest level cluster
        # (the one with the fewest groups)
        num_groups <- sapply(bars, function(b) length(unique(self$data[[as.character(b[[3]])]])))
        best_bar_idx <- which.min(num_groups)
        
        grp_var <- as.character(bars[[best_bar_idx]][[3]])
        clusters <- unique(self$data[[grp_var]])
        resampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
        
        # Fast resampling using split and lapply
        split_data <- split(self$data, self$data[[grp_var]])
        resampled_list <- lapply(resampled_clusters, function(c) {
          split_data[[as.character(c)]]
        })
        return(do.call(rbind, resampled_list))
      } else {
        idx <- sample(1:N, N, replace = TRUE)
        return(self$data[idx, , drop = FALSE])
      }
    },

    #' @description (Internal) Perform a single fit.
    #' @param data The data to fit.
    .perform_fit = function(data) {
      if (self$type == "lm") {
        fit <- stats::lm(self$formula, data = data)
        return(list(df_combined = fit, V_beta = vcov(fit)))
      } else if (self$type %in% c("lmer", "lmm", "glmm")) {
        obj <- self$obj
        
        fe_params <- c()
        fixed_names <- c()
        if ("Intercept" %in% names(obj$par_list)) {
           fe_params <- c(fe_params, "Intercept")
           fixed_names <- c(fixed_names, "Intercept")
        }
        if ("Intercept_c" %in% names(obj$par_list)) {
           fe_params <- c(fe_params, "Intercept_c")
           fixed_names <- c(fixed_names, "Intercept") # Display as Intercept
        }
        if ("b" %in% names(obj$par_list)) {
           fe_params <- c(fe_params, "b")
           b_names <- if (!is.null(obj$par_names$b)) obj$par_names$b else paste0("b[", 1:obj$par_list$b$length, "]")
           fixed_names <- c(fixed_names, b_names)
        }
        K <- length(fixed_names)

        ad_setup <- obj$build_ad_obj(laplace = TRUE, jacobian_target = "none")
        ad_obj <- ad_setup$ad_obj
        opt <- nlminb(ad_obj$par, ad_obj$fn, ad_obj$gr)
        
        # Use sdreport with joint precision to get covariance of random-treated fixed effects
        sd_rep <- RTMB::sdreport(ad_obj, getJointPrecision = TRUE)
        
        # Get the full joint covariance if precision is available
        V_joint <- if (!is.null(sd_rep$jointPrecision)) {
          solve(sd_rep$jointPrecision)
        } else {
          sd_rep$cov.fixed # Fallback (should not happen for classic=TRUE)
        }
        
        smry_ran <- summary(sd_rep, select = "random")
        ran_names <- rownames(smry_ran)
        
        # Match fixed effects by prefix (e.g., "b[1]" matches "b")
        fe_idx_in_ran <- which(gsub("\\[.*\\]$", "", ran_names) %in% fe_params)
        
        # The joint covariance includes all random effects. 
        # Fixed effects treated as random are in fe_idx_in_ran.
        V_beta <- V_joint[fe_idx_in_ran, fe_idx_in_ran, drop = FALSE]
        colnames(V_beta) <- rownames(V_beta) <- fixed_names
        
        reml_res <- obj$calculate_reml_satterthwaite_df(ad_obj, opt$par, fe_idx_in_ran)
        
        # Rename parameters using obj$par_names
        raw_names <- ran_names[fe_idx_in_ran]
        new_names <- character(length(raw_names))
        b_count <- 0
        for (i in seq_along(raw_names)) {
          if (raw_names[i] == "Intercept" || grepl("^Intercept\\[", raw_names[i])) {
            new_names[i] <- "Intercept"
          } else if (raw_names[i] == "b" || grepl("^b\\[", raw_names[i])) {
            b_count <- b_count + 1
            p_match <- regmatches(raw_names[i], regexec("^b\\[([0-9]+)\\]$", raw_names[i]))[[1]]
            idx <- if (length(p_match) > 1) as.numeric(p_match[2]) else b_count
            if (!is.null(obj$par_names$b) && idx <= length(obj$par_names$b)) {
              new_names[i] <- paste0("b[", obj$par_names$b[idx], "]")
            } else {
              new_names[i] <- raw_names[i]
            }
          } else {
            new_names[i] <- raw_names[i]
          }
        }
        if (any(duplicated(new_names))) new_names <- make.unique(new_names)

        df_combined <- data.frame(
          Estimate = smry_ran[fe_idx_in_ran, "Estimate"],
          `Std. Error` = reml_res$se,
          df = reml_res$df,
          row.names = new_names,
          check.names = FALSE
        )
        
        smry_fix <- summary(sd_rep, select = "fixed")
        
        # 1. sigma
        sig_idx <- which(rownames(smry_fix) == "sigma")
        if (length(sig_idx) > 0) {
           sig_smry <- smry_fix[sig_idx, , drop = FALSE]
           new_sig_names <- rep("sigma", length(sig_idx))
           if (!is.null(obj$par_names$sigma)) {
              for (j in seq_along(new_sig_names)) {
                if (j <= length(obj$par_names$sigma)) {
                   new_sig_names[j] <- paste0("sigma[", obj$par_names$sigma[j], "]")
                }
              }
           }
           if (length(new_sig_names) > 1 && any(duplicated(new_sig_names))) {
             new_sig_names <- make.unique(new_sig_names)
           }
           
           df_combined <- rbind(df_combined, data.frame(
             Estimate = exp(sig_smry[, "Estimate"]),
             `Std. Error` = NA, df = NA, row.names = new_sig_names, check.names = FALSE
           ))
        }
        
        # 2. sd
        sd_idx <- grepl("^sd[0-9]*$", rownames(smry_fix))
        if (any(sd_idx)) {
          sd_smry <- smry_fix[sd_idx, , drop = FALSE]
          sd_names <- rownames(sd_smry)
          new_sd_names <- sd_names
          for (p in names(obj$par_names)) {
             if (grepl("^sd", p)) {
                labels <- obj$par_names[[p]]
                idx_in_sd <- which(sd_names == p)
                for (j in seq_along(idx_in_sd)) {
                   if (j <= length(labels)) {
                      new_sd_names[idx_in_sd[j]] <- paste0("sd[", labels[j], "]")
                   }
                }
             }
          }
          df_combined <- rbind(df_combined, data.frame(
             Estimate = exp(sd_smry[, "Estimate"]),
             `Std. Error` = NA, df = NA, row.names = new_sd_names, check.names = FALSE
          ))
        }
        
        # 3. Other positive parameters
        for (p in c("shape", "phi", "nu")) {
          if (p %in% rownames(smry_fix)) {
             df_combined <- rbind(df_combined, data.frame(
               Estimate = exp(smry_fix[p, "Estimate"]),
               `Std. Error` = NA, df = NA, row.names = p, check.names = FALSE
             ))
          }
        }
        
        # 4. Residual correlation (interval -1 to 1)
        if ("rho_resid" %in% rownames(smry_fix)) {
           rho_idx <- which(rownames(smry_fix) == "rho_resid")
           rho_smry <- smry_fix[rho_idx, , drop = FALSE]
           rho_est <- (1 / (1 + exp(-rho_smry[, "Estimate"]))) * 2 - 1
           
           # Handle multiple lags (TOEP)
           new_rho_names <- if (nrow(rho_smry) > 1) {
              paste0("rho_resid[lag", 1:nrow(rho_smry), "]")
           } else {
              "rho_resid"
           }
           
           df_combined <- rbind(df_combined, data.frame(
             Estimate = rho_est,
             `Std. Error` = NA, df = NA, row.names = new_rho_names, check.names = FALSE
           ))
        }

        df_combined$`t value` <- df_combined$Estimate / df_combined$`Std. Error`
        df_combined$Pr <- 2 * pt(-abs(df_combined$`t value`), df = df_combined$df)
        df_combined$`Lower 95%` <- df_combined$Estimate + qt(0.025, df = df_combined$df) * df_combined$`Std. Error`
        df_combined$`Upper 95%` <- df_combined$Estimate + qt(0.975, df = df_combined$df) * df_combined$`Std. Error`
        
        return(list(df_combined = df_combined, V_beta = V_beta))
        
      } else if (self$type == "corr") {
        P_y <- if (!is.null(self$extra$P_y)) self$extra$P_y else ncol(data)
        P_x <- if (!is.null(self$extra$P_x)) self$extra$P_x else 0
        method <- if (!is.null(self$extra$method)) self$extra$method else "pearson"
        P <- ncol(data)
        var_names <- colnames(data)
        target_names <- var_names[1:P_y]
        N <- nrow(data)
        
        df_list <- list()
        full_cor_mat <- cor(data, use = "pairwise.complete.obs", method = method)
        
        if (P_x > 0) {
           # Partial correlation derivation
           R_yy <- full_cor_mat[1:P_y, 1:P_y, drop = FALSE]
           R_yx <- full_cor_mat[1:P_y, (P_y+1):P, drop = FALSE]
           R_xx <- full_cor_mat[(P_y+1):P, (P_y+1):P, drop = FALSE]
           
           # Solve for partial covariance
           P_cov <- R_yy - R_yx %*% MASS::ginv(R_xx) %*% t(R_yx)
           # Convert to correlation
           D <- diag(1 / sqrt(pmax(diag(P_cov), 1e-12)))
           cor_mat <- D %*% P_cov %*% D
           N_eff <- N - P_x
        } else {
           cor_mat <- full_cor_mat
           N_eff <- N
        }
        
        if (P_y < 2) stop("Correlation requires at least 2 target variables.")
        
        for (i in 1:(P_y-1)) {
          for (j in (i+1):P_y) {
            r <- cor_mat[i, j]
            # Fisher z-transformation for CI and SE
            z <- 0.5 * log((1 + r) / (1 - r))
            se_z <- 1 / sqrt(N_eff - 3)
            ci_lower_z <- z - 1.96 * se_z
            ci_upper_z <- z + 1.96 * se_z
            ci_lower <- (exp(2 * ci_lower_z) - 1) / (exp(2 * ci_lower_z) + 1)
            ci_upper <- (exp(2 * ci_upper_z) - 1) / (exp(2 * ci_upper_z) + 1)
            
            df_val <- N_eff - 2
            t_val <- r * sqrt(df_val / (1 - r^2))
            p_val <- 2 * pt(-abs(t_val), df = df_val)
            
            label_type <- if (P_x > 0) "pcorr" else "corr"
            p_name <- paste0(label_type, "[", var_names[i], ", ", var_names[j], "]")
            df_list[[p_name]] <- data.frame(
              Estimate = r,
              `Std. Error` = NA,
              `Lower 95%` = ci_lower,
              `Upper 95%` = ci_upper,
              `t value` = t_val,
              df = df_val,
              Pr = p_val,
              row.names = p_name,
              check.names = FALSE
            )
          }
        }
        df_combined <- do.call(rbind, df_list)
        return(list(df_combined = df_combined, V_beta = NULL))
      } else if (self$type == "ttest") {
        Y1 <- self$extra$Y1
        Y2 <- self$extra$Y2
        var.equal <- if (!is.null(self$extra$var.equal)) self$extra$var.equal else FALSE
        paired <- if (!is.null(self$extra$paired)) self$extra$paired else FALSE
        
        # Standard t-test
        fit <- stats::t.test(Y1, Y2, var.equal = var.equal, paired = paired)
        
        # Calculate Cohen's d (Estimate for 'delta')
        if (paired) {
          diffs <- Y1 - Y2
          n <- length(diffs)
          # For paired t-test, d is mean(diff) / sd(diff)
          d_est <- mean(diffs, na.rm = TRUE) / stats::sd(diffs, na.rm = TRUE)
          se_d <- sqrt(1/n + d_est^2 / (2*n))
          diff_est <- as.numeric(fit$estimate)
        } else {
          n1 <- length(Y1); n2 <- length(Y2)
          v1 <- stats::var(Y1); v2 <- stats::var(Y2)
          s_pooled <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
          # Align with t.test: mean(x) - mean(y)
          d_est <- (mean(Y1) - mean(Y2)) / s_pooled
          se_d <- sqrt((n1 + n2) / (n1 * n2) + d_est^2 / (2 * (n1 + n2)))
          diff_est <- as.numeric(fit$estimate[1] - fit$estimate[2])
        }
        
        crit_t <- stats::qt(0.975, df = fit$parameter)
        
        df_combined <- data.frame(
          Estimate = d_est,
          `Std. Error` = se_d,
          `Lower 95%` = d_est - crit_t * se_d,
          `Upper 95%` = d_est + crit_t * se_d,
          `t value` = fit$statistic,
          df = fit$parameter,
          Pr = fit$p.value,
          row.names = "delta",
          check.names = FALSE
        )
        
        # Add the raw difference as well
        df_combined <- rbind(df_combined, data.frame(
          Estimate = diff_est,
          `Std. Error` = fit$stderr,
          `Lower 95%` = fit$conf.int[1],
          `Upper 95%` = fit$conf.int[2],
          `t value` = fit$statistic,
          df = fit$parameter,
          Pr = fit$p.value,
          row.names = "diff",
          check.names = FALSE
        ))
        
        return(list(df_combined = df_combined, V_beta = NULL))
      } else if (self$type == "mediation") {
        obj <- self$obj
        obj$data <- data
        ad_setup <- obj$build_ad_obj(jacobian_target = "none", include_generate = TRUE)
        ad_obj <- ad_setup$ad_obj
        opt <- nlminb(ad_obj$par, ad_obj$fn, ad_obj$gr)
        sd_rep <- RTMB::sdreport(ad_obj)
        smry_fix <- summary(sd_rep, select = "fixed")
        smry_gen <- summary(sd_rep, select = "report")
        
        df_res <- as.data.frame(rbind(smry_fix, smry_gen))
        colnames(df_res) <- c("Estimate", "Std. Error")
        
        sig_idx <- grepl("^sigma[0-9]*$", rownames(df_res))
        if (any(sig_idx)) {
           df_res[sig_idx, "Estimate"] <- exp(df_res[sig_idx, "Estimate"])
        }
        
        new_names <- rownames(df_res)
        if (!is.null(obj$par_names)) {
          for (p in names(obj$par_names)) {
            idx <- which(grepl(paste0("^", p, "($|\\[|\\.)"), new_names))
            if (length(idx) > 0) {
              labels <- obj$par_names[[p]]
              for (j in seq_along(idx)) {
                if (j <= length(labels)) {
                  label <- labels[j]
                  new_names[idx[j]] <- paste0(p, "[", label, "]")
                }
              }
            }
          }
        }
        rownames(df_res) <- new_names
        
        df_res$df <- nrow(data) - 1
        if (!is.null(self$extra$df_map)) {
          df_map <- self$extra$df_map
          for (p in names(df_map)) {
            idx <- which(grepl(paste0("^", p, "($|\\[|\\.)"), rownames(df_res)))
            if (length(idx) > 0) df_res$df[idx] <- df_map[[p]]
          }
        }
        
        df_res$`t value` <- df_res$Estimate / df_res$`Std. Error`
        df_res$Pr <- 2 * pt(-abs(df_res$`t value`), df = df_res$df)
        df_res$`Lower 95%` <- df_res$Estimate + qt(0.025, df = df_res$df) * df_res$`Std. Error`
        df_res$`Upper 95%` <- df_res$Estimate + qt(0.975, df = df_res$df) * df_res$`Std. Error`
        
        # Reorder based on self$view if present
        if (!is.null(self$view)) {
          v_order <- self$view
          row_names <- rownames(df_res)
          new_idx <- integer(0)
          
          # Helper to match parameter names (with or without [labels])
          for (v in v_order) {
            # Match "b1" or "b1[perf]"
            matches <- which(row_names == v | grepl(paste0("^", v, "\\["), row_names))
            # Filter out already added indices
            matches <- setdiff(matches, new_idx)
            if (length(matches) > 0) {
              new_idx <- c(new_idx, matches)
            }
          }
          # Add any remaining variables (not in view_order) at the end
          remaining_idx <- setdiff(seq_len(nrow(df_res)), new_idx)
          if (length(remaining_idx) > 0) {
            new_idx <- c(new_idx, remaining_idx[order(row_names[remaining_idx])])
          }
          df_res <- df_res[new_idx, ]
        } else {
          df_res <- df_res[order(rownames(df_res)), ]
        }
        return(df_res)
      }
    }
  )
)

#' Classic fit object
#'
#' @description
#' An R6 class representing the results of a classical (frequentist) estimation.
#'
#' @field model The `Classic_Model` object used for estimation.
#' @field fit The result of the estimation (dataframe or lm object).
#' @field par A named list of parameter estimates.
#' @field vcov Variance-covariance matrix of fixed effects.
#' @field bootstrap_results A matrix containing bootstrap samples, if applicable.
#'
#' @export
Classic_Fit <- R6::R6Class(
  classname = "Classic_Fit",
  public = list(
    model = NULL,
    fit = NULL,
    par = NULL,
    vcov = NULL,
    bootstrap_results = NULL,

    #' @description Create a new `Classic_Fit` object.
    #' @param model The `Classic_Model` object.
    #' @param fit The result of the estimation.
    #' @param vcov Variance-covariance matrix of fixed effects.
    initialize = function(model, fit, vcov = NULL) {
      self$model <- model
      self$fit <- fit
      self$vcov <- vcov
      self$par <- self$.construct_par_list(fit)
    },

    #' @description Display a summary of the estimation results.
    #' @param digits Number of digits to print for estimates.
    summary = function(digits = 5) {
      cat("\nCall:\n")
      cat(paste("Classical estimation via", self$model$type), "\n")
      if (!is.null(self$bootstrap_results)) {
        cat("Based on non-parametric bootstrap (", nrow(self$bootstrap_results), " samples)\n")
      }
      cat("\nPoint Estimates and Confidence Intervals:\n")
      if (self$model$type == "ttest") {
        levs <- self$model$extra$levs
        cat(sprintf("(Comparison: %s - %s)\n", levs[1], levs[2]))
      }
      
      if (is.data.frame(self$fit) || inherits(self$fit, "lm")) {
        # Convert lm to dataframe if needed
        if (inherits(self$fit, "lm")) {
          s_lm <- summary(self$fit)
          df_print <- as.data.frame(s_lm$coefficients)
          colnames(df_print) <- c("Estimate", "Std. Error", "t value", "Pr")
          df_print$df <- self$fit$df.residual
          ci <- confint(self$fit)
          df_print$`Lower 95%` <- ci[, 1]
          df_print$`Upper 95%` <- ci[, 2]
          rownames(df_print)[rownames(df_print) == "(Intercept)"] <- "Intercept"
        } else {
          df_print <- self$fit
        }
        
        # --- Internal to Requested Contrast Transformation for Display ---
        req_ct <- if (!is.null(self$model$obj$requested_contrasts)) self$model$obj$requested_contrasts else "treatment"
        curr_ct <- if (!is.null(self$model$obj$contrasts)) self$model$obj$contrasts else "sum"
        
        if (req_ct != curr_ct) {
          # Use names for safe identification of fixed effects
          beta_full <- df_print$Estimate
          names(beta_full) <- rownames(df_print)
          V_full <- self$vcov
          
          fe_idx <- which(grepl("^(Intercept|Intercept_c|b\\[)", names(beta_full)))
          if (length(fe_idx) > 0) {
            beta <- beta_full[fe_idx]
            V <- V_full[fe_idx, fe_idx]
            
            formula <- nobars(self$model$formula)
            predictor_vars <- all.vars(delete.response(terms(formula)))
            relevant_data <- self$model$data[, predictor_vars, drop = FALSE]
            levs <- lapply(relevant_data, function(x) if(is.factor(x)) levels(x) else unique(x))
            grid <- expand.grid(levs)
            
            # X_from (current internal: usually sum)
            old_opts <- options(contrasts = if (curr_ct == "sum") c("contr.sum", "contr.poly") else c("contr.treatment", "contr.poly"))
            X_from <- model.matrix(delete.response(terms(formula)), grid)
            options(old_opts)
            
            # X_to (requested: usually treatment)
            old_opts <- options(contrasts = if (req_ct == "treatment") c("contr.treatment", "contr.poly") else c("contr.sum", "contr.poly"))
            X_to <- model.matrix(delete.response(terms(formula)), grid)
            options(old_opts)
            
            if (ncol(X_from) == length(beta)) {
              # Use MASS::ginv for robust transformation even with rank-deficient models
              M <- MASS::ginv(X_to) %*% X_from
              beta_new <- as.numeric(M %*% beta)
              V_new <- M %*% V %*% t(M)
              se_new <- sqrt(diag(V_new))
              
              # Update df_print with new values for fixed effects
              df_print$Estimate[fe_idx] <- beta_new
              df_print$`Std. Error`[fe_idx] <- se_new
              df_print$`t value`[fe_idx] <- beta_new / se_new
              if (!is.null(df_print$Pr)) {
                df_print$Pr[fe_idx] <- 2 * pt(-abs(df_print$`t value`[fe_idx]), df = df_print$df[fe_idx])
              }
              if ("Lower 95%" %in% names(df_print)) {
                df_print$`Lower 95%`[fe_idx] <- beta_new - 1.96 * se_new
                df_print$`Upper 95%`[fe_idx] <- beta_new + 1.96 * se_new
              }
              # Update Row Names to match requested contrast coding
              new_names <- colnames(X_to)
              new_names[new_names == "(Intercept)"] <- "Intercept"
              # Wrap other fixed effects in b[...]
              is_fe <- new_names != "Intercept"
              new_names[is_fe] <- paste0("b[", new_names[is_fe], "]")
              
              rownames(df_print)[fe_idx] <- new_names
            }
          }
        }
        
        # 1. Round main numeric columns
        cols_to_round <- setdiff(names(df_print), c("Pr", "df"))
        for (col in cols_to_round) {
          if (is.numeric(df_print[[col]])) {
            df_print[[col]] <- round(df_print[[col]], digits)
          }
        }
        
        # 2. Round df to 1 decimal place
        if ("df" %in% names(df_print)) {
          df_print$df <- round(df_print$df, 1)
        }
        
        # 3. Format Pr
        if (!is.null(df_print$Pr)) {
          # Standard significance symbols
          sig <- symnum(df_print$Pr, corr = FALSE, na = FALSE,
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols = c("***", "**", "*", ".", " "))
          
          # Format Pr values (e.g., .000, .012)
          formatted_pr <- sapply(df_print$Pr, function(p) {
            if (is.na(p)) return("")
            if (is.character(p)) return(p) # Already formatted
            if (p < 0.001) return(".000")
            res <- format(round(p, 3), nsmall = 3)
            if (startsWith(res, "0")) return(substring(res, 2)) # Remove leading zero
            return(res)
          })
          
          df_print$Pr <- formatted_pr
          df_print$sig <- as.character(sig)
        }
        
        # 4. Reorder columns for consistency
        desired_order <- c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "t value", "df", "Pr")
        current_names <- names(df_print)
        
        # Filter existing columns
        order_to_use <- intersect(desired_order, current_names)
        other_names <- setdiff(current_names, c(desired_order, "sig"))
        
        # Construct final selection
        final_cols <- c(order_to_use, other_names)
        if ("sig" %in% current_names) {
          final_cols <- c(final_cols, "sig")
        }
        
        df_final <- df_print[, final_cols, drop = FALSE]
        
        # Rename sig to blank just before printing
        if ("sig" %in% names(df_final)) {
           names(df_final)[names(df_final) == "sig"] <- ""
        }
        
        # Print the data frame cleanly
        print(df_final, quote = FALSE, right = TRUE)

        # --- Enhanced Output for LM/GLM ---
        if (inherits(self$fit, "lm")) {
          s_lm <- summary(self$fit)
          cat("\n---\n")
          cat(sprintf("Residual standard error: %s on %d degrees of freedom\n", 
                      format(round(s_lm$sigma, digits), nsmall = digits), s_lm$df[2]))
          cat(sprintf("Multiple R-squared: %s, Adjusted R-squared: %s\n", 
                      format(round(s_lm$r.squared, 4), nsmall = 4), 
                      format(round(s_lm$adj.r.squared, 4), nsmall = 4)))
          if (!is.null(s_lm$fstatistic)) {
            f <- s_lm$fstatistic
            p_f <- pf(f[1], f[2], f[3], lower.tail = FALSE)
            cat(sprintf("F-statistic: %s on %d and %d DF, p-value: %s\n", 
                        format(round(f[1], 2), nsmall = 2), f[2], f[3], 
                        if (p_f < 0.001) "< .001" else format(round(p_f, 4), nsmall = 4)))
          }
        }
      } else {
        print(self$fit)
      }
    },

    #' @description Perform ANOVA (Wald F-tests) on the fitted model.
    #' @param method Character; "reml" (standard) or "ls" (experimental).
    #' @param type Integer; Type of Sum of Squares (only Type III supported currently).
    #' @return A data frame containing the ANOVA table.
    anova = function(method = c("reml", "ls"), type = 3) {
      # Ensure consistent contrasts for term matching if needed
      old_opts <- NULL
      if (!is.null(self$model$contrasts)) {
        if (self$model$contrasts == "sum") {
          old_opts <- options(contrasts = c("contr.sum", "contr.poly"))
        } else if (self$model$contrasts == "treatment") {
          old_opts <- options(contrasts = c("contr.treatment", "contr.poly"))
        }
      }
      if (!is.null(old_opts)) on.exit(options(old_opts), add = TRUE)

      method <- match.arg(method)
      if (is.null(self$model$extra$X_assign)) stop("ANOVA requires term assignments (X_assign).")
      V_full <- self$vcov
      beta_full <- if (is.data.frame(self$fit)) self$fit$Estimate else stats::coef(self$fit)
      
      # Handle names for lm objects specifically
      full_names <- if (is.data.frame(self$fit)) rownames(self$fit) else names(beta_full)
      if (inherits(self$fit, "lm")) {
        full_names[full_names == "(Intercept)"] <- "Intercept"
      }
      names(beta_full) <- full_names
      
      assign_idx <- self$model$extra$X_assign
      
      # Match fixed effects based on model type
      if (inherits(self$fit, "lm")) {
        # For lm, all parameters are fixed effects
        fe_idx <- seq_along(beta_full)
      } else {
        # For RTMB models, they start with Intercept or b
        fe_idx <- which(grepl("^(Intercept|Intercept_c|b($|\\[))", names(beta_full)))
      }
      
      if (length(fe_idx) == 0) stop("Could not identify fixed effects for ANOVA.")
      
      beta <- beta_full[fe_idx]
      V <- V_full[fe_idx, fe_idx]
      
      # full_assign for fixed effects (Intercept=0, others 1..T)
      full_assign <- c(0, assign_idx)
      
      # --- Automatic Contrast Transformation to 'sum' ---
      # Type III ANOVA requires orthogonal (sum) contrasts for meaningful main effects.
      ct_setting <- if (!is.null(self$model$obj$contrasts)) self$model$obj$contrasts else "sum"
      if (ct_setting != "sum") {
        formula <- nobars(self$model$formula)
        # Use a minimal grid of unique PREDICTOR combinations for perfect M calculation
        predictor_vars <- all.vars(delete.response(terms(formula)))
        relevant_data <- self$model$data[, predictor_vars, drop = FALSE]
        levs <- lapply(relevant_data, function(x) if(is.factor(x)) levels(x) else unique(x))
        grid <- expand.grid(levs)
        
        # 1. X with current contrasts
        old_opts <- options(contrasts = if (ct_setting == "treatment") 
                            c("contr.treatment", "contr.poly") else options()$contrasts)
        X_curr <- model.matrix(delete.response(terms(formula)), grid)
        options(old_opts)
        
        # 2. X with sum contrasts
        old_opts <- options(contrasts = c("contr.sum", "contr.poly"))
        X_sum <- model.matrix(delete.response(terms(formula)), grid)
        options(old_opts)
        
        # 3. Calculate transformation matrix M using ginv for robustness
        if (ncol(X_curr) == length(beta)) {
          M <- MASS::ginv(X_sum) %*% X_curr
          beta <- as.numeric(M %*% beta)
          names(beta) <- colnames(X_sum) 
          V <- M %*% V %*% t(M)
          # Update assign_idx to match the new contrast coding
          assign_idx <- attr(X_sum, "assign")
        }
      }

      terms <- self$model$extra$X_terms
      res_list <- list()
      
      # Mapping coefficients to terms
      has_int <- any(grepl("^Intercept", names(beta))) || any(full_assign == 0)
      
      # Handle Intercept separately or include in terms
      # Users typically don't want Intercept in ANOVA table for mixed models
      all_assign <- if (has_int) seq_along(terms) else seq_along(terms)
      all_term_names <- terms

      for (i in seq_along(all_assign)) {
        a_id <- all_assign[i]
        t_name <- all_term_names[i]
        
        # Find indices in beta/V corresponding to this term
        # assign_idx matches the order in model.matrix (which includes intercept as 0)
        idx <- which(assign_idx == a_id)
        if (length(idx) == 0) next
        
        L <- matrix(0, nrow = length(idx), ncol = length(beta))
        for (j in seq_along(idx)) L[j, idx[j]] <- 1
        
        # Wald statistic: F = (L*beta)' * inv(L*V*L') * (L*beta) / rank(L)
        LVL <- L %*% V %*% t(L)
        inv_LVL <- try(solve(LVL), silent = TRUE)
        if (inherits(inv_LVL, "try-error")) inv_LVL <- MASS::ginv(LVL)
        
        W <- as.numeric(t(L %*% beta) %*% inv_LVL %*% (L %*% beta))
        df1 <- length(idx)
        f_val <- W / df1
        
        # Denominator DF (Satterthwaite)
        # We take the average of the DFs for the involved parameters as an approximation
        # or use the first one if they are all from the same level.
        if (is.data.frame(self$fit) && !is.null(self$fit$df)) {
           df2 <- min(self$fit$df[idx], na.rm = TRUE)
        } else if (inherits(self$fit, "lm")) {
           df2 <- self$fit$df.residual
        } else {
           df2 <- Inf
        }
        
        p_val <- pf(f_val, df1, df2, lower.tail = FALSE)
        
        res_list[[t_name]] <- data.frame(
          `NumDF` = df1,
          `DenDF` = df2,
          `F value` = f_val,
          `Pr(>F)` = p_val,
          row.names = t_name,
          check.names = FALSE
        )
      }
      
      res_df <- do.call(rbind, res_list)
      
      # Add significance symbols
      sig <- symnum(res_df$`Pr(>F)`, corr = FALSE, na = FALSE,
                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                    symbols = c("***", "**", "*", ".", " "))
      res_df$signif <- as.character(sig)

      # --- Enhanced ANOVA Output with R2 ---
      class(res_df) <- c("anova", "data.frame")
      heading <- "ANOVA Table (Wald F-tests)"
      if (inherits(self$fit, "lm")) {
        s_lm <- summary(self$fit)
        heading <- paste0(heading, sprintf("\nMultiple R-squared: %s, Adjusted R-squared: %s", 
                          format(round(s_lm$r.squared, 4), nsmall = 4), 
                          format(round(s_lm$adj.r.squared, 4), nsmall = 4)))
      }
      attr(res_df, "heading") <- heading
      return(res_df)
    },

    #' @description Calculate Least Squares Means (Marginal Means).
    #' @param specs Character vector of factors to calculate means for.
    #' @return A data frame containing the marginal means, SEs, and CIs.
    lsmeans = function(specs) {
      # Ensure consistent contrasts for reference grid construction
      old_opts <- NULL
      if (!is.null(self$model$contrasts)) {
        if (self$model$contrasts == "sum") {
          old_opts <- options(contrasts = c("contr.sum", "contr.poly"))
        } else if (self$model$contrasts == "treatment") {
          old_opts <- options(contrasts = c("contr.treatment", "contr.poly"))
        }
      }
      if (!is.null(old_opts)) on.exit(options(old_opts), add = TRUE)

      if (is.null(self$model$data)) stop("lsmeans requires the original data.")
      data <- self$model$data
      formula <- self$model$formula
      beta <- if (is.data.frame(self$fit)) self$fit$Estimate else stats::coef(self$fit)
      V <- self$vcov
      
      # 1. Identify all factors in the model
      vars <- all.vars(nobars(formula))[-1]
      is_cat <- sapply(data[vars], function(x) is.factor(x) || is.character(x))
      cat_vars <- vars[is_cat]
      
      if (!all(specs %in% cat_vars)) stop("specs must be categorical factors in the model.")
      
      # 2. Create a reference grid for all categorical factors
      grid_list <- lapply(data[cat_vars], function(x) levels(as.factor(x)))
      ref_grid <- expand.grid(grid_list)
      
      # 3. Handle continuous covariates by setting them to their mean
      cont_vars <- vars[!is_cat]
      for (v in cont_vars) {
        ref_grid[[v]] <- mean(data[[v]], na.rm = TRUE)
      }
      
      # 4. Generate model matrix for the reference grid
      # We need to use the same terms and contrast settings as the original fit
      mf_orig <- model.frame(nobars(formula), data)
      orig_terms <- delete.response(terms(mf_orig))
      
      # Explicitly construct contrasts list based on model preference
      ct_setting <- if (!is.null(self$model$obj$contrasts)) self$model$obj$contrasts else "sum"
      ct_list <- list()
      # Only apply to factors present in the variables
      for (v in cat_vars) {
        if (ct_setting == "treatment") {
          ct_list[[v]] <- "contr.treatment"
        } else {
          ct_list[[v]] <- "contr.sum"
        }
      }
      
      # Use raw model matrix columns to ensure exact parameter matching
      X_grid_raw <- model.matrix(orig_terms, ref_grid, contrasts.arg = ct_list)
      
      # Match columns with the original model (handling dropped coefficients)
      orig_cols <- self$model$extra$X_colnames
      if (!is.null(orig_cols)) {
        # Intercept handling
        grid_cols <- colnames(X_grid_raw)
        grid_cols[grid_cols == "(Intercept)"] <- "Intercept" # Internal consistency
        
        col_idx <- match(orig_cols, colnames(X_grid_raw))
        # Handle cases where "(Intercept)" vs "Intercept" naming differs
        if (any(is.na(col_idx)) && "Intercept" %in% orig_cols) {
           col_idx[orig_cols == "Intercept"] <- which(colnames(X_grid_raw) == "(Intercept)")
        }
        
        if (any(is.na(col_idx))) {
           stop("Could not match lsmeans grid columns with original model parameters.")
        }
        X_grid <- X_grid_raw[, col_idx, drop = FALSE]
      } else {
        X_grid <- X_grid_raw
      }
      
      # Positional matching is more robust than name matching for RTMB fixed effects
      beta_vals <- if (is.data.frame(self$fit)) self$fit$Estimate else stats::coef(self$fit)
      fe_idx <- which(grepl("^(Intercept|Intercept_c|b\\[)", rownames(self$fit)))
      
      # Prediction matrix X_grid must match the order and count of beta_match
      beta_match <- beta_vals[fe_idx]
      V_match <- self$vcov[fe_idx, fe_idx]
      
      if (ncol(X_grid) != length(beta_match)) {
         stop("Dimension mismatch in lsmeans: X_grid columns (", ncol(X_grid), 
              ") do not match fixed effects (", length(beta_match), ").")
      }
      # Target factors are in 'specs'
      groups <- interaction(ref_grid[specs], drop = TRUE, sep = ":")
      unique_groups <- levels(groups)
      
      res_list <- list()
      for (grp in unique_groups) {
        idx <- which(groups == grp)
        # Average the rows of X_grid for this group
        L <- colMeans(X_grid[idx, , drop = FALSE])
        
        est <- as.numeric(L %*% beta_match)
        se <- sqrt(as.numeric(t(L) %*% V_match %*% L))
        
        # DF (approximate)
        # Try to find a representative DF for the specified factors
        df_val <- Inf
        if (is.data.frame(self$fit) && !is.null(self$fit$df)) {
          # Look for fixed effect terms (b[...]) that match any of the specs
          match_pattern <- paste0("^b\\[(", paste(specs, collapse="|"), ")")
          match_idx <- grepl(match_pattern, rownames(self$fit))
          if (any(match_idx)) {
            # Use the minimum DF among relevant terms (conservative)
            df_val <- min(self$fit$df[match_idx], na.rm = TRUE)
          } else {
            # Fallback to Intercept DF or mean
            if ("Intercept" %in% rownames(self$fit)) {
              df_val <- self$fit["Intercept", "df"]
            } else {
              df_val <- mean(self$fit$df, na.rm = TRUE)
            }
          }
        } else if (inherits(self$fit, "lm")) {
          df_val <- self$fit$df.residual
        }
        
        res_list[[grp]] <- data.frame(
          lsmean = est,
          `Std. Error` = se,
          df = df_val,
          `Lower 95%` = est + qt(0.025, df_val) * se,
          `Upper 95%` = est + qt(0.975, df_val) * se,
          row.names = grp,
          check.names = FALSE
        )
      }
      
      return(do.call(rbind, res_list))
    },

    #' @description Print the fit results.
    print = function() {
      self$summary()
    },

    #' @description (Internal) Construct a list of parameters from the fit.
    #' @param fit The fit result (dataframe or lm object).
    .construct_par_list = function(fit) {
      par_list <- list()
      if (inherits(fit, "lm")) {
        coefs <- stats::coef(fit)
        par_list$Intercept <- coefs["(Intercept)"]
        par_list$b <- coefs[names(coefs) != "(Intercept)"]
        return(par_list)
      }
      if (is.data.frame(fit)) {
        est <- fit$Estimate; names(est) <- rownames(fit)
        par_list$Intercept <- est["Intercept"]
        par_list$b <- est[grepl("^b\\[", names(est))]
        par_list$sd <- est[grepl("^sd\\[", names(est))]
        par_list$sigma <- est[grepl("^sigma", names(est))]
        if (self$model$type == "corr") {
          par_list$corr <- est[grepl("^corr\\[", names(est))]
        }
        par_list$IE <- est[grepl("^IE_", names(est))]
      }
      return(par_list)
    }
  )
)
