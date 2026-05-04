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

    #' @description Create a new `Classic_Model` object.
    #' @param type Character string specifying the model type (e.g., "lm", "lmer").
    #' @param formula The formula used for the model.
    #' @param data The data used for estimation.
    #' @param family The distribution family.
    #' @param view Parameter names to prioritize in summary display.
    #' @param obj Optional underlying model object (e.g., RTMB_Model for lmer).
    #' @param refit_fn Optional function to re-fit the model on new data.
    initialize = function(type, formula, data, family = "gaussian", view = NULL, obj = NULL, refit_fn = NULL) {
      self$type <- type
      self$formula <- formula
      self$data <- data
      self$family <- family
      self$view <- view
      self$obj <- obj
      self$refit_fn <- refit_fn
    },

    #' @description Perform estimation using classical methods.
    #' @param bootstrap Logical; whether to perform non-parametric bootstrap.
    #' @param n_boot Integer; number of bootstrap samples.
    #' @return A `Classic_Fit` object containing the results.
    estimate = function(bootstrap = FALSE, n_boot = 1000) {
      # 1. Perform original fit
      fit_res <- self$.perform_fit(self$data)
      
      if (!bootstrap) {
        return(Classic_Fit$new(self, fit_res))
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
      
      res <- Classic_Fit$new(self, res_df)
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
        
        grp_var <- as.character(bars[[1]][[3]])
        clusters <- unique(self$data[[grp_var]])
        resampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
        
        resampled_list <- lapply(resampled_clusters, function(c) {
          self$data[self$data[[grp_var]] == c, , drop = FALSE]
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
        return(stats::lm(self$formula, data = data))
      } else if (self$type %in% c("lmer", "lmm", "glmm")) {
        obj <- self$obj
        
        fe_params <- c()
        if ("Intercept" %in% names(obj$par_list)) fe_params <- c(fe_params, "Intercept")
        if ("b" %in% names(obj$par_list)) fe_params <- c(fe_params, "b")

        ad_setup <- obj$build_ad_obj(laplace = TRUE, jacobian_target = "none")
        ad_obj <- ad_setup$ad_obj
        opt <- nlminb(ad_obj$par, ad_obj$fn, ad_obj$gr)
        
        sd_rep <- RTMB::sdreport(ad_obj)
        smry_ran <- summary(sd_rep, select = "random")
        ran_names <- rownames(smry_ran)
        fe_idx_in_ran <- which(gsub("\\[.*\\]$", "", ran_names) %in% fe_params)
        
        reml_res <- obj$calculate_reml_satterthwaite_df(ad_obj, opt$par, fe_idx_in_ran)
        
        df_combined <- data.frame(
          Estimate = smry_ran[fe_idx_in_ran, "Estimate"],
          `Std. Error` = reml_res$se,
          df = reml_res$df,
          check.names = FALSE
        )
        
        # Rename parameters using obj$par_names
        raw_names <- ran_names[fe_idx_in_ran]
        new_names <- character(length(raw_names))
        b_count <- 0
        for (i in seq_along(raw_names)) {
          if (raw_names[i] == "Intercept" || grepl("^Intercept\\[", raw_names[i])) {
            new_names[i] <- "Intercept"
          } else if (raw_names[i] == "b" || grepl("^b\\[", raw_names[i])) {
            b_count <- b_count + 1
            # Try to get index from regex first, e.g., b[1]
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
        
        # Final safety check for duplicate names
        if (any(duplicated(new_names))) {
           new_names <- make.unique(new_names)
        }
        rownames(df_combined) <- new_names
        
        smry_fix <- summary(sd_rep, select = "fixed")
        
        # 1. sigma
        sig_idx <- grepl("^sigma[0-9]*$", rownames(smry_fix))
        if (any(sig_idx)) {
           sig_smry <- smry_fix[sig_idx, , drop = FALSE]
           df_combined <- rbind(df_combined, data.frame(
             Estimate = exp(sig_smry[, "Estimate"]),
             `Std. Error` = NA, df = NA, row.names = rownames(sig_smry), check.names = FALSE
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

        df_combined$`t value` <- df_combined$Estimate / df_combined$`Std. Error`
        df_combined$Pr <- 2 * pt(-abs(df_combined$`t value`), df = df_combined$df)
        df_combined$`Lower 95%` <- df_combined$Estimate + qt(0.025, df = df_combined$df) * df_combined$`Std. Error`
        df_combined$`Upper 95%` <- df_combined$Estimate + qt(0.975, df = df_combined$df) * df_combined$`Std. Error`
        
        return(df_combined)
        
      } else if (self$type == "corr") {
        cor_mat <- cor(data, use = "pairwise.complete.obs")
        P <- ncol(data); var_names <- colnames(data)
        N <- nrow(data)
        df_list <- list()
        for (i in 1:(P-1)) {
          for (j in (i+1):P) {
            r <- cor_mat[i, j]
            z <- 0.5 * log((1 + r) / (1 - r))
            se_z <- 1 / sqrt(N - 3)
            ci_lower_z <- z - 1.96 * se_z
            ci_upper_z <- z + 1.96 * se_z
            ci_lower <- (exp(2 * ci_lower_z) - 1) / (exp(2 * ci_lower_z) + 1)
            ci_upper <- (exp(2 * ci_upper_z) - 1) / (exp(2 * ci_upper_z) + 1)
            
            df_val <- N - 2
            t_val <- r * sqrt(df_val / (1 - r^2))
            p_val <- 2 * pt(-abs(t_val), df = df_val)
            
            p_name <- paste0("corr[", var_names[i], ", ", var_names[j], "]")
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
        return(do.call(rbind, df_list))
      } else if (self$type == "mediation") {
        obj <- self$obj
        obj$data <- data
        ad_setup <- obj$build_ad_obj(jacobian_target = "none")
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
                  new_names[idx[j]] <- if (label == "Intercept") p else paste0(p, "[", label, "]")
                }
              }
            }
          }
        }
        rownames(df_res) <- new_names
        
        df_res$df <- nrow(data) - 1
        df_res$`t value` <- df_res$Estimate / df_res$`Std. Error`
        df_res$Pr <- 2 * pt(-abs(df_res$`t value`), df = df_res$df)
        df_res$`Lower 95%` <- df_res$Estimate + qt(0.025, df = df_res$df) * df_res$`Std. Error`
        df_res$`Upper 95%` <- df_res$Estimate + qt(0.975, df = df_res$df) * df_res$`Std. Error`
        
        return(df_res[order(rownames(df_res)), ])
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
#' @field bootstrap_results A matrix containing bootstrap samples, if applicable.
#'
#' @export
Classic_Fit <- R6::R6Class(
  classname = "Classic_Fit",
  public = list(
    model = NULL,
    fit = NULL,
    par = NULL,
    bootstrap_results = NULL,

    #' @description Create a new `Classic_Fit` object.
    #' @param model The `Classic_Model` object.
    #' @param fit The result of the estimation.
    initialize = function(model, fit) {
      self$model <- model
      self$fit <- fit
      self$par <- self$.construct_par_list(fit)
    },

    #' @description Display a summary of the estimation results.
    #' @param digits Number of digits to print.
    summary = function(digits = 5) {
      cat("\nCall:\n")
      cat(paste("Classical estimation via", self$model$type), "\n")
      if (!is.null(self$bootstrap_results)) {
        cat("Based on non-parametric bootstrap (", nrow(self$bootstrap_results), " samples)\n")
      }
      cat("\nPoint Estimates and Confidence Intervals:\n")
      
      if (is.data.frame(self$fit)) {
        df_print <- self$fit
        if (!is.null(df_print$Pr)) {
          sig <- symnum(df_print$Pr, corr = FALSE, na = FALSE,
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols = c("***", "**", "*", ".", " "))
          df_print$sig <- as.character(sig)
        }
        print(round(df_print[ , !names(df_print) %in% "sig", drop=FALSE], digits))
      } else {
        print(self$fit)
      }
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
