#' Classic fit object
#'
#' @description
#' An R6 class representing the results of a classical (frequentist) estimation.
#'
#' @field model The `RTMB_Model` object used for estimation.
#' @field fit The result of the estimation (dataframe or lm object).
#' @field par a named list of parameter estimates.
#' @field vcov Variance-covariance matrix of fixed effects.
#' @field se_method Character string specifying the method used for standard errors.
#' @field cluster Character string specifying the cluster variable name, if any.
#' @field bootstrap_results A matrix containing bootstrap samples, if applicable.
#' @field test_results List of additional test results (e.g., chisq.test).
#' @field view Character vector of parameter names to prioritize in summary.
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
    test_results = list(),
    se_method = "wald",
    cluster = NULL,
    view = NULL,

    #' @description Create a new `Classic_Fit` object.
    #' @param model The `RTMB_Model` object.
    #' @param fit The result of the estimation.
    #' @param vcov Variance-covariance matrix of fixed effects.
    #' @param se_method Character; "wald", "robust", or "bootstrap".
    #' @param cluster Character; cluster variable name.
    #' @param test_results List of additional test results (e.g., chisq.test).
    #' @param view Character vector of parameter names to prioritize in summary.
    initialize = function(model, fit, vcov = NULL, se_method = "wald", cluster = NULL, test_results = list(), view = NULL) {
      self$model <- model
      self$fit <- fit
      self$vcov <- vcov
      self$se_method <- se_method
      self$cluster <- cluster
      self$test_results <- test_results
      self$view <- view
      self$par <- self$.construct_par_list(fit)
    },

    #' @description Compute robust standard errors (sandwich estimator).
    #' @param cluster Character; variable name for clustering.
    #' @return Self.
    compute_robust = function(cluster = NULL) {
      self$se_method <- "robust"
      self$cluster <- cluster

      formula <- nobars(self$model$formula)
      if (!is.null(self$model$raw_data)) {
        dat <- as.data.frame(self$model$raw_data)
      } else {
        dat <- as.data.frame(self$model$data)
      }
      mf <- model.frame(formula, dat)
      X <- model.matrix(formula, mf)
      y <- model.response(mf)

      if (is.data.frame(self$fit)) {
        beta_all <- self$fit$Estimate
        names(beta_all) <- rownames(self$fit)
        fe_names_in_fit <- names(beta_all)[grepl("^(Intercept|Intercept_c|b\\[)", names(beta_all))]
      } else if (inherits(self$fit, "lm")) {
        beta_all <- stats::coef(self$fit)
        fe_names_in_fit <- names(beta_all)
      } else {
        stop(paste("Unsupported fit type:", class(self$fit)[1]))
      }

      if (length(fe_names_in_fit) == 0) stop("No fixed effects found in fit object.")

      X_cols <- colnames(X)
      idx_map <- integer(length(fe_names_in_fit))
      for (i in seq_along(fe_names_in_fit)) {
        fname <- fe_names_in_fit[i]
        pos <- which(X_cols == fname)
        if (length(pos) == 0) {
          if (tolower(fname) == "intercept" || fname == "intercept_c" || fname == "(intercept)") {
            pos <- which(tolower(X_cols) == "intercept" | tolower(X_cols) == "(intercept)")
          } else {
            vname <- gsub("^b\\[(.*)\\]$", "\\1", fname)
            pos <- which(X_cols == vname)
          }
        }
        if (length(pos) > 0) idx_map[i] <- pos[1]
      }

      keep <- idx_map > 0
      if (!any(keep)) stop("None of the fixed effects could be matched.")

      fe_names_in_fit <- fe_names_in_fit[keep]
      beta <- beta_all[fe_names_in_fit]
      idx_map <- idx_map[keep]
      X_subset <- X[, idx_map, drop = FALSE]

      XtX <- t(X_subset) %*% X_subset
      bread_ols <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))
      res <- as.numeric(y - X_subset %*% beta)

      if (is.null(cluster)) {
        h <- tryCatch(diag(X_subset %*% bread_ols %*% t(X_subset)), error = function(e) rep(0, length(res)))
        omega <- res^2 / (1 - h)^2
        meat <- t(X_subset) %*% (as.numeric(omega) * X_subset)
        V_final <- bread_ols %*% meat %*% bread_ols
      } else {
        if (!cluster %in% names(dat)) stop(paste("Cluster variable", cluster, "not found."))
        grp <- as.factor(dat[[cluster]])
        scores <- X_subset * as.numeric(res)
        cluster_scores <- aggregate(scores, by = list(grp), sum)[, -1, drop = FALSE]
        meat <- t(as.matrix(cluster_scores)) %*% as.matrix(cluster_scores)
        m <- length(unique(grp)); n <- nrow(X_subset); k <- ncol(X_subset)
        adj <- (m / (m - 1)) * ((n - 1) / (n - k))
        V_final <- adj * bread_ols %*% meat %*% bread_ols
      }

      self$vcov <- V_final
      colnames(self$vcov) <- rownames(self$vcov) <- fe_names_in_fit
      self$.update_fit_with_vcov()
      return(self)
    },

    #' @description Compute bootstrap standard errors.
    #' @param n_boot Integer; number of samples.
    #' @param cluster Character; clustering variable.
    #' @return Self.
    compute_bootstrap = function(n_boot = 1000, cluster = NULL) {
      self$se_method <- "bootstrap"
      cat(sprintf("Performing non-parametric bootstrap (%d samples)...\n", n_boot))

      boot_results <- vector("list", n_boot)
      pb <- txtProgressBar(min = 0, max = n_boot, style = 3)

      old_opt <- options(BayesRTMB.silent = TRUE)
      on.exit(options(old_opt), add = TRUE)

      for (i in 1:n_boot) {
        resampled_data <- self$model$.resample_data()

        tryCatch({
          capture.output({
            boot_fit <- if (!is.null(self$model$refit_fn)) {
               self$model$refit_fn(resampled_data)
            } else {
               self$model$.perform_fit(resampled_data)
            }
          })

          if (is.data.frame(boot_fit)) {
            boot_results[[i]] <- boot_fit$Estimate
          } else if (inherits(boot_fit, "Classic_Fit")) {
            boot_results[[i]] <- boot_fit$fit$Estimate
          } else if (inherits(boot_fit, "lm")) {
            boot_results[[i]] <- stats::coef(boot_fit)
          } else if (is.list(boot_fit) && !is.null(boot_fit$df_combined)) {
            # Handle list return from .perform_fit
            boot_results[[i]] <- if (is.data.frame(boot_fit$df_combined)) boot_fit$df_combined$Estimate else stats::coef(boot_fit$df_combined)
          }
        }, error = function(e) { })
        setTxtProgressBar(pb, i)
      }
      close(pb)

      boot_results <- Filter(Negate(is.null), boot_results)
      if (length(boot_results) == 0) stop("All bootstrap fits failed.")

      boot_mat <- do.call(rbind, boot_results)
      self$bootstrap_results <- boot_mat

      boot_se <- apply(boot_mat, 2, sd, na.rm = TRUE)
      boot_lower <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
      boot_upper <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

      # Update self$fit
      if (is.data.frame(self$fit)) {
        res_df <- self$fit
      } else {
        # Convert lm to df
        s_lm <- summary(self$fit)
        res_df <- data.frame(
          Estimate = stats::coef(self$fit),
          `Std. Error` = s_lm$coefficients[, 2],
          check.names = FALSE
        )
      }

      n_rows <- min(nrow(res_df), length(boot_se))
      res_df$`Std. Error`[1:n_rows] <- boot_se[1:n_rows]
      res_df$`Lower 95%`[1:n_rows] <- boot_lower[1:n_rows]
      res_df$`Upper 95%`[1:n_rows] <- boot_upper[1:n_rows]
      res_df$`t value`[1:n_rows] <- res_df$Estimate[1:n_rows] / pmax(res_df$`Std. Error`[1:n_rows], 1e-12)

      if (!is.null(res_df$df)) {
        res_df$Pr[1:n_rows] <- 2 * pt(-abs(res_df$`t value`[1:n_rows]), df = res_df$df[1:n_rows])
      }

      self$fit <- res_df
      return(self)
    },

    #' @description (Internal) Update fit data frame with current vcov.
    .update_fit_with_vcov = function() {
      if (is.null(self$vcov)) return()
      new_se <- sqrt(diag(self$vcov))

      if (is.data.frame(self$fit)) {
        n <- min(nrow(self$fit), length(new_se))
        self$fit$`Std. Error`[1:n] <- new_se[1:n]
        self$fit$`t value`[1:n] <- self$fit$Estimate[1:n] / pmax(self$fit$`Std. Error`[1:n], 1e-12)
        if (!is.null(self$fit$df)) {
          self$fit$Pr[1:n] <- 2 * pt(-abs(self$fit$`t value`[1:n]), df = self$fit$df[1:n])
        }
        # Update CIs
        crit <- if (!is.null(self$fit$df)) qt(0.975, df = self$fit$df) else 1.96
        self$fit$`Lower 95%`[1:n] <- self$fit$Estimate[1:n] - crit[1:n] * self$fit$`Std. Error`[1:n]
        self$fit$`Upper 95%`[1:n] <- self$fit$Estimate[1:n] + crit[1:n] * self$fit$`Std. Error`[1:n]
      }
      # For lm objects, we don't update the object itself,
      # but ensure summary() uses self$vcov.
    },

    #' @description Get the AIC of the fitted model.
    AIC = function() {
      ll <- self$logLik()
      k <- as.numeric(attr(ll, "df"))
      if (length(k) == 0 || is.na(k)) k <- 0
      return(-2 * as.numeric(ll) + 2 * k[1])
    },

    #' @description Get the BIC of the fitted model.
    BIC = function() {
      ll <- self$logLik()
      k <- as.numeric(attr(ll, "df"))
      n <- as.numeric(attr(ll, "nobs"))
      if (length(k) == 0 || is.na(k)) k <- 0
      if (length(n) == 0 || is.na(n) || n <= 0) n <- 1
      return(-2 * as.numeric(ll) + log(n[1]) * k[1])
    },

    #' @description Print the fit results.
    #' @param ... Additional arguments passed to `summary()`.
    print = function(...) {
      print(self$summary(...))
      invisible(self)
    },

    #' @description Get the Log-Likelihood of the fitted model.
    logLik = function() {
      if (inherits(self$fit, "lm")) {
        return(stats::logLik(self$fit))
      }

      val <- if (!is.null(self$model$extra[["loglik"]])) self$model$extra[["loglik"]] else NA
      df_val <- if (!is.null(self$model$extra[["df"]])) self$model$extra[["df"]] else nrow(self$fit)
      nobs_val <- if (!is.null(self$model$extra[["nobs"]])) self$model$extra[["nobs"]] else 0

      res <- val
      attr(res, "df") <- df_val
      attr(res, "nobs") <- nobs_val
      class(res) <- "logLik"
      return(res)
    },

    #' @description Display a summary of the estimation results.
    #' @param digits Number of digits to print for estimates.
    #' @param max_rows Maximum number of rows to display in the coefficient table.
    summary = function(digits = 5, max_rows = 10) {
      res <- list(
        type = self$model$type,
        family = self$model$family,
        se_method = self$se_method,
        cluster = self$cluster,
        bootstrap = if (!is.null(self$bootstrap_results)) nrow(self$bootstrap_results) else NULL,
        logLik = self$logLik(),
        AIC = self$AIC(),
        BIC = self$BIC(),
        extra = self$model$extra,
        test_results = self$test_results,
        digits = digits
      )

      if (is.data.frame(self$fit) || inherits(self$fit, "lm")) {
        # Convert lm to dataframe if needed
        if (inherits(self$fit, "lm")) {
          s_lm <- summary(self$fit)
          is_asymptotic <- inherits(self$model, c("rtmb_loglinear", "rtmb_table"))
          df_print <- as.data.frame(s_lm$coefficients)

          # Update with robust vcov if available
          if (!is.null(self$vcov)) {
            new_se <- sqrt(diag(self$vcov))
            n_match <- min(nrow(df_print), length(new_se))
            df_print[1:n_match, 2] <- new_se[1:n_match]
            # Recalculate t/z values and p-values
            df_print[1:n_match, 3] <- df_print$Estimate[1:n_match] / pmax(df_print[1:n_match, 2], 1e-12)
          }

          if (is_asymptotic) {
            colnames(df_print) <- c("Estimate", "Std. Error", "z value", "Pr")
            df_print$df <- Inf
            crit <- 1.96
            if (!is.null(self$vcov)) {
               df_print[1:n_match, 4] <- 2 * stats::pnorm(-abs(df_print[1:n_match, 3]))
            }
          } else {
            colnames(df_print) <- c("Estimate", "Std. Error", "t value", "Pr")
            df_print$df <- self$fit$df.residual
            crit <- stats::qt(0.975, df = pmax(df_print$df, 1e-6))
            if (!is.null(self$vcov)) {
               df_print[1:n_match, 4] <- 2 * stats::pt(-abs(df_print[1:n_match, 3]), df = df_print$df[1:n_match])
            }
          }

          df_print$`Lower 95%` <- df_print$Estimate - crit * df_print$`Std. Error`
          df_print$`Upper 95%` <- df_print$Estimate + crit * df_print$`Std. Error`
          rownames(df_print)[rownames(df_print) == "(Intercept)"] <- "Intercept"
        } else {
          df_print <- self$fit
        }

        # --- View Prioritization Logic ---
        view_to_use <- if (!is.null(self$view)) self$view else self$model$view
        if (!is.null(view_to_use)) {
          var_names <- rownames(df_print)
          base_names <- gsub("\\[.*\\]$", "", var_names)
          priority_idx <- integer(0)
          for (v in view_to_use) {
            match_idx <- which(var_names == v | base_names == v)
            priority_idx <- c(priority_idx, match_idx)
          }
          priority_idx <- unique(priority_idx)
          other_idx <- setdiff(seq_len(nrow(df_print)), priority_idx)
          df_print <- df_print[c(priority_idx, other_idx), , drop = FALSE]
        }

        # --- Truncation Logic ---
        if (!is.null(max_rows) && nrow(df_print) > max_rows) {
           res$truncated <- TRUE
           res$total_rows <- nrow(df_print)
           res$max_rows <- max_rows
        } else {
           res$truncated <- FALSE
        }

        res$coefficients <- df_print

        # --- Internal to Requested Contrast Transformation for Display ---
        req_ct <- if (!is.null(self$model$requested_contrasts)) self$model$requested_contrasts else "treatment"
        curr_ct <- if (!is.null(self$model$contrasts)) self$model$contrasts else "sum"

        if (req_ct != curr_ct) {
          # Use names for safe identification of fixed effects
          beta_full <- df_print$Estimate
          names(beta_full) <- rownames(df_print)
          V_full <- self$vcov

          # Identify fixed effects by name pattern (Strict matching for Intercept and b)
          fe_names <- names(beta_full)[grepl("^(Intercept|Intercept_c|b($|\\[))", names(beta_full))]
          # Exclude non-fixed effect parameters
          fe_names <- setdiff(fe_names, c("sigma", "sd", "rho", "phi", "nu"))
          fe_names <- fe_names[!grepl("^(r_re|u|z|lambda|tau|p|prob)($|\\[)", fe_names)]

          if (length(fe_names) > 0) {
            beta <- beta_full[fe_names]

            # Extract corresponding block from V_full by name
            common_fe <- intersect(fe_names, rownames(V_full))
            if (length(common_fe) > 0) {
               V <- V_full[common_fe, common_fe, drop = FALSE]
               beta <- beta[common_fe]
               current_fe_names <- common_fe

               formula <- nobars(self$model$formula)
               predictor_vars <- all.vars(delete.response(terms(formula)))
               relevant_data <- self$model$raw_data[, predictor_vars, drop = FALSE]
               # For minimal grid, use levels for factors and mean for numeric
               levs <- lapply(relevant_data, function(x) if(is.factor(x)) levels(x) else mean(as.numeric(x), na.rm = TRUE))
               grid <- expand.grid(levs)

               # X_from (current internal: usually sum)
               old_opts <- options(contrasts = if (curr_ct == "sum") c("contr.sum", "contr.poly") else c("contr.treatment", "contr.poly"))
               X_from <- model.matrix(delete.response(terms(formula)), grid)
               options(old_opts)

               # X_to (requested: usually treatment)
               old_opts <- options(contrasts = if (req_ct == "treatment") c("contr.treatment", "contr.poly") else c("contr.sum", "contr.poly"))
               X_to <- model.matrix(delete.response(terms(formula)), grid)
               options(old_opts)

               if (ncol(X_from) == length(beta) && ncol(X_to) == length(beta)) {
                 M <- MASS::ginv(X_to) %*% X_from
                 beta_new <- as.numeric(M %*% beta)
                 V_new <- M %*% V %*% t(M)
                 se_new <- sqrt(diag(V_new))

                 # Update df_print with new values for fixed effects
                 df_print$Estimate[current_fe_names] <- beta_new
                 df_print$`Std. Error`[current_fe_names] <- se_new

                 stat_col <- if ("z value" %in% names(df_print)) "z value" else "t value"
                 df_print[[stat_col]][current_fe_names] <- beta_new / pmax(se_new, 1e-12)

                 if (!is.null(df_print$Pr)) {
                   df_print$Pr[current_fe_names] <- 2 * stats::pt(-abs(df_print[[stat_col]][current_fe_names]), df = df_print$df[current_fe_names])
                 }
                 if ("Lower 95%" %in% names(df_print)) {
                   df_print$`Lower 95%`[current_fe_names] <- beta_new - 1.96 * se_new
                   df_print$`Upper 95%`[current_fe_names] <- beta_new + 1.96 * se_new
                 }
                 # Update Row Names to match requested contrast coding
                 new_names <- colnames(X_to)
                 new_names[new_names == "(Intercept)"] <- "Intercept"
                 is_fe <- new_names != "Intercept"
                 new_names[is_fe] <- paste0("b[", new_names[is_fe], "]")

                 rownames(df_print)[which(rownames(df_print) %in% current_fe_names)] <- new_names
               }
            }
          }
        }

        res$coefficients <- df_print
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
        is_asymptotic <- inherits(self$model, "rtmb_loglinear") || (!is.null(self$model$type) && self$model$type == "table")
        if (is_asymptotic) {
          desired_order <- c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "z value", "Pr")
        } else {
          desired_order <- c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "z value", "t value", "df", "Pr")
        }
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

        # Remove df for asymptotic models if requested
        if (is_asymptotic && "df" %in% names(df_final)) {
          df_final$df <- NULL
        }

        # Rename sig to blank just before printing
        if ("sig" %in% names(df_final)) {
           names(df_final)[names(df_final) == "sig"] <- ""
        }

        res$coefficients <- df_final

        # --- Enhanced Output for LM/GLM ---
        if (inherits(self$fit, c("lm", "glm"))) {
          s_lm <- summary(self$fit)
          res$lm_info <- list(
            dispersion = if (inherits(self$fit, "glm")) s_lm$dispersion else NULL,
            null_deviance = if (inherits(self$fit, "glm")) s_lm$null.deviance else NULL,
            df_null = if (inherits(self$fit, "glm")) s_lm$df.null else NULL,
            deviance = if (inherits(self$fit, "glm")) s_lm$deviance else NULL,
            df_residual = s_lm$df.residual,
            sigma = if (inherits(self$fit, "lm")) s_lm$sigma else NULL,
            df_lm = if (inherits(self$fit, "lm")) s_lm$df else NULL,
            r_squared = if (inherits(self$fit, "lm")) s_lm$r.squared else NULL,
            adj_r_squared = if (inherits(self$fit, "lm")) s_lm$adj.r.squared else NULL,
            fstatistic = if (inherits(self$fit, "lm")) s_lm$fstatistic else NULL
          )
        }
      } else {
        res$fit_summary <- self$fit
      }

      class(res) <- "summary_Classic_Fit"
      return(res)
    },


    #' @description Perform ANOVA (Wald F-tests / Chisq-tests) on the fitted model.
    #' @param method Character; "reml" (standard) or "ls" (experimental).
    #' @param type Integer; Type of Sum of Squares (only Type III supported currently).
    #' @return A data frame containing the ANOVA table.
    anova = function(method = c("reml", "ls"), type = 3) {
      if (!is.null(self$model$type) && self$model$type == "table") {
        obs_tab <- self$model$extra$tab
        N_tot <- sum(obs_tab)
        R_dim <- nrow(obs_tab)
        C_dim <- ncol(obs_tab)

        df_mu <- if (is.data.frame(self$fit)) self$fit[grepl("^mu\\[", rownames(self$fit)), , drop = FALSE] else data.frame()

        if (nrow(df_mu) == R_dim * C_dim) {
          p_est <- df_mu$Estimate / N_tot
          p_mat <- matrix(p_est, nrow = R_dim, ncol = C_dim)
          p_row <- rowSums(p_mat)
          p_col <- colSums(p_mat)
          E_mat <- matrix(0, nrow = R_dim, ncol = C_dim)
          for (i in 1:R_dim) for (j in 1:C_dim) E_mat[i, j] <- p_row[i] * p_col[j] * N_tot
          chisq_stat <- sum((as.vector(obs_tab) - as.vector(E_mat))^2 / as.vector(E_mat))
          df_val <- (R_dim - 1) * (C_dim - 1)
          p_val <- stats::pchisq(chisq_stat, df = df_val, lower.tail = FALSE)

          chisq_res <- data.frame(`X-squared` = chisq_stat, df = df_val, `p-value` = p_val, check.names = FALSE)
          rownames(chisq_res) <- "Pearson's Chi-squared (from est. prob)"
        } else {
          chisq_res <- NULL
        }

        fisher_res <- if (!is.null(self$test_results$fisher)) self$test_results$fisher else self$model$extra$fisher
        res <- list(chisq = chisq_res, fisher = fisher_res)
        class(res) <- "anova_rtmb_table"
        return(res)
      }

      # --- 2. Special Case: Correlation Analysis (Overall Independence Test) ---
      if (!is.null(self$model$type) && self$model$type == "corr") {
        corr_patterns <- c("corr", "rho", "W_corr", "B_corr", "W_pcorr", "B_pcorr")
        idx_corr <- which(sapply(rownames(self$fit), function(n) {
          any(sapply(corr_patterns, function(p) grepl(paste0("^", p, "($|\\[)"), n)))
        }))

        if (length(idx_corr) > 0) {
          base_names <- gsub("\\[.*\\]$", "", rownames(self$fit)[idx_corr])
          unique_bases <- unique(base_names)

          res_list <- list()
          for (b in unique_bases) {
            sub_idx <- idx_corr[base_names == b]
            rho <- self$fit$Estimate[sub_idx]
            V_rho <- self$vcov[sub_idx, sub_idx, drop = FALSE]

            # Wald test: H0: all rho = 0
            chisq_stat <- as.numeric(t(rho) %*% MASS::ginv(V_rho) %*% rho)
            df_val <- length(sub_idx)
            p_val <- stats::pchisq(chisq_stat, df = df_val, lower.tail = FALSE)

            res_list[[b]] <- data.frame(
              Chisq = chisq_stat,
              Df = df_val,
              `Pr(>Chisq)` = p_val,
              check.names = FALSE
            )
          }
          res <- do.call(rbind, res_list)
          rownames(res) <- paste0("Overall test for ", unique_bases, " (H0: all zero)")
          class(res) <- c("anova", "data.frame")
          return(res)
        }
      }

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

      full_names <- if (is.data.frame(self$fit)) rownames(self$fit) else names(beta_full)
      if (inherits(self$fit, "lm")) {
        full_names[full_names == "(Intercept)"] <- "Intercept"
      }
      names(beta_full) <- full_names

      assign_idx <- self$model$extra$X_assign

      if (inherits(self$fit, "lm")) {
        fe_idx <- seq_along(beta_full)
      } else {
        fe_idx <- which(grepl("^(Intercept|Intercept_c|b($|\\[))", names(beta_full)))
      }

      if (length(fe_idx) == 0) stop("Could not identify fixed effects for ANOVA.")

      beta <- beta_full[fe_idx]
      V <- V_full[fe_idx, fe_idx]

      full_assign <- c(0, assign_idx)

      ct_setting <- if (!is.null(self$model$contrasts)) self$model$contrasts else "sum"
      if (ct_setting != "sum") {
        formula <- nobars(self$model$formula)
        predictor_vars <- all.vars(delete.response(terms(formula)))
        relevant_data <- self$model$raw_data[, predictor_vars, drop = FALSE]
        if (is.null(relevant_data)) relevant_data <- self$model$data[, predictor_vars, drop = FALSE]
        levs <- lapply(relevant_data, function(x) if(is.factor(x) || is.character(x)) unique(x) else mean(x, na.rm=TRUE))
        grid <- expand.grid(levs)

        old_opts2 <- options(contrasts = if (ct_setting == "treatment")
          c("contr.treatment", "contr.poly") else options()$contrasts)
        X_curr <- model.matrix(delete.response(terms(formula)), grid)
        options(old_opts2)

        old_opts3 <- options(contrasts = c("contr.sum", "contr.poly"))
        X_sum <- model.matrix(delete.response(terms(formula)), grid)
        options(old_opts3)

        if (ncol(X_curr) == length(beta)) {
          M <- MASS::ginv(X_sum) %*% X_curr
          beta <- as.numeric(M %*% beta)
          names(beta) <- colnames(X_sum)
          V <- M %*% V %*% t(M)
          assign_idx <- attr(X_sum, "assign")
        }
      }

      terms <- self$model$extra$X_terms
      res_list <- list()

      has_int <- any(grepl("^Intercept", names(beta))) || any(full_assign == 0)
      all_assign <- seq_along(terms)
      all_term_names <- terms

      for (i in seq_along(all_assign)) {
        a_id <- all_assign[i]
        t_name <- all_term_names[i]

        idx <- which(assign_idx == a_id)
        if (length(idx) == 0) next

        L <- matrix(0, nrow = length(idx), ncol = length(beta))
        for (j in seq_along(idx)) L[j, idx[j]] <- 1

        LVL <- L %*% V %*% t(L)
        inv_LVL <- try(solve(LVL), silent = TRUE)
        if (inherits(inv_LVL, "try-error")) inv_LVL <- MASS::ginv(LVL)

        W <- as.numeric(t(L %*% beta) %*% inv_LVL %*% (L %*% beta))
        df1 <- length(idx)
        f_val <- W / df1

        if (is.data.frame(self$fit) && !is.null(self$fit$df)) {
          df2 <- min(self$fit$df[fe_idx[idx]], na.rm = TRUE)
        } else if (inherits(self$fit, "lm")) {
          df2 <- self$fit$df.residual
        } else {
          df2 <- Inf
        }

        # --- 追加: 漸近的なモデルの場合は強制的に Chi-squared 検定にする ---
        is_asymp <- (!is.null(self$model$family) && self$model$family %in% c("poisson", "binomial", "bernoulli", "multinomial", "ordered")) || (!is.null(self$model$type) && self$model$type %in% c("loglinear", "table"))
        if (is_asymp) df2 <- Inf

        if (df2 == 0 || is.infinite(df2)) {
          p_val <- stats::pchisq(W, df1, lower.tail = FALSE)
          res_list[[t_name]] <- data.frame(`Df` = df1, `Chisq` = W, `Pr(>Chisq)` = p_val, row.names = t_name, check.names = FALSE)
        } else {
          p_val <- stats::pf(f_val, df1, df2, lower.tail = FALSE)
          res_list[[t_name]] <- data.frame(`NumDF` = df1, `DenDF` = df2, `F value` = f_val, `Pr(>F)` = p_val, row.names = t_name, check.names = FALSE)
        }
      }

      res_df <- do.call(rbind, res_list)
      p_col_name <- if ("Pr(>Chisq)" %in% names(res_df)) "Pr(>Chisq)" else "Pr(>F)"
      sig <- symnum(res_df[[p_col_name]], corr = FALSE, na = FALSE,
                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                    symbols = c("***", "**", "*", ".", " "))
      res_df$signif <- as.character(sig)

      class(res_df) <- c("anova", "data.frame")
      heading <- if ("Chisq" %in% names(res_df)) "ANOVA Table (Wald Chisq tests)" else "ANOVA Table (Wald F-tests)"
      attr(res_df, "heading") <- heading
      return(res_df)
    },

    #' @description Calculate Least Squares Means (Marginal Means) and contrasts.
    #' @param specs Character vector of factors to calculate means for.
    #' @param pairwise Logical; whether to perform pairwise comparisons.
    #' @param simple Character vector of factors to hold constant for simple main effects.
    #' @param adjust Character; p-value adjustment method (e.g., "bonferroni", "holm", "none").
    #' @return A data frame containing the marginal means or contrasts.
    lsmeans = function(specs, pairwise = FALSE, simple = NULL, adjust = "bonferroni") {
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

      if (!is.null(self$model$raw_data)) {
        data <- as.data.frame(self$model$raw_data)
      } else if (!is.null(self$model$data)) {
        data <- as.data.frame(self$model$data)
      } else {
        stop("lsmeans requires the original data.")
      }
      formula <- self$model$formula

      # 1. Identify all factors in the model
      vars <- all.vars(nobars(formula))[-1]
      is_cat <- sapply(data[vars], function(x) is.factor(x) || is.character(x))
      cat_vars <- vars[is_cat]

      full_specs <- unique(c(specs, simple))
      if (!all(full_specs %in% cat_vars)) stop("specs/simple must be categorical factors in the model.")

      # 2. Create a reference grid for all categorical factors
      grid_list <- lapply(data[cat_vars], function(x) levels(as.factor(x)))
      ref_grid <- expand.grid(grid_list)

      # 3. Handle continuous covariates by setting them to their mean
      cont_vars <- vars[!is_cat]
      for (v in cont_vars) {
        ref_grid[[v]] <- mean(data[[v]], na.rm = TRUE)
      }

      # 4. Generate model matrix for the reference grid
      mf_orig <- model.frame(nobars(formula), data)
      orig_terms <- delete.response(terms(mf_orig))

      ct_setting <- if (!is.null(self$model$obj$contrasts)) self$model$obj$contrasts else "sum"
      ct_list <- list()
      for (v in cat_vars) {
        ct_list[[v]] <- if (ct_setting == "treatment") "contr.treatment" else "contr.sum"
      }

      X_grid_raw <- model.matrix(orig_terms, ref_grid, contrasts.arg = ct_list)

      # Match columns with the original model
      orig_cols <- self$model$extra$X_colnames
      if (!is.null(orig_cols)) {
        grid_cols <- colnames(X_grid_raw)
        grid_cols[grid_cols == "(Intercept)"] <- "Intercept"
        col_idx <- match(orig_cols, colnames(X_grid_raw))
        if (any(is.na(col_idx)) && "Intercept" %in% orig_cols) {
           col_idx[orig_cols == "Intercept"] <- which(colnames(X_grid_raw) == "(Intercept)")
        }
        if (any(is.na(col_idx))) stop("Could not match lsmeans grid columns with original model parameters.")
        X_grid <- X_grid_raw[, col_idx, drop = FALSE]
      } else {
        X_grid <- X_grid_raw
      }

      beta_vals <- if (is.data.frame(self$fit)) self$fit$Estimate else stats::coef(self$fit)
      fe_idx <- which(grepl("^(Intercept|Intercept_c|b\\[)", rownames(self$fit)))
      beta_match <- beta_vals[fe_idx]
      V_match <- self$vcov[fe_idx, fe_idx]

      # Group by specs + simple with descriptive labels
      label_grid <- ref_grid[full_specs]
      for (v in full_specs) {
        label_grid[[v]] <- paste0(v, "=", label_grid[[v]])
      }
      groups_full <- interaction(label_grid, drop = TRUE, sep = ":")
      unique_groups_full <- levels(groups_full)

      # Calculate L-vectors for each group
      L_list <- list()
      for (grp in unique_groups_full) {
        idx <- which(groups_full == grp)
        L_list[[grp]] <- colMeans(X_grid[idx, , drop = FALSE])
      }

      if (!pairwise) {
        # --- Marginal Means ---
        res_list <- list()
        for (grp in unique_groups_full) {
          L <- L_list[[grp]]
          est <- as.numeric(L %*% beta_match)
          se <- sqrt(as.numeric(t(L) %*% V_match %*% L))
          df_val <- self$.get_lsmeans_df(specs)

          res_list[[grp]] <- data.frame(
            estimate = est, `Std. Error` = se, df = df_val,
            `Lower 95%` = est + qt(0.025, df_val) * se,
            `Upper 95%` = est + qt(0.975, df_val) * se,
            row.names = grp, check.names = FALSE
          )
        }
        res <- do.call(rbind, res_list)
        class(res) <- c("rtmb_lsmeans", class(res))
        return(res)
      } else {
        # --- Pairwise Comparisons ---
        contrast_list <- list()
        # If simple is provided, we only compare within levels of 'simple'
        if (!is.null(simple)) {
          mod_grid <- unique(ref_grid[simple])
          focal_grid <- unique(ref_grid[specs])

          for (row_m in 1:nrow(mod_grid)) {
            m_vals <- mod_grid[row_m, , drop = FALSE]
            for (i in 1:(nrow(focal_grid)-1)) {
              for (j in (i+1):nrow(focal_grid)) {
                # Find indices in unique_groups_full that match focal + moderator
                # This is more robust than string pasting
                find_idx <- function(f_row, m_row) {
                  vals <- character(length(full_specs))
                  vals[1:length(specs)] <- as.character(unlist(f_row))
                  vals[(length(specs)+1):length(full_specs)] <- as.character(unlist(m_row))
                  target_label <- paste(vals, collapse = ":")
                  idx <- which(unique_groups_full == target_label)
                  if (length(idx) == 0) return(NULL)
                  return(idx)
                }

                idx_i <- find_idx(focal_grid[i,], m_vals)
                idx_j <- find_idx(focal_grid[j,], m_vals)

                if (!is.null(idx_i) && !is.null(idx_j)) {
                  name_i <- unique_groups_full[idx_i]
                  name_j <- unique_groups_full[idx_j]
                  L_diff <- L_list[[name_i]] - L_list[[name_j]]
                  contrast_label <- paste0("(", name_i, ") - (", name_j, ")")
                  contrast_list[[contrast_label]] <- self$.calc_contrast(L_diff, specs)
                }
              }
            }
          }
        } else {
          # All pairwise
          for (i in 1:(length(unique_groups_full)-1)) {
            for (j in (i+1):length(unique_groups_full)) {
              name_i <- unique_groups_full[i]
              name_j <- unique_groups_full[j]
              L_diff <- L_list[[name_i]] - L_list[[name_j]]
              contrast_label <- paste0("(", name_i, ") - (", name_j, ")")
              contrast_list[[contrast_label]] <- self$.calc_contrast(L_diff, specs)
            }
          }
        }

        res <- do.call(rbind, contrast_list)
        if (adjust != "none") {
          res$Pr <- p.adjust(res$Pr, method = adjust)
          attr(res, "adjustment") <- adjust
        }
        class(res) <- c("rtmb_contrasts", class(res))
        return(res)
      }
    },

    #' @description (Internal) Calculate metrics for a contrast.
    #' @param L Contrast matrix.
    #' @param specs Variable names for lsmeans.
    #' @return A data frame with estimate, SE, df, t-value, and p-value.
    .calc_contrast = function(L, specs) {
      beta_vals <- if (is.data.frame(self$fit)) self$fit$Estimate else stats::coef(self$fit)
      fe_idx <- which(grepl("^(Intercept|Intercept_c|b\\[)", rownames(self$fit)))
      beta_match <- beta_vals[fe_idx]
      V_match <- self$vcov[fe_idx, fe_idx]

      est <- as.numeric(L %*% beta_match)
      se <- sqrt(as.numeric(t(L) %*% V_match %*% L))
      df_val <- self$.get_lsmeans_df(specs)
      t_val <- est / pmax(se, 1e-12)
      p_val <- 2 * pt(-abs(t_val), df = df_val)

      data.frame(
        estimate = est, `Std. Error` = se, df = df_val,
        `t value` = t_val, Pr = p_val,
        check.names = FALSE
      )
    },

    #' @description (Internal) Get representative DF for lsmeans.
    #' @param specs Variable names for lsmeans.
    #' @return Degrees of freedom.
    .get_lsmeans_df = function(specs) {
      df_val <- Inf
      if (is.data.frame(self$fit) && !is.null(self$fit$df)) {
        match_pattern <- paste0("^b\\[(", paste(specs, collapse="|"), ")")
        match_idx <- grepl(match_pattern, rownames(self$fit))
        if (any(match_idx)) {
          df_val <- min(self$fit$df[match_idx], na.rm = TRUE)
        } else if (inherits(self$fit, "lm")) {
          df_val <- self$fit$df.residual
        }
      } else if (inherits(self$fit, "lm")) {
        df_val <- self$fit$df.residual
      }
      return(df_val)
    },

    #' @description (Internal) Construct a list of parameters from the fit.
    #' @param fit The fit result (dataframe or lm object).
    #' @return A named list of parameters.
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
        if (!is.null(self$model$type) && self$model$type == "corr") {
          par_list$corr <- est[grepl("^corr\\[", names(est))]
        }
        par_list$IE <- est[grepl("^IE_", names(est))]
      }
      return(par_list)
    }
  )
)

#' @export
summary.Classic_Fit <- function(object, ...) {
  object$summary(...)
}

#' @export
#' @export
print.summary_Classic_Fit <- function(x, ...) {
  digits <- if (!is.null(x$digits)) x$digits else 5

  if (!is.null(x$type) && x$type == "table") {
    cat("\nContingency Table Analysis\n")
  } else if (!is.null(x$type) && x$type == "loglinear") {
    cat("\nLog-Linear Model Analysis\n")
  } else {
    cat("\nCall:\n")
    type_label <- if (is.null(x$type)) "Generic Model" else x$type
    cat(paste("Classical estimation via", type_label), "\n")
  }

  if (!is.null(x$se_method) && x$se_method != "wald") {
    method_label <- switch(x$se_method,
                           robust = if (!is.null(x$cluster)) paste("Robust (Cluster:", x$cluster, ")") else "Robust (HC3)",
                           bootstrap = paste0("Bootstrap (", x$bootstrap, " samples)"),
                           "Standard")
    cat(sprintf("Standard Errors: %s\n", method_label))
  }

  if (!is.null(x$logLik) && !is.na(x$logLik)) {
    cat(sprintf("\nLog-Likelihood: %.3f, AIC: %.3f, BIC: %.3f\n",
                as.numeric(x$logLik), x$AIC, x$BIC))
  }

  if (!is.null(x$type) && x$type == "ttest") {
    levs <- x$extra$levs
    cat(sprintf("(Comparison: %s - %s)\n", levs[1], levs[2]))
  } else if (!is.null(x$type) && x$type == "table" && !is.null(x$extra$tab)) {
    cat("\nContingency Table (Observed):\n")
    print(x$extra$tab)
  }

  if (!is.null(x$coefficients)) {
    # ---------------------------------------------------------
    # 分岐 1: TABLE 出力
    # ---------------------------------------------------------
    if (!is.null(x$type) && x$type == "table") {
      cat("\n---\n")
      obs_tab <- x$extra$tab
      N_tot <- sum(obs_tab)
      R_dim <- nrow(obs_tab)
      C_dim <- ncol(obs_tab)

      df_mu <- x$coefficients[grepl("^mu\\[", rownames(x$coefficients)), , drop = FALSE]

      if (nrow(df_mu) == R_dim * C_dim) {
        mu_est <- df_mu$Estimate
        p_est <- mu_est / N_tot

        p_mat <- matrix(p_est, nrow = R_dim, ncol = C_dim)
        p_row <- rowSums(p_mat)
        p_col <- colSums(p_mat)

        E_mat <- matrix(0, nrow = R_dim, ncol = C_dim)
        for (i in 1:R_dim) for (j in 1:C_dim) E_mat[i, j] <- p_row[i] * p_col[j] * N_tot

        E_vec <- as.vector(E_mat)
        obs_vec <- as.vector(obs_tab)

        cat("Cell Probabilities (p) and Confidence Intervals:\n")
        df_p <- df_mu
        df_p$Estimate <- df_p$Estimate / N_tot
        df_p$`Std. Error` <- df_p$`Std. Error` / N_tot
        if ("Lower 95%" %in% names(df_p)) {
          df_p$`Lower 95%` <- df_p$`Lower 95%` / N_tot
          df_p$`Upper 95%` <- df_p$`Upper 95%` / N_tot
        }
        rownames(df_p) <- gsub("^mu\\[", "p[", rownames(df_p))
        cols_to_show <- intersect(c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "z value", "t value", "Pr"), names(df_p))
        print(df_p[, cols_to_show, drop = FALSE], quote = FALSE, right = TRUE)

        cat("\nExpected Counts (Independence) and Pearson Residuals:\n")
        residuals <- (obs_vec - E_vec) / pmax(sqrt(E_vec), 1e-8)
        rnames <- gsub("^mu\\[", "E\\[", rownames(df_mu))

        df_resid <- data.frame(
          Expected = round(E_vec, digits),
          Residual = round(residuals, digits),
          row.names = rnames
        )
        print(df_resid, quote = FALSE, right = TRUE)
      }

      # ---------------------------------------------------------
      # 分岐 2: LOGLINEAR 出力
      # ---------------------------------------------------------
    } else if (!is.null(x$type) && x$type == "loglinear") {
      cat("\nLog-linear Parameters and Confidence Intervals:\n")
      df_params <- x$coefficients[!grepl("^mu\\[", rownames(x$coefficients)), , drop = FALSE]
      cols_to_show <- intersect(c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "z value", "t value", "Pr"), names(df_params))
      print(df_params[, cols_to_show, drop = FALSE], quote = FALSE, right = TRUE)

      df_mu <- x$coefficients[grepl("^mu\\[", rownames(x$coefficients)), , drop = FALSE]
      if (nrow(df_mu) > 0 && !is.null(x$extra$obs_Y)) {
        cat("\n---\nExpected Counts (Fitted) and Pearson Residuals:\n")
        E_vec <- df_mu$Estimate
        obs_vec <- x$extra$obs_Y
        residuals <- (obs_vec - E_vec) / pmax(sqrt(E_vec), 1e-8)

        df_resid <- data.frame(
          Observed = obs_vec,
          Expected = round(E_vec, digits),
          Residual = round(residuals, digits),
          row.names = rownames(df_mu)
        )

        # 行数が多い場合は省略表示する
        max_r <- if (!is.null(x$max_rows)) x$max_rows else 20
        if (nrow(df_resid) > max_r) {
          print(df_resid[1:max_r, , drop=FALSE], quote=FALSE, right=TRUE)
          cat(sprintf("... (omitted %d rows)\n", nrow(df_resid) - max_r))
        } else {
          print(df_resid, quote = FALSE, right = TRUE)
        }
      }

      # ---------------------------------------------------------
      # 分岐 3: 通常の出力 (lm, glm, lmer, corr 等)
      # ---------------------------------------------------------
    } else {
      cat("\nPoint Estimates and Confidence Intervals:\n")
      df_to_print <- x$coefficients
      if (isTRUE(x$truncated)) df_to_print <- df_to_print[1:x$max_rows, , drop = FALSE]
      print(df_to_print, quote = FALSE, right = TRUE)
      if (isTRUE(x$truncated)) {
        cat(sprintf("... (omitted %d parameters; use summary(max_rows = ...) to show more)\n",
                    x$total_rows - x$max_rows))
      }
    }

    if (!is.null(x$lm_info)) {
      cat("\n---\n")
      info <- x$lm_info
      if (!is.null(info$dispersion)) {
        cat(sprintf("Dispersion parameter for %s family taken to be %s\n",
                    x$family, format(round(info$dispersion, digits), nsmall = digits)))
        cat(sprintf("Null deviance: %s on %d degrees of freedom\n",
                    format(round(info$null_deviance, 2), nsmall = 2), info$df_null))
        cat(sprintf("Residual deviance: %s on %d degrees of freedom\n",
                    format(round(info$deviance, 2), nsmall = 2), info$df_residual))
      } else if (!is.null(info$sigma)) {
        cat(sprintf("Residual standard error: %s on %d degrees of freedom\n",
                    format(round(info$sigma, digits), nsmall = digits), info$df_residual))
        cat(sprintf("Multiple R-squared: %s, Adjusted R-squared: %s\n",
                    format(round(info$r_squared, 4), nsmall = 4),
                    format(round(info$adj_r_squared, 4), nsmall = 4)))
        if (!is.null(info$fstatistic)) {
          f <- info$fstatistic
          p_f <- stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
          cat(sprintf("F-statistic: %s on %d and %d DF, p-value: %s\n",
                      format(round(f[1], 2), nsmall = 2), f[2], f[3],
                      if (p_f < 0.001) "< .001" else format(round(p_f, 4), nsmall = 4)))
        }
      }
    }
  } else if (!is.null(x$fit_summary)) {
    print(x$fit_summary)
  }

  invisible(x)
}

#' @export
print.anova_rtmb_table <- function(x, ...) {
  cat("\nContingency Table Analysis Tests\n")
  cat(rep("-", 40), "\n", sep = "")

  if (!is.null(x$chisq)) {
    print(x$chisq, row.names = TRUE, right = TRUE)
  }

  if (!is.null(x$fisher) && !inherits(x$fisher, "try-error")) {
    cat("\nFisher's Exact Test for Count Data\n")
    cat(sprintf("p-value = %.5f\n", x$fisher$p.value))
    if (!is.null(x$fisher$estimate)) {
      cat(sprintf("alternative hypothesis: %s\n", x$fisher$alternative))
      cat(sprintf("odds ratio: %.4f\n", x$fisher$estimate))
    }
  }
  invisible(x)
}

#' @export
anova.Classic_Fit <- function(object, ...) {
  object$anova(...)
}

#' @export
print.Classic_Fit <- function(x, ...) {
  x$print(...)
}

#' Least-squares means (marginal means)
#'
#' Generic function to calculate least-squares means (marginal means) for fitted models.
#'
#' @param object A fitted model object.
#' @param specs Variable names to calculate marginal means for.
#' @param ... Additional arguments.
#' @export
lsmeans <- function(object, specs, ...) {
  UseMethod("lsmeans")
}

#' @export
lsmeans.Classic_Fit <- function(object, specs, ...) {
  object$lsmeans(specs, ...)
}

#' @export
logLik.Classic_Fit <- function(object, ...) {
  object$logLik()
}

#' @export
AIC.Classic_Fit <- function(object, ..., k = 2) {
  object$AIC()
}

#' @export
BIC.Classic_Fit <- function(object, ...) {
  object$BIC()
}
