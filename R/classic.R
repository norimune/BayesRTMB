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
#' @field robust_type Character string specifying the robust standard error correction type, if used.
#' @field bootstrap_results A matrix containing bootstrap samples, if applicable.
#' @field test_results List of additional test results (e.g., chisq.test).
#' @field view Character vector of parameter names to prioritize in summary.
#' @field par_vec Numeric vector of parameter estimates.
#' @field objective Final objective value.
#' @field log_lik Log-likelihood.
#' @field restricted_log_lik Restricted log-likelihood, if classical estimation used REML-style marginalization.
#' @field rss Residual sum of squares for classical linear models, if available.
#' @field df_residual Residual degrees of freedom for classical linear models, if available.
#' @field convergence Convergence code.
#' @field sd_rep TMB sdreport object.
#' @field df_fixed Dataframe of fixed effects results.
#' @field random_effects Dataframe of random effects results.
#' @field df_transform Dataframe of transformed parameters.
#' @field df_generate Dataframe of generated quantities.
#' @field opt_history Dataframe of optimization history.
#' @field transform List of transformed parameters.
#' @field generate List of generated quantities.
#' @field se_samples List of simulated samples for SE estimation.
#' @field par_unc Numeric vector of unconstrained parameter estimates.
#' @field vcov_unc Variance-covariance matrix of parameters in unconstrained space.
#' @field ci_method Method used for CI estimation.
#' @field laplace Whether Laplace approximation was used.
#' @field map Parameter mapping used.
#' @field df_method Character string specifying the degrees of freedom calculation method.
#' @field idx_fix_active Numeric vector; mapping between active parameters and full unconstrained vector.
#' @field show_df Logical; whether to display degrees of freedom in the summary output.
#' @importFrom stats AIC BIC logLik anova
#' @export
Classic_Fit <- R6::R6Class(
  classname = "Classic_Fit",
  inherit = RTMB_Fit_Base,
  public = list(
    model = NULL,
    fit = NULL,
    par = NULL,
    par_vec = NULL,
    vcov = NULL,

    #' @description Get point estimate for a target parameter.
    #' @param target Target parameter name.
    #' @param ... Additional arguments, ignored for classic fits.
    #' @return Matrix, array, vector, or scalar point estimate.
    get_point_estimate = function(target, ...) {
      if (!is.null(self$par[[target]])) return(self$par[[target]])
      if (!is.null(self$transform[[target]])) return(self$transform[[target]])
      if (!is.null(self$generate[[target]])) return(self$generate[[target]])
      stop("Parameter not found: ", target)
    },
    objective = NULL,
    log_lik = NULL,
    restricted_log_lik = NULL,
    rss = NULL,
    df_residual = NULL,
    convergence = NULL,
    sd_rep = NULL,
    df_fixed = NULL,
    random_effects = NULL,
    df_transform = NULL,
    df_generate = NULL,
    opt_history = NULL,
    transform = NULL,
    generate = NULL,
    se_samples = NULL,
    par_unc = NULL,
    vcov_unc = NULL,
    bootstrap_results = NULL,
    test_results = list(),
    se_method = "wald",
    ci_method = "wald",
    laplace = TRUE,
    map = NULL,
    cluster = NULL,
    view = NULL,
    df_method = "auto",
    idx_fix_active = NULL,
    show_df = TRUE,
    robust_type = NULL,

    #' @description Create a new `Classic_Fit` object.
    #' @param model The `RTMB_Model` object.
    #' @param par_vec Numeric vector of parameter estimates.
    #' @param par List of parameter estimates.
    #' @param objective Final objective value.
    #' @param log_lik Full log-likelihood used for information criteria.
    #' @param restricted_log_lik Restricted log-likelihood, if classical estimation used REML-style marginalization.
    #' @param rss Residual sum of squares for classical linear models, if available.
    #' @param df_residual Residual degrees of freedom for classical linear models, if available.
    #' @param convergence Convergence code.
    #' @param sd_rep TMB sdreport object.
    #' @param df_fixed Dataframe of fixed effects results.
    #' @param random_effects Dataframe of random effects results.
    #' @param df_transform Dataframe of transformed parameters.
    #' @param df_generate Dataframe of generated quantities.
    #' @param opt_history Dataframe of optimization history.
    #' @param transform List of transformed parameters.
    #' @param generate List of generated quantities.
    #' @param se_samples List of simulated samples for SE estimation.
    #' @param par_unc Parameter vector on unconstrained scale.
    #' @param vcov_unc Covariance matrix on unconstrained scale.
    #' @param ci_method Method used for CI estimation.
    #' @param laplace Whether Laplace approximation was used.
    #' @param map Parameter mapping used.
    #' @param test_results List of additional test results (e.g., chisq.test).
    #' @param view Character vector of parameter names to prioritize in summary.
    #' @param fit Legacy argument for backward compatibility (maps to df_fixed).
    #' @param vcov Variance-covariance matrix of parameters.
    #' @param df_method Method for degrees of freedom calculation.
    #' @param idx_fix_active Numeric vector; mapping between active parameters and full unconstrained vector.
    #' @param show_df Logical; whether to display degrees of freedom in the summary output.
    #' @param ... Additional arguments passed to the constructor.
    initialize = function(model, par_vec = NULL, par = NULL, objective = NULL, log_lik = NULL,
                          restricted_log_lik = NULL,
                          convergence = NULL, sd_rep = NULL, df_fixed = NULL, random_effects = NULL,
                          df_transform = NULL, df_generate = NULL, opt_history = NULL,
                          transform = NULL, generate = NULL, se_samples = NULL, par_unc = NULL,
                          vcov_unc = NULL, ci_method = "wald", laplace = TRUE, map = NULL,
                          test_results = list(), view = NULL, fit = NULL, vcov = NULL, 
                          df_method = "auto", idx_fix_active = NULL, show_df = TRUE,
                          rss = NULL, df_residual = NULL, ...) {
      self$model <- model
      self$par_vec <- par_vec
      self$par <- if (!is.null(par)) par else if (!is.null(df_fixed)) self$.construct_par_list(df_fixed) else if (!is.null(fit)) self$.construct_par_list(fit) else NULL
      self$objective <- objective
      
      # Handle log_ml fallback from ...
      dots <- list(...)
      if (is.null(log_lik) && !is.null(dots$log_ml)) {
        log_lik <- dots$log_ml
      }
      self$log_lik <- log_lik
      self$restricted_log_lik <- restricted_log_lik
      self$rss <- rss
      self$df_residual <- df_residual
      
      self$convergence <- convergence
      self$sd_rep <- sd_rep
      self$df_fixed <- if (!is.null(df_fixed)) df_fixed else fit
      self$fit <- self$df_fixed # Legacy compatibility
      self$random_effects <- random_effects
      self$df_transform <- df_transform
      self$df_generate <- df_generate
      self$opt_history <- opt_history
      self$transform <- transform
      self$generate <- generate
      self$se_samples <- se_samples
      self$par_unc <- par_unc
      self$vcov_unc <- vcov_unc
      self$ci_method <- ci_method
      self$se_method <- ci_method # Legacy compatibility
      self$laplace <- laplace
      self$map <- map
      self$test_results <- test_results
      self$view <- view
      self$vcov <- vcov
      self$df_method <- df_method
      self$idx_fix_active <- idx_fix_active
      self$show_df <- show_df
    },

    #' @description Compute robust standard errors (sandwich estimator).
    #' @param cluster Character; variable name for clustering.
    #' @param type Character; "HC3" (default), "HC0", or "HC1" for non-clustered
    #'   robust SE; "CR1" or "CR0" for clustered robust SE.
    #' @param inplace Logical; if FALSE, return a robust-SE copy without modifying the original fit.
    #'   If TRUE, update the object in place.
    #' @param ... Additional arguments.
    #' @return Self.
    robust_se = function(cluster = NULL, type = c("HC3", "HC0", "HC1", "CR1", "CR0"), inplace = FALSE, ...) {
      type <- match.arg(type)
      cluster_context <- !is.null(cluster) || !is.null(self$model$extra$glmer_info)
      if (isTRUE(cluster_context) && type %in% c("HC3", "HC0", "HC1")) type <- "CR1"
      if (!isTRUE(cluster_context) && type %in% c("CR1", "CR0")) type <- "HC3"

      if (!isTRUE(inplace)) {
        out <- self$clone(deep = TRUE)
        return(out$robust_se(cluster = cluster, type = type, inplace = TRUE, ...))
      }

      # --- 1. Detect if GLMM ---
      if (!is.null(self$model$extra$glmer_info)) {
        return(private$.compute_glmer_robust_se(cluster = cluster, type = type, update = TRUE))
      }

      # --- 2. Fallback to OLS-based sandwich for LM ---
      if (!is.null(self$model$type) && self$model$type == "lm") {
        return(private$.compute_lm_robust_se(cluster = cluster, type = type, update = TRUE))
      }

      # --- 3. Score-based sandwich for GLM ---
      if (!is.null(self$model$type) && self$model$type == "glm") {
        return(private$.compute_glm_robust_se(cluster = cluster, type = type, update = TRUE))
      }

      stop("robust_se() is currently supported for 'lm', 'glm', and mixed models ('lmer'/'glmer').", call. = FALSE)
    },

    #' @description (Deprecated) Use robust_se() instead.
    #' @param ... Arguments passed to robust_se(), including inplace.
    #' @return Self.
    compute_robust = function(...) {
      warning("compute_robust() is deprecated. Use robust_se() instead.", call. = FALSE)
      self$robust_se(...)
    },

    #' @description Compute nonparametric bootstrap standard errors and confidence intervals.
    #' Currently implemented only for mediation models.
    #' Bootstrap refits mediation equations with base R lm.fit()/glm.fit() when possible,
    #' and falls back to RTMB refits for unsupported families.
    #' @param n_boot Integer; number of bootstrap samples.
    #' @param seed Optional integer random seed.
    #' @param inplace Logical; if FALSE, return a bootstrap-SE copy without modifying the original fit.
    #'   If TRUE, update the object in place.
    #' @param ... Additional arguments.
    #' @return A Classic_Fit object.
    bootstrap = function(n_boot = 1000, seed = NULL, inplace = FALSE, ...) {
      self$compute_bootstrap(
        n_boot = n_boot,
        seed = seed,
        inplace = inplace,
        ...
      )
    },

    #' @description Compute nonparametric bootstrap standard errors and confidence intervals.
    #' Currently implemented only for mediation models.
    #' Bootstrap refits mediation equations with base R lm.fit()/glm.fit() when possible,
    #' and falls back to RTMB refits for unsupported families.
    #' @param n_boot Integer; number of bootstrap samples.
    #' @param seed Optional integer random seed.
    #' @param inplace Logical; if FALSE, return a bootstrap-SE copy without modifying the original fit.
    #'   If TRUE, update the object in place.
    #' @param ... Additional arguments.
    #' @return A Classic_Fit object.
    compute_bootstrap = function(n_boot = 1000, seed = NULL, inplace = FALSE, ...) {
      if (!identical(self$model$type, "mediation")) {
        stop(
          "bootstrap() is currently implemented only for mediation models. ",
          "Use robust_se() for lm/glmer robust standard errors.",
          call. = FALSE
        )
      }

      if (!isTRUE(inplace)) {
        out <- self$clone(deep = TRUE)
        return(out$compute_bootstrap(
          n_boot = n_boot,
          seed = seed,
          inplace = TRUE,
          ...
        ))
      }

      if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) || n_boot < 1) {
        stop("n_boot must be a positive integer.", call. = FALSE)
      }

      n_boot <- as.integer(n_boot)

      extract_estimates <- function(fit) {
        dfs <- list(
          fixed = fit$df_fixed,
          random = fit$random_effects,
          transform = fit$df_transform,
          generate = fit$df_generate
        )
        keep <- vapply(dfs, is.data.frame, logical(1))
        dfs <- dfs[keep]
        if (length(dfs) == 0L) return(stats::setNames(numeric(0), character(0)))
        out_df <- private$.rbind_result_dfs(dfs)
        out_df <- out_df[!duplicated(rownames(out_df)), , drop = FALSE]
        out <- out_df$Estimate
        names(out) <- rownames(out_df)
        out
      }

      update_df_with_boot <- function(df, boot_se, boot_lower, boot_upper) {
        if (!is.data.frame(df)) return(df)

        common <- intersect(rownames(df), names(boot_se))
        if (length(common) == 0L) return(df)

        df[common, "Std. Error"] <- boot_se[common]
        df[common, "Lower 95%"] <- boot_lower[common]
        df[common, "Upper 95%"] <- boot_upper[common]

        stat_col <- if ("z value" %in% names(df)) {
          "z value"
        } else if ("t value" %in% names(df)) {
          "t value"
        } else {
          NULL
        }

        if (!is.null(stat_col)) {
          df[common, stat_col] <- df[common, "Estimate"] /
            pmax(df[common, "Std. Error"], 1e-12)
        }

        if ("Pr" %in% names(df) && !is.null(stat_col)) {
          if ("df" %in% names(df)) {
            df[common, "Pr"] <- 2 * stats::pt(
              -abs(df[common, stat_col]),
              df = df[common, "df"]
            )
          } else {
            df[common, "Pr"] <- 2 * stats::pnorm(
              -abs(df[common, stat_col])
            )
          }
        }

        df
      }

      compute_mediation_boot_estimates <- function(dat_b, target_names, formula, family) {
        n_eq <- length(formula)
        if (!is.list(family)) family <- rep(list(family), n_eq)
        if (length(family) != n_eq) stop("Length of family list must match length of formula list.")

        coefs <- vector("list", n_eq)
        x_names <- vector("list", n_eq)
        resp_names <- character(n_eq)
        out <- stats::setNames(rep(NA_real_, length(target_names)), target_names)

        for (j in seq_len(n_eq)) {
          f <- formula[[j]]
          fam <- family[[j]]
          fam_name <- if (inherits(fam, "family")) fam$family else as.character(fam)[1]
          mf <- stats::model.frame(f, data = dat_b)
          y <- stats::model.response(mf)
          X <- stats::model.matrix(f, data = mf)
          cn <- colnames(X)
          cn[cn == "(Intercept)"] <- "Intercept"
          x_names[[j]] <- cn
          resp_names[j] <- as.character(f[[2]])

          if (fam_name == "gaussian") {
            fit_j <- stats::lm.fit(x = X, y = as.numeric(y))
            beta <- as.numeric(fit_j$coefficients)
            names(beta) <- cn
            coefs[[j]] <- beta
            b_names <- paste0("b", j, "[", cn, "]")
            out[intersect(b_names, target_names)] <- beta[match(intersect(b_names, target_names), b_names)]
            sigma_name <- paste0("sigma", j)
            if (sigma_name %in% target_names) {
              out[[sigma_name]] <- sqrt(mean(fit_j$residuals^2, na.rm = TRUE))
            }
          } else if (fam_name %in% c("binomial", "bernoulli")) {
            fit_j <- stats::glm.fit(x = X, y = y, family = stats::binomial())
            beta <- as.numeric(fit_j$coefficients)
            names(beta) <- cn
            coefs[[j]] <- beta
            b_names <- paste0("b", j, "[", cn, "]")
            out[intersect(b_names, target_names)] <- beta[match(intersect(b_names, target_names), b_names)]
          } else if (fam_name == "poisson") {
            fit_j <- stats::glm.fit(x = X, y = y, family = stats::poisson())
            beta <- as.numeric(fit_j$coefficients)
            names(beta) <- cn
            coefs[[j]] <- beta
            b_names <- paste0("b", j, "[", cn, "]")
            out[intersect(b_names, target_names)] <- beta[match(intersect(b_names, target_names), b_names)]
          } else {
            stop(sprintf("Unsupported mediation bootstrap family: %s", fam_name), call. = FALSE)
          }
        }

        mediators <- intersect(resp_names, unique(unlist(x_names)))
        indeps <- setdiff(unique(unlist(x_names)), c(resp_names, "Intercept"))
        if (length(indeps) == 0L && length(mediators) > 0L) {
          indeps <- setdiff(x_names[[1]], "Intercept")
        }

        for (iv in indeps) {
          for (m in mediators) {
            idx_m_resp <- which(resp_names == m)
            idx_m_pred <- which(vapply(x_names, function(x) m %in% x, logical(1)))
            if (length(idx_m_resp) == 0L || length(idx_m_pred) == 0L) next

            for (m_resp in idx_m_resp) {
              if (!(iv %in% names(coefs[[m_resp]]))) next
              a_val <- coefs[[m_resp]][[iv]]

              for (dv_idx in idx_m_pred) {
                if (dv_idx == m_resp) next
                dv_name <- resp_names[[dv_idx]]
                if (!(m %in% names(coefs[[dv_idx]]))) next
                b_val <- coefs[[dv_idx]][[m]]

                ie_name <- paste0("IE_", iv, "_", m, "_", dv_name)
                if (ie_name %in% target_names) out[[ie_name]] <- a_val * b_val

                if (iv %in% names(coefs[[dv_idx]])) {
                  de_val <- coefs[[dv_idx]][[iv]]
                  de_name <- paste0("DE_", iv, "_", dv_name)
                  te_name <- paste0("TE_", iv, "_", m, "_", dv_name)
                  if (de_name %in% target_names) out[[de_name]] <- de_val
                  if (te_name %in% target_names) out[[te_name]] <- a_val * b_val + de_val
                }
              }
            }
          }
        }

        out
      }

      use_fast_mediation_bootstrap <- function(family) {
        fam_list <- if (is.list(family)) family else list(family)
        fam_names <- vapply(fam_list, function(fam) {
          if (inherits(fam, "family")) fam$family else as.character(fam)[1]
        }, character(1))
        all(fam_names %in% c("gaussian", "binomial", "bernoulli", "poisson"))
      }

      base_est <- extract_estimates(self)
      target_names <- names(base_est)
      if (length(target_names) == 0L) {
        stop("No estimates are available for bootstrap.", call. = FALSE)
      }

      med_info <- self$model$extra$mediation %||% list()
      formula <- med_info$formula %||% self$model$formula
      family <- med_info$family %||% self$model$family %||% "gaussian"
      view <- med_info$view %||% self$model$view
      dat <- self$model$raw_data

      if (is.null(formula) || is.null(dat)) {
        stop(
          "Mediation bootstrap requires the original formula and raw data stored in the model.",
          call. = FALSE
        )
      }
      dat <- as.data.frame(dat)
      n <- nrow(dat)
      if (is.na(n) || n <= 0L) stop("Original data are empty.", call. = FALSE)

      if (!is.null(seed)) set.seed(seed)

      boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(target_names))
      colnames(boot_mat) <- target_names
      ok <- rep(FALSE, n_boot)
      fast_boot <- use_fast_mediation_bootstrap(family)

      progress_mode <- .rtmb_resolve_progress("auto")
      .rtmb_progress_line(
        sprintf("Performing mediation bootstrap estimation (n = %d)...", n_boot),
        progress_mode
      )
      meter <- .rtmb_progress_meter(n_boot, progress_mode, label = "bootstrap")

      old_silent <- options(BayesRTMB.silent = TRUE)
      on.exit(options(old_silent), add = TRUE)

      for (i in seq_len(n_boot)) {
        idx <- sample.int(n, size = n, replace = TRUE)
        dat_b <- dat[idx, , drop = FALSE]

        boot_est <- tryCatch({
          if (fast_boot) {
            compute_mediation_boot_estimates(
              dat_b = dat_b,
              target_names = target_names,
              formula = formula,
              family = family
            )
          } else {
            boot_mdl <- rtmb_mediation(
              formula = formula,
              data = dat_b,
              family = family,
              prior = prior_flat(),
              view = view
            )
            boot_fit <- boot_mdl$classic()
            extract_estimates(boot_fit)
          }
        }, error = function(e) NULL)

        if (!is.null(boot_est) && length(boot_est) > 0L) {
          common <- intersect(target_names, names(boot_est))
          if (length(common) > 0L) {
            boot_mat[i, common] <- boot_est[common]
            ok[i] <- TRUE
          }
        }

        meter$advance(1L)
      }
      meter$finish()

      if (!any(ok)) {
        stop("All bootstrap fits failed.", call. = FALSE)
      }
      if (sum(ok) < n_boot) {
        warning(
          sprintf("Only %d of %d bootstrap fits succeeded.", sum(ok), n_boot),
          call. = FALSE
        )
      }

      boot_mat_ok <- boot_mat[ok, , drop = FALSE]
      boot_se <- apply(boot_mat_ok, 2, stats::sd, na.rm = TRUE)
      boot_lower <- apply(boot_mat_ok, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
      boot_upper <- apply(boot_mat_ok, 2, stats::quantile, probs = 0.975, na.rm = TRUE)

      boot_lower <- as.numeric(boot_lower)
      boot_upper <- as.numeric(boot_upper)
      names(boot_lower) <- names(boot_upper) <- colnames(boot_mat_ok)

      self$bootstrap_results <- boot_mat_ok
      self$se_method <- "bootstrap"
      self$ci_method <- "bootstrap"
      self$cluster <- NULL
      self$robust_type <- NULL
      self$vcov <- stats::cov(boot_mat_ok, use = "pairwise.complete.obs")

      self$df_fixed <- update_df_with_boot(self$df_fixed, boot_se, boot_lower, boot_upper)
      self$random_effects <- update_df_with_boot(self$random_effects, boot_se, boot_lower, boot_upper)
      self$df_transform <- update_df_with_boot(self$df_transform, boot_se, boot_lower, boot_upper)
      self$df_generate <- update_df_with_boot(self$df_generate, boot_se, boot_lower, boot_upper)
      self$fit <- self$df_fixed

      return(self)
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

    #' @description WAIC is not defined for classical fits.
    WAIC = function() {
      stop(
        "WAIC requires posterior draws of pointwise log_lik. ",
        "Use sample() or variational() on a model created with WAIC = TRUE.",
        call. = FALSE
      )
    },

    #' @description Run basic diagnostics for the classical fit.
    #' @param ... Additional arguments passed to `diagnose_classic_fit()`.
    #' @return A `diagnose_BayesRTMB` object.
    diagnose = function(...) {
      diagnose_classic_fit(self, ...)
    },

    #' @description Print the fit results.
    #' @param ... Additional arguments passed to `summary()`.
    print = function(...) {
      print(self$summary(...))
      invisible(self)
    },

    #' @description Get the Log-Likelihood of the fitted model.
    logLik = function() {
      if (inherits(self$fit, c("lm", "glm"))) {
        return(stats::logLik(self$fit))
      }

      val <- if (!is.null(self$log_lik)) {
        self$log_lik
      } else if (!is.null(self$model$extra[["loglik"]])) {
        self$model$extra[["loglik"]]
      } else {
        NA_real_
      }

      df_val <- if (!is.null(self$model$extra[["df"]])) {
        self$model$extra[["df"]]
      } else if (is.data.frame(self$fit)) {
        nrow(self$fit)
      } else {
        length(self$par_vec)
      }

      nobs_val <- if (!is.null(self$model$extra[["nobs"]])) {
        self$model$extra[["nobs"]]
      } else if (!is.null(self$model$raw_data)) {
        nrow(as.data.frame(self$model$raw_data))
      } else {
        NA_real_
      }

      res <- val
      attr(res, "df") <- df_val
      attr(res, "nobs") <- nobs_val
      class(res) <- "logLik"
      res
    },

    #' @description EAP estimates are not available for Classic_Fit.
    #' @param ... Ignored.
    EAP = function(...) {
      stop(
        "EAP() is not available for Classic_Fit. Use estimate() instead.",
        call. = FALSE
      )
    },

    #' @description MAP estimates are not available for Classic_Fit.
    #' @param ... Ignored.
    MAP = function(...) {
      stop(
        "MAP() is not available for Classic_Fit. Use estimate() instead.",
        call. = FALSE
      )
    },

    #' @description Display a summary of the estimation results.
    #' @param view Character vector of parameter names to prioritize or filter by.
    #' @param digits Number of digits to print for estimates.
    #' @param max_rows Maximum number of rows to display in the coefficient table.
    #' @param ranef Logical; whether to show random effects in the summary.
    summary = function(view = NULL, digits = 5, max_rows = 10, ranef = FALSE) {
      if (!is.numeric(digits)) digits <- 5
      res <- list(
        type = self$model$type,
        family = self$model$family,
        se_method = self$se_method,
        cluster = self$cluster,
        robust_type = self$robust_type,
        bootstrap = if (!is.null(self$bootstrap_results)) nrow(self$bootstrap_results) else NULL,
        logLik = self$logLik(),
        AIC = self$AIC(),
        BIC = self$BIC(),
        extra = self$model$extra,
        test_results = self$test_results,
        digits = digits
      )

      if (is.data.frame(self$fit) || inherits(self$fit, "lm")) {
        # Combine all result dataframes (fixed, random, transform, generate)
        all_dfs <- list()
        if (is.data.frame(self$fit)) all_dfs$fixed <- self$fit
        
        if (!is.null(self$random_effects) && is.data.frame(self$random_effects)) {
           all_dfs$random <- self$random_effects
        }
        if (!is.null(self$df_transform) && is.data.frame(self$df_transform)) {
           all_dfs$transform <- self$df_transform
        }
        if (!is.null(self$df_generate) && is.data.frame(self$df_generate)) {
           all_dfs$generate <- self$df_generate
        }

        if (length(all_dfs) > 0) {
           df_all <- private$.rbind_result_dfs(all_dfs)
           # Remove duplicates (preferring earlier ones)
           df_all <- df_all[!duplicated(rownames(df_all)), , drop = FALSE]
           # Track which rows are random effects
           is_re <- if (!is.null(self$random_effects)) rownames(df_all) %in% rownames(self$random_effects) else rep(FALSE, nrow(df_all))
        } else if (inherits(self$fit, "lm")) {
           # Handle lm case specifically
           s_lm <- summary(self$fit)
           df_all <- as.data.frame(s_lm$coefficients)
        } else {
           df_all <- self$fit
        }

        # --- Correlation Model Filtering ---
        if (!is.null(self$model$type) && self$model$type == "corr") {
           keep <- rep(TRUE, nrow(df_all))
           seen_pairs <- list()
           for (i in 1:nrow(df_all)) {
             nm <- rownames(df_all)[i]
             est <- as.numeric(df_all$Estimate[i])
             se <- as.numeric(df_all$`Std. Error`[i])
             if (!is.na(est) && abs(est - 1) < 1e-7 && (is.na(se) || se < 1e-10)) {
               keep[i] <- FALSE
               next
             }
             if (grepl("\\[.*,.*\\]", nm)) {
               parts <- gsub(".*\\[(.*)\\]", "\\1", nm); bits <- strsplit(parts, ",")[[1]]
               if (length(bits) == 2) {
                 pair_id <- paste(sort(trimws(bits)), collapse = "_")
                 prefix <- gsub("\\[.*", "", nm); pair_key <- paste0(prefix, "_", pair_id)
                 if (pair_key %in% names(seen_pairs)) { keep[i] <- FALSE } else {
                   seen_pairs[[pair_key]] <- TRUE
                   if (trimws(bits[1]) == trimws(bits[2])) keep[i] <- FALSE
                 }
               }
             }
           }
           df_all <- df_all[keep, , drop = FALSE]
           if (exists("is_re")) is_re <- is_re[keep]
        }

        df_print <- df_all
        
        # --- Handle lm specifically if needed ---
        if (inherits(self$fit, "lm")) {
           s_lm <- summary(self$fit)
           is_asymptotic <- inherits(self$model, c("rtmb_loglinear", "rtmb_table"))
           df_print <- as.data.frame(s_lm$coefficients)
           if (!is.null(self$vcov)) {
             new_se <- sqrt(diag(self$vcov))
             n_match <- min(nrow(df_print), length(new_se))
             df_print[1:n_match, 2] <- new_se[1:n_match]
             df_print[1:n_match, 3] <- df_print$Estimate[1:n_match] / pmax(df_print[1:n_match, 2], 1e-12)
           }
           if (is_asymptotic) {
             colnames(df_print) <- c("Estimate", "Std. Error", "z value", "Pr")
             df_print$df <- Inf
           } else {
             colnames(df_print) <- c("Estimate", "Std. Error", "t value", "Pr")
             df_print$df <- self$fit$df.residual
           }
           rownames(df_print)[rownames(df_print) == "(Intercept)"] <- "Intercept"
        }

        # --- View Prioritization Logic ---
        # Distinction: if 'view' is explicitly passed, we apply STRICT filtering.
        view_explicit <- !is.null(view)
        view_to_use <- if (view_explicit) view else if (!is.null(self$view)) self$view else self$model$view
        
        # If it's a correlation model, default is c("corr", "mean", "sd")
        if (is.null(view_to_use) && !is.null(self$model$type) && self$model$type == "corr") {
           view_to_use <- c("corr", "mean", "sd")
        }
        
        priority_idx <- integer(0)
        cf_idx <- integer(0)
        
        if (!is.null(view_to_use)) {
          var_names <- rownames(df_print)
          base_names <- gsub("\\[.*\\]$", "", var_names)
          
          # Pass 1: Exact matches and base names for ALL keywords
          for (v in view_to_use) {
            m1 <- which(var_names == v | base_names == v)
            priority_idx <- c(priority_idx, setdiff(m1, priority_idx))
          }
          
          # Pass 2: Standard prefix matches (e.g. corr[...]) for ALL keywords
          for (v in view_to_use) {
            # Specifically exclude CF_ prefix here to prevent it from jumping ahead
            m2 <- which(grepl(paste0("^", v), var_names) & !grepl("^CF_", var_names))
            priority_idx <- c(priority_idx, setdiff(m2, priority_idx))
          }
          
          # Pass 3: Internal CF_ matches (e.g. CF_corr[...]) for ALL keywords - ALWAYS LAST
          for (v in view_to_use) {
            m3 <- which(grepl(paste0("^CF_", v), var_names))
            cf_idx <- c(cf_idx, setdiff(m3, cf_idx))
          }
          
          priority_idx <- unique(priority_idx)
          cf_idx <- setdiff(unique(cf_idx), priority_idx)
          
          # Reorder based on hierarchy: Standard matches > CF matches > Others
          final_order <- c(priority_idx, cf_idx)
          
          if (view_explicit) {
             # STRICT FILTERING: Only keep matched rows
             other_idx <- integer(0)
          } else {
             # PRIORITIZATION ONLY: Keep everything else at the bottom
             other_idx <- setdiff(seq_len(nrow(df_print)), final_order)
          }
          final_order <- c(final_order, other_idx)
          
          # Update is_re and df_print
          if (exists("is_re")) is_re <- is_re[final_order]
          matched_indices <- seq_along(priority_idx)
          if (length(cf_idx) > 0) matched_indices <- c(matched_indices, (length(priority_idx)+1):(length(priority_idx)+length(cf_idx)))
          
          df_print <- df_print[final_order, , drop = FALSE]
        }
        
        # --- Random Effect Filtering ---
        if (!isTRUE(ranef) && exists("is_re")) {
           # Keep if NOT a random effect OR it was explicitly matched by view
           keep_mask <- !is_re
           if (exists("matched_indices") && length(matched_indices) > 0) {
              keep_mask[matched_indices] <- TRUE
           }
           df_print <- df_print[keep_mask, , drop = FALSE]
        }
        
        # --- Truncation Logic ---
        if (!is.null(max_rows) && nrow(df_print) > max_rows) {
           res$truncated <- TRUE
           res$total_rows <- nrow(df_print)
           res$max_rows <- max_rows
        } else {
           res$truncated <- FALSE
        }

        is_asymptotic <- (!is.null(self$model$type) &&
                            self$model$type %in% c("glm", "glmer", "table", "loglinear", "fa", "irt", "mixture")) ||
          inherits(self$model, "rtmb_loglinear")
        hide_tests <- !is.null(self$model$type) && self$model$type %in% c("mdu", "fa")

        # --- Add t/z-test results for classic mode ---
        if (!is.null(self$sd_rep)) {
           u_est <- if (!is.null(self$sd_rep$par.fixed)) self$sd_rep$par.fixed else numeric(0)
           u_se <- if (!is.null(self$sd_rep$cov.fixed)) sqrt(pmax(diag(self$sd_rep$cov.fixed), 0)) else rep(NA, length(u_est))
           u_t <- ifelse(is.na(u_se), NA_real_, u_est / pmax(u_se, 1e-12))
           
           # Map to df_print rows
           u_idx_map <- list()
           if (length(u_est) > 0) {
             for (nm in unique(names(u_est))) {
                u_idx_map[[nm]] <- which(names(u_est) == nm)
             }
           }
           
           u_counter <- list()
           row_t <- rep(NA, nrow(df_print))
           
           for (i in 1:nrow(df_print)) {
              row_nm <- rownames(df_print)[i]
              base_nm <- gsub("\\[.*\\]$", "", row_nm)
              
              if (base_nm %in% names(u_idx_map)) {
                 if (is.null(u_counter[[base_nm]])) u_counter[[base_nm]] <- 1
                 curr_idx <- u_counter[[base_nm]]
                 if (curr_idx <= length(u_idx_map[[base_nm]])) {
                    real_idx <- u_idx_map[[base_nm]][curr_idx]
                    row_t[i] <- u_t[real_idx]
                    u_counter[[base_nm]] <- curr_idx + 1
                 }
              }
           }
           
           # Fallback to constrained ratio for parameters not in par.fixed (e.g. transformed)
           missing_t <- is.na(row_t)
           if (any(missing_t)) {
              row_t[missing_t] <- ifelse(is.na(df_print$`Std. Error`[missing_t]), NA_real_,
                df_print$Estimate[missing_t] / pmax(df_print$`Std. Error`[missing_t], 1e-12))
           }

           if (identical(self$se_method, "bootstrap")) {
              row_t <- ifelse(is.na(df_print$`Std. Error`), NA_real_,
                df_print$Estimate / pmax(df_print$`Std. Error`, 1e-12))
           } else if (self$se_method %in% c("robust", "cluster-robust") && !is.null(self$vcov)) {
              robust_rows <- intersect(rownames(df_print), rownames(self$vcov))
              if (length(robust_rows) > 0L) {
                robust_idx <- match(robust_rows, rownames(df_print))
                row_t[robust_idx] <- ifelse(is.na(df_print$`Std. Error`[robust_idx]), NA_real_,
                  df_print$Estimate[robust_idx] / pmax(df_print$`Std. Error`[robust_idx], 1e-12))
              }
           } else if (!identical(self$se_method, "wald")) {
              row_t <- ifelse(is.na(df_print$`Std. Error`), NA_real_,
                df_print$Estimate / pmax(df_print$`Std. Error`, 1e-12))
           }
           
           if (is_asymptotic) {
                df_print$`z value` <- row_t
                if ("df" %in% names(df_print)) df_print$df <- NULL # Hide DF for asymptotic models
             } else {
                df_print$`t value` <- row_t
             }
             
             if (!is.null(df_print$df)) {
                df_print$Pr <- 2 * pt(-abs(row_t), df = df_print$df)
             } else {
                df_print$Pr <- 2 * pnorm(-abs(row_t))
             }

             if (!is.null(self$model$type) && self$model$type == "ttest") {
                test_rows <- rownames(df_print) == "diff"
                if ("t value" %in% names(df_print)) df_print$`t value`[!test_rows] <- NA_real_
                if ("z value" %in% names(df_print)) df_print$`z value`[!test_rows] <- NA_real_
                df_print$Pr[!test_rows] <- NA_real_
             }

             if (!is.null(self$model$type) && self$model$type %in% c("lm", "glm", "lmer", "glmer", "mediation")) {
                row_base <- gsub("\\[.*\\]$", "", rownames(df_print))
                scale_rows <- grepl("^sigma", row_base) | grepl("^sd", row_base)
                if (any(scale_rows)) {
                   if ("t value" %in% names(df_print)) df_print$`t value`[scale_rows] <- NA_real_
                   if ("z value" %in% names(df_print)) df_print$`z value`[scale_rows] <- NA_real_
                   df_print$Pr[scale_rows] <- NA_real_
                }
             }

             if (hide_tests) {
                df_print$`z value` <- NULL
                df_print$`t value` <- NULL
                df_print$Pr <- NULL
             }
        }

        # --- Hide DF column if requested ---
        if (!isTRUE(self$show_df) && "df" %in% names(df_print)) {
           df_print$df <- NULL
        }

        # 1. Round main numeric columns
        cols_to_round <- setdiff(names(df_print), c("Pr", "df"))
        for (col in cols_to_round) {
          if (is.numeric(df_print[[col]])) {
            df_print[[col]] <- round(df_print[[col]], digits)
          }
        }

        # 2. Round df to 1 decimal place
        if ("df" %in% names(df_print) && is.numeric(df_print$df)) {
          df_print$df <- round(df_print$df, 1)
        }

        # 3. Format Pr
        if (!is.null(df_print$Pr)) {
          # Standard significance symbols
          sig <- symnum(as.numeric(df_print$Pr), corr = FALSE, na = FALSE,
                        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                        symbols = c("***", "**", "*", ".", " "))

          # Format Pr values (e.g., .0000, .0123)
          formatted_pr <- sapply(df_print$Pr, function(p) {
            if (is.na(p)) return("")
            p_val <- as.numeric(p)
            if (p_val < 0.0001) return("<.0001")
             fmt <- paste0("%.", digits, "f")
             s <- sprintf(fmt, p_val)
             sub("^0", "", s) # Remove leading zero
          })

          df_print$Pr <- formatted_pr
          df_print$sig <- as.character(sig)
        }

        # 4. Reorder columns for consistency
        if (is_asymptotic) {
          desired_order <- c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "df", "z value", "Pr")
        } else {
          desired_order <- c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "df", "t value", "Pr")
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
      # --- 1. Restricted Model Types ---
      restricted_types <- c("ttest", "corr", "fa", "irt")
      if (!is.null(self$model$type) && self$model$type %in% restricted_types) {
        stop(sprintf("anova() is not supported for '%s' models.", self$model$type), call. = FALSE)
      }

      if (!is.null(self$model$type) && self$model$type == "table") {
        obs_tab <- self$model$extra$tab
        N_tot <- sum(obs_tab)
        R_dim <- nrow(obs_tab)
        C_dim <- ncol(obs_tab)

        # Combine all result dataframes to find 'mu'
        all_dfs <- list()
        if (!is.null(self$df_fixed)) all_dfs$fixed <- self$df_fixed
        if (!is.null(self$df_transform)) all_dfs$transform <- self$df_transform
        if (!is.null(self$df_generate)) all_dfs$generate <- self$df_generate
        
        df_combined <- if (length(all_dfs) > 0) {
          tmp <- private$.rbind_result_dfs(all_dfs)
          tmp[!duplicated(rownames(tmp)), , drop = FALSE]
        } else data.frame()
        
        df_mu <- df_combined[grepl("^mu(\\[|$)", rownames(df_combined)), , drop = FALSE]

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

           chisq_res <- data.frame(
             `X-squared` = round(chisq_stat, 5), 
             df = round(df_val, 1), 
             `p-value` = round(p_val, 5), 
             check.names = FALSE
           )
           rownames(chisq_res) <- "Pearson's Chi-squared (from est. prob)"
        } else {
          chisq_res <- NULL
        }

        fisher_res <- if (!is.null(self$test_results$fisher)) self$test_results$fisher else self$model$extra$fisher
        res <- list(chisq = chisq_res, fisher = fisher_res)
        class(res) <- "anova_rtmb_table"
        return(res)
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
        # Expand patterns for identifying fixed effects (adding beta, mu, delta, etc.)
        fix_pats <- c("Intercept", "Intercept_c", "b", "mean", "prob", "beta", "mu", "delta", "diff", "mean_diff")
        fe_regex <- paste0("^(", paste(fix_pats, collapse="|"), ")($|\\[)")
        fe_idx <- which(grepl(fe_regex, names(beta_full)))
      }

      if (length(fe_idx) == 0) stop("Could not identify fixed effects for ANOVA.")

      beta <- beta_full[fe_idx]
      fe_names_actual <- names(beta)
      V <- V_full[fe_names_actual, fe_names_actual, drop = FALSE]

      full_assign <- c(0, assign_idx)

      ct_setting <- if (!is.null(self$model$contrasts)) self$model$contrasts else "sum"
      if (ct_setting != "sum") {
        formula <- nobars(self$model$formula)
        predictor_vars <- all.vars(delete.response(terms(formula)))
        relevant_data <- self$model$raw_data[, predictor_vars, drop = FALSE]
        if (is.null(relevant_data)) relevant_data <- self$model$data[, predictor_vars, drop = FALSE]
        levs <- lapply(relevant_data, function(x) {
          if (is.factor(x) || is.character(x)) {
            levels(as.factor(x))
          } else {
            mean(x, na.rm = TRUE)
          }
        })
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

        # --- Special Case: Calculate F-value by squaring t-value for 1-DF terms (continuous or 2-level factors) ---
        # Attempt direct extraction from original beta_full and V_full to avoid contrast transformation issues
        is_single <- length(idx) == 1
        W <- NULL
        if (is_single) {
          term_name_in_beta <- names(beta)[idx]
          # Look for the corresponding coefficient in the original estimation results (beta_full)
          # Try both b[name] format and the name itself
          orig_match <- which(names(beta_full) == term_name_in_beta | names(beta_full) == paste0("b[", term_name_in_beta, "]"))
          if (length(orig_match) == 1) {
            b_val <- beta_full[orig_match]
            se_val <- sqrt(V_full[orig_match, orig_match])
            W <- (b_val / se_val)^2
          }
        }

        # Standard Wald test if not 1-DF or if the above failed
        if (is.null(W)) {
          L <- matrix(0, nrow = length(idx), ncol = length(beta))
          for (j in seq_along(idx)) L[j, idx[j]] <- 1

          LVL <- L %*% V %*% t(L)
          inv_LVL <- try(solve(LVL), silent = TRUE)
          if (inherits(inv_LVL, "try-error")) inv_LVL <- MASS::ginv(LVL)

          W <- as.numeric(t(L %*% beta) %*% inv_LVL %*% (L %*% beta))
        }

        df1 <- length(idx)
        f_val <- W / df1

        fit_df_col <- if ("df" %in% names(self$fit)) self$fit$df else if ("DF" %in% names(self$fit)) self$fit$DF else NULL
        if (is.data.frame(self$fit) && !is.null(fit_df_col)) {
          # Align with beta_full indices
          term_name_in_beta <- names(beta)[idx]
          orig_match_indices <- sapply(term_name_in_beta, function(nm) {
            m <- which(names(beta_full) == nm | names(beta_full) == paste0("b[", nm, "]"))
            if (length(m) > 0) m[1] else NA
          })
          df2 <- if (any(!is.na(orig_match_indices))) min(fit_df_col[orig_match_indices], na.rm = TRUE) else Inf
        } else if (inherits(self$fit, "lm")) {
          df2 <- self$fit$df.residual
        } else {
          df2 <- Inf
        }

        # --- Decision: Use F-test or Chi-squared test format ---
        # Unify to F-test format unless it is an asymptotic model (Poisson, Binomial, Table, etc.)
        is_asymp <- (!is.null(self$model$family) && self$model$family %in% c("poisson", "binomial", "bernoulli", "multinomial", "ordered")) || 
                    (!is.null(self$model$type) && self$model$type %in% c("loglinear", "table"))
        
        if (is_asymp) {
          # Asymptotic model: Chi-squared test format
          p_val <- stats::pchisq(W, df1, lower.tail = FALSE)
          res_list[[t_name]] <- data.frame(`df` = df1, `Chisq` = W, `Pr(>Chisq)` = p_val, row.names = t_name, check.names = FALSE)
        } else {
          # Regression models: Unify to F-test format (including cases where denominator DF is Inf)
          p_val <- if (is.infinite(df2) || df2 <= 0) stats::pchisq(W, df1, lower.tail = FALSE) 
                   else stats::pf(f_val, df1, df2, lower.tail = FALSE)
          res_list[[t_name]] <- data.frame(`num_df` = df1, `den_df` = df2, `F value` = f_val, `Pr(>F)` = p_val, row.names = t_name, check.names = FALSE)
        }
      }

      res_df <- do.call(rbind, res_list)

      class(res_df) <- c("anova_rtmb", "anova", "data.frame")
      heading <- if ("Chisq" %in% names(res_df)) "ANOVA Table (Wald Chisq tests)" else "ANOVA Table (Wald F-tests)"
      attr(res_df, "heading") <- heading

      # --- Added: R-squared and overall model test (for rtmb_lm) ---
      is_lm <- !is.null(self$model$type) && self$model$type == "lm"
      if (is_lm && !is.null(self$model$data$Y)) {
        Y <- self$model$data$Y
        TSS <- sum((Y - mean(Y))^2)
        sigma_est <- beta_full["sigma"]
        if (!is.na(sigma_est) && TSS > 0) {
          # Obtain residual degrees of freedom
          df_resid <- if (is.null(res_df$den_df)) (length(Y) - sum(res_df$num_df) - 1) else min(res_df$den_df, na.rm = TRUE)
          if (is.infinite(df_resid)) df_resid <- length(Y) - sum(res_df$num_df) - 1
          
          RSS <- (sigma_est^2) * df_resid
          r2 <- 1 - (RSS / TSS)
          adj_r2 <- 1 - (1 - r2) * (length(Y) - 1) / pmax(df_resid, 1)

          # Overall model Wald test (all predictors vs. Intercept)
          all_pred_idx <- which(assign_idx > 0)
          if (length(all_pred_idx) > 0) {
            L_all <- matrix(0, nrow = length(all_pred_idx), ncol = length(beta))
            for (j in seq_along(all_pred_idx)) L_all[j, all_pred_idx[j]] <- 1
            LVL_all <- L_all %*% V %*% t(L_all)
            W_all <- as.numeric(t(L_all %*% beta) %*% MASS::ginv(LVL_all) %*% (L_all %*% beta))
            f_all <- W_all / length(all_pred_idx)
            p_all <- stats::pf(f_all, length(all_pred_idx), df_resid, lower.tail = FALSE)
            
            attr(res_df, "lm_info") <- list(
              r_squared = r2,
              adj_r_squared = adj_r2,
              fstatistic = c(value = f_all, numdf = length(all_pred_idx), dendf = df_resid),
              p_value = p_all,
              sigma = sigma_est,
              df_resid = df_resid
            )
          }
        }
      }
      return(res_df)
    },

    #' @description Calculate Least Squares Means (Marginal Means) and contrasts.
    #' @param specs Character vector of factors to calculate means for.
    #' @param pairwise Logical; whether to perform pairwise comparisons.
    #' @param simple Character vector of factors to hold constant for simple main effects.
    #' @param adjust Character; p-value adjustment method (e.g., "bonferroni", "holm", "none").
    #' @param protect Logical; whether to use hierarchical (protected) testing.
    #' @return A data frame containing the marginal means or contrasts.
    lsmeans = function(specs = NULL, pairwise = FALSE, simple = NULL, adjust = "holm", protect = FALSE) {
      # --- 1. Restricted Model Types ---
      restricted_types <- c("ttest", "corr", "fa", "irt", "table")
      if (!is.null(self$model$type) && self$model$type %in% restricted_types) {
        stop(sprintf("lsmeans() is not supported for '%s' models.", self$model$type), call. = FALSE)
      }

      # --- 2. Check for missing specs ---
      if (is.null(specs)) {
        formula <- self$model$formula
        data <- if (!is.null(self$model$raw_data)) as.data.frame(self$model$raw_data) else as.data.frame(self$model$data)
        vars <- all.vars(nobars(formula))[-1]
        is_cat <- sapply(data[vars], function(x) is.factor(x) || is.character(x))
        cat_vars <- vars[is_cat]
        
        if (length(cat_vars) == 0) {
          stop("lsmeans() is only available for models with categorical factors.", call. = FALSE)
        } else {
          stop(sprintf("Please specify 'specs' (e.g., fit$lsmeans('%s')).\nAvailable factors: %s", 
                       cat_vars[1], paste(cat_vars, collapse = ", ")), call. = FALSE)
        }
      }

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
      # Reorder variables so that specs vary fastest and simple varies slowest among the targeted specs
      other_vars <- setdiff(cat_vars, full_specs)
      ordered_grid_vars <- c(other_vars, full_specs)
      grid_list <- lapply(data[ordered_grid_vars], function(x) levels(as.factor(x)))
      ref_grid <- expand.grid(grid_list)

      # 3. Handle continuous covariates by setting them to their mean
      cont_vars <- vars[!is_cat]
      for (v in cont_vars) {
        ref_grid[[v]] <- mean(data[[v]], na.rm = TRUE)
      }

      # 4. Generate model matrix for the reference grid
      mf_orig <- model.frame(nobars(formula), data)
      orig_terms <- delete.response(terms(mf_orig))

      ct_setting <- if (!is.null(self$model$contrasts)) self$model$contrasts else "sum"
      ct_list <- list()
      for (v in cat_vars) {
        ct_list[[v]] <- if (ct_setting == "treatment") "contr.treatment" else "contr.sum"
      }

      X_grid_raw <- model.matrix(orig_terms, ref_grid, contrasts.arg = ct_list)

      # Match columns with the original model
      orig_cols <- self$model$extra$X_colnames
      if (!is.null(orig_cols)) {
        grid_cols <- colnames(X_grid_raw)
        # Normalize BOTH to 'Intercept' for matching
        grid_cols[grid_cols == "(Intercept)"] <- "Intercept"
        orig_cols_norm <- orig_cols
        orig_cols_norm[orig_cols_norm == "(Intercept)"] <- "Intercept"
        
        col_idx <- match(orig_cols_norm, grid_cols)
        if (any(is.na(col_idx))) {
          missing_cols <- orig_cols[is.na(col_idx)]
          stop(sprintf("Could not match lsmeans grid columns.\nMissing from grid: %s\nAvailable in grid: %s", 
                       paste(missing_cols, collapse = ", "), 
                       paste(grid_cols, collapse = ", ")), call. = FALSE)
        }
        X_grid <- X_grid_raw[, col_idx, drop = FALSE]
      } else {
        X_grid <- X_grid_raw
      }

      beta_vals <- if (is.data.frame(self$fit)) self$fit$Estimate else stats::coef(self$fit)
      fe_idx <- which(grepl("^(Intercept|Intercept_c|b\\[)", rownames(self$fit)))
      beta_match <- beta_vals[fe_idx]
      fe_names_match <- rownames(self$fit)[fe_idx]
      V_match <- self$vcov[fe_names_match, fe_names_match, drop = FALSE]

      # Group by specs + simple with descriptive labels
      label_grid <- ref_grid[full_specs]
      for (v in full_specs) {
        label_grid[[v]] <- paste0(v, "=", label_grid[[v]])
      }
      groups_full <- interaction(label_grid, drop = TRUE, sep = ":")
      levels(groups_full) <- unique(as.character(groups_full))
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
        # --- Pairwise Comparisons & Simple Main Effects ---
        results_by_group <- list()
        
        # Moderator levels
        if (!is.null(simple)) {
          mod_grid <- unique(ref_grid[simple])
          focal_grid <- unique(ref_grid[specs])

          for (row_m in 1:nrow(mod_grid)) {
            m_vals <- mod_grid[row_m, , drop = FALSE]
            m_label <- paste(sapply(1:ncol(m_vals), function(k) paste0(names(m_vals)[k], "=", m_vals[1, k])), collapse = ":")
            
            # --- Simple Main Effect (SME) ---
            L_group <- list()
            group_names <- c()
            for (i in 1:nrow(focal_grid)) {
                f_label <- paste(sapply(1:ncol(focal_grid), function(k) paste0(names(focal_grid)[k], "=", focal_grid[i, k])), collapse=":")
                nm <- paste(c(f_label, m_label), collapse = ":")
                L_group[[i]] <- L_list[[nm]]
                group_names[i] <- nm
            }
            
            if (length(L_group) > 1) {
              L_sme <- as.matrix(do.call(rbind, lapply(2:length(L_group), function(i) L_group[[i]] - L_group[[1]])))
              V_match_mat <- as.matrix(V_match)
              
              est_sme <- as.numeric(L_sme %*% beta_match)
              V_sme <- L_sme %*% V_match_mat %*% t(L_sme)
              W_val <- as.numeric(t(est_sme) %*% MASS::ginv(V_sme) %*% est_sme)
              num_df <- nrow(L_sme)
              den_df <- self$.get_lsmeans_df(specs)
              f_val <- W_val / num_df
              p_sme <- pf(f_val, num_df, den_df, lower.tail = FALSE)
            } else {
              f_val <- NA; num_df <- 0; den_df <- self$.get_lsmeans_df(specs); p_sme <- NA
            }
            
            sme_row <- data.frame(
               Source = paste("Simple Main Effect of", paste(specs, collapse=":")),
               `F value` = f_val, num_df = num_df, den_df = den_df, Pr = p_sme,
               check.names = FALSE
            )
            
            group_contrasts <- list()
            for (i in 1:(nrow(focal_grid)-1)) {
              for (j in (i+1):nrow(focal_grid)) {
                L_diff <- L_group[[i]] - L_group[[j]]
                contrast_label <- paste0("(", group_names[i], ") - (", group_names[j], ")")
                group_contrasts[[contrast_label]] <- self$.calc_contrast(L_diff, specs)
              }
            }
            pairwise_df <- do.call(rbind, group_contrasts)
            if (adjust != "none") pairwise_df$Pr <- p.adjust(pairwise_df$Pr, method = adjust)
            if (protect && !is.na(p_sme) && p_sme > 0.05) pairwise_df$Pr <- 1.0
            
            results_by_group[[m_label]] <- list(sme = sme_row, pairwise = pairwise_df)
          }
        } else {
           group_contrasts <- list()
           for (i in 1:(length(unique_groups_full)-1)) {
             for (j in (i+1):length(unique_groups_full)) {
               name_i <- unique_groups_full[i]
               name_j <- unique_groups_full[j]
               L_diff <- L_list[[name_i]] - L_list[[name_j]]
               contrast_label <- paste0("(", name_i, ") - (", name_j, ")")
               group_contrasts[[contrast_label]] <- self$.calc_contrast(L_diff, specs)
             }
           }
           pairwise_df <- do.call(rbind, group_contrasts)
           if (adjust != "none") pairwise_df$Pr <- p.adjust(pairwise_df$Pr, method = adjust)
           results_by_group[["All"]] <- list(pairwise = pairwise_df)
        }

        class(results_by_group) <- c("rtmb_lsmeans_grouped", "list")
        attr(results_by_group, "adjustment") <- adjust
        attr(results_by_group, "protect") <- protect
        return(results_by_group)
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
      data.frame(estimate = est, `Std. Error` = se, df = df_val, `t value` = t_val, Pr = p_val, check.names = FALSE)
    },

    #' @description (Internal) Get representative DF for lsmeans.
    #' @param specs Variable names for lsmeans.
    #' @return Degrees of freedom.
    .get_lsmeans_df = function(specs) {
      df_val <- Inf
      if (is.data.frame(self$fit) && !is.null(self$fit$df)) {
        match_pattern <- paste0("^b\\[(", paste(specs, collapse="|"), ")")
        match_idx <- grepl(match_pattern, rownames(self$fit))
        if (any(match_idx)) df_val <- min(self$fit$df[match_idx], na.rm = TRUE)
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
      if (is.data.frame(fit) && "Estimate" %in% names(fit)) {
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
  ),
  private = list(
    .rbind_result_dfs = function(dfs) {
      dfs <- dfs[vapply(dfs, is.data.frame, logical(1))]
      if (length(dfs) == 0L) return(data.frame())

      all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
      dfs <- lapply(dfs, function(df) {
        missing_cols <- setdiff(all_cols, names(df))
        for (col in missing_cols) df[[col]] <- NA
        df[, all_cols, drop = FALSE]
      })

      do.call(rbind, unname(dfs))
    },
    .update_fixed_effects_with_vcov = function(V) {
      se <- sqrt(diag(V))

      for (nm in names(se)) {
        if (!nm %in% rownames(self$fit)) next

        self$fit[nm, "Std. Error"] <- se[[nm]]

        stat_col <- if ("z value" %in% names(self$fit)) "z value" else "t value"
        self$fit[nm, stat_col] <- self$fit[nm, "Estimate"] / pmax(se[[nm]], 1e-12)

        if ("df" %in% names(self$fit) && is.finite(self$fit[nm, "df"])) {
          crit <- qt(0.975, df = self$fit[nm, "df"])
          self$fit[nm, "Pr"] <- 2 * pt(-abs(self$fit[nm, stat_col]), df = self$fit[nm, "df"])
        } else {
          crit <- qnorm(0.975)
          self$fit[nm, "Pr"] <- 2 * pnorm(-abs(self$fit[nm, stat_col]))
        }

        self$fit[nm, "Lower 95%"] <- self$fit[nm, "Estimate"] - crit * se[[nm]]
        self$fit[nm, "Upper 95%"] <- self$fit[nm, "Estimate"] + crit * se[[nm]]
      }

      self$df_fixed <- self$fit
      invisible(self)
    },
    .compute_lm_robust_se = function(cluster = NULL, type = c("HC3", "HC0", "HC1", "CR1", "CR0"), update = TRUE) {
      type <- match.arg(type)
      if (is.null(cluster) && type %in% c("CR1", "CR0")) type <- "HC3"
      if (!is.null(cluster) && type %in% c("HC3", "HC0", "HC1")) type <- "CR1"
      self$se_method <- if (is.null(cluster)) "robust" else "cluster-robust"
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
        omega <- res^2
        if (type == "HC1") {
          n <- nrow(X_subset); k <- ncol(X_subset)
          omega <- omega * (n / pmax(n - k, 1))
        } else if (type == "HC3") {
          omega <- omega / pmax(1 - h, 1e-8)^2
        }
        meat <- t(X_subset) %*% (as.numeric(omega) * X_subset)
        V_final <- bread_ols %*% meat %*% bread_ols
      } else {
        if (!cluster %in% names(dat)) stop(paste("Cluster variable", cluster, "not found."))
        grp <- as.factor(dat[[cluster]])
        scores <- X_subset * as.numeric(res)
        cluster_scores <- aggregate(scores, by = list(grp), sum)[, -1, drop = FALSE]
        meat <- t(as.matrix(cluster_scores)) %*% as.matrix(cluster_scores)
        m <- length(unique(grp)); n <- nrow(X_subset); k <- ncol(X_subset)
        adj <- if (type == "CR1") (m / (m - 1)) * ((n - 1) / (n - k)) else 1
        V_final <- adj * bread_ols %*% meat %*% bread_ols
      }

      self$vcov <- V_final
      colnames(self$vcov) <- rownames(self$vcov) <- fe_names_in_fit

      self$se_method <- if (is.null(cluster)) "robust" else "cluster-robust"
      self$ci_method <- self$se_method
      self$cluster <- cluster
      self$robust_type <- type

      if (update) private$.update_fixed_effects_with_vcov(self$vcov)
      
      message(sprintf("%s standard errors calculated.", tools::toTitleCase(self$se_method)))
      return(self)
    },
      .compute_glm_robust_se = function(cluster = NULL, type = c("HC3", "HC0", "HC1", "CR1", "CR0"), update = TRUE) {
        type <- match.arg(type)
        if (is.null(cluster) && type %in% c("CR1", "CR0")) type <- "HC3"
        if (!is.null(cluster) && type %in% c("HC3", "HC0", "HC1")) type <- "CR1"

        formula <- nobars(self$model$formula)
        dat <- if (!is.null(self$model$raw_data)) {
          as.data.frame(self$model$raw_data)
        } else {
          as.data.frame(self$model$data)
        }
        mf <- stats::model.frame(formula, dat)
        X <- stats::model.matrix(formula, mf)
        y_raw <- stats::model.response(mf)

        beta_all <- if (is.data.frame(self$fit)) {
          est <- self$fit$Estimate
          names(est) <- rownames(self$fit)
          est
        } else {
          stats::coef(self$fit)
        }
        fe_names_in_fit <- names(beta_all)[grepl("^(Intercept|Intercept_c|b\\[)", names(beta_all))]
        if (length(fe_names_in_fit) == 0L) stop("No fixed effects found in fit object.", call. = FALSE)

        X_cols <- colnames(X)
        idx_map <- integer(length(fe_names_in_fit))
        for (i in seq_along(fe_names_in_fit)) {
          fname <- fe_names_in_fit[i]
          pos <- which(X_cols == fname)
          if (length(pos) == 0L) {
            if (tolower(fname) == "intercept" || fname == "intercept_c" || fname == "(intercept)") {
              pos <- which(tolower(X_cols) == "intercept" | tolower(X_cols) == "(intercept)")
            } else {
              vname <- gsub("^b\\[(.*)\\]$", "\\1", fname)
              pos <- which(X_cols == vname)
            }
          }
          if (length(pos) > 0L) idx_map[i] <- pos[1]
        }

        keep <- idx_map > 0
        if (!any(keep)) stop("None of the fixed effects could be matched.", call. = FALSE)
        fe_names_in_fit <- fe_names_in_fit[keep]
        beta <- beta_all[fe_names_in_fit]
        idx_map <- idx_map[keep]
        X_subset <- X[, idx_map, drop = FALSE]

        family <- self$model$family %||% "gaussian"
        family <- if (inherits(family, "family")) family$family else as.character(family)[1]
        if (!family %in% c("binomial", "bernoulli", "poisson", "gamma")) {
          stop(
            "robust_se() for glm currently supports binomial/bernoulli, poisson, and gamma families.",
            call. = FALSE
          )
        }

        prior_weights <- stats::model.weights(mf)
        if (is.null(prior_weights)) prior_weights <- rep(1, nrow(X_subset))

        if (family %in% c("binomial", "bernoulli")) {
          if (is.matrix(y_raw) && ncol(y_raw) == 2L) {
            binom_n <- rowSums(y_raw)
            y <- y_raw[, 1] / pmax(binom_n, 1)
            prior_weights <- prior_weights * binom_n
          } else {
            y <- as.numeric(y_raw)
            if (is.factor(y_raw)) y <- as.numeric(y_raw) - 1
          }
        } else {
          y <- as.numeric(y_raw)
        }

        eta <- as.numeric(X_subset %*% beta)
        if (family %in% c("binomial", "bernoulli")) {
          mu <- stats::plogis(eta)
          dmu_deta <- mu * (1 - mu)
          var_mu <- mu * (1 - mu)
        } else if (family == "poisson") {
          mu <- exp(eta)
          dmu_deta <- mu
          var_mu <- mu
        } else {
          mu <- exp(eta)
          dmu_deta <- mu
          var_mu <- mu^2
        }
        var_mu <- pmax(var_mu, 1e-12)

        W <- as.numeric(prior_weights) * (dmu_deta^2) / var_mu
        bread <- tryCatch(solve(t(X_subset) %*% (W * X_subset)), error = function(e) MASS::ginv(t(X_subset) %*% (W * X_subset)))
        score <- X_subset * as.numeric(prior_weights * (y - mu) * dmu_deta / var_mu)

        if (is.null(cluster)) {
          if (type %in% c("HC1", "HC3")) {
            n <- nrow(X_subset); k <- ncol(X_subset)
            if (type == "HC1") {
              score <- score * sqrt(n / pmax(n - k, 1))
            } else {
              Xw <- X_subset * sqrt(W)
              h <- tryCatch(diag(Xw %*% bread %*% t(Xw)), error = function(e) rep(0, n))
              score <- score / pmax(1 - h, 1e-8)
            }
          }
          meat <- t(score) %*% score
          V_final <- bread %*% meat %*% bread
        } else {
          if (!cluster %in% names(dat)) stop(paste("Cluster variable", cluster, "not found."), call. = FALSE)
          grp <- as.factor(dat[[cluster]])
          cluster_scores <- stats::aggregate(score, by = list(grp), sum)[, -1, drop = FALSE]
          meat <- t(as.matrix(cluster_scores)) %*% as.matrix(cluster_scores)
          m <- length(unique(grp)); n <- nrow(X_subset); k <- ncol(X_subset)
          adj <- if (type == "CR1") (m / (m - 1)) * ((n - 1) / pmax(n - k, 1)) else 1
          V_final <- adj * bread %*% meat %*% bread
        }

        self$vcov <- V_final
        colnames(self$vcov) <- rownames(self$vcov) <- fe_names_in_fit
        self$se_method <- if (is.null(cluster)) "robust" else "cluster-robust"
        self$ci_method <- self$se_method
        self$cluster <- cluster
        self$robust_type <- type

        if (update) private$.update_fixed_effects_with_vcov(self$vcov)

        message(sprintf("%s standard errors calculated.", tools::toTitleCase(self$se_method)))
        return(self)
      },
      .compute_glmer_robust_se = function(cluster = NULL, type = c("CR1", "CR0"), update = TRUE) {
        type <- match.arg(type)
        info <- self$model$extra$glmer_info
        
        if (!is.null(self$map) && length(self$map) > 0L) {
          stop("robust_se() for mixed models with mapped/fixed parameters is not supported yet.", call. = FALSE)
        }
        
        if (info$num_bars != 1L) {
          stop("robust_se() for mixed models currently supports exactly one grouping factor.", call. = FALSE)
        }

        gvar <- info$group_labels[[1]]
        if (is.null(cluster)) cluster <- gvar
        if (!identical(cluster, gvar)) {
          stop(sprintf("robust_se() for mixed models currently requires `cluster` to be the random-effect grouping variable '%s'.", gvar), call. = FALSE)
        }

        raw_data <- if (!is.null(self$model$raw_data)) self$model$raw_data else self$model$data
        if (!cluster %in% names(raw_data)) stop(paste("Cluster variable", cluster, "not found."))
        
        grp_factor <- as.factor(raw_data[[cluster]])
        clusters <- levels(grp_factor)
        m <- length(clusters)
        n <- length(grp_factor)

        # Verify one-to-one mapping
        idx_name <- info$group_idx_names[1]
        re_grp_idx <- self$model$data[[idx_name]]
        tab <- table(re_grp_idx, grp_factor)
        if (any(rowSums(tab > 0) != 1L)) {
          stop("Each random-effect level must belong to exactly one robust cluster.", call. = FALSE)
        }

        obs_to_cluster <- as.integer(grp_factor)
        re_level_to_cluster <- as.integer(apply(tab > 0, 1, function(x) which(x > 0)[1]))

        # 3. Build AD object for score calculation
        message("Preparing parameters for robust AD object...")
        
        tmp_par_list <- self$model$par_list
        tmp_par_list[[info$obs_weight_name]] <- list(
          type = "vector", length = n, dim = n, unc_length = n,
          bounds = "none", random = FALSE
        )
        for (b in 1:info$num_bars) {
          L_w <- info$num_groups_list[[b]]
          tmp_par_list[[info$weight_names[b]]] <- list(
            type = "vector", length = L_w, dim = L_w, unc_length = L_w,
            bounds = "none", random = FALSE
          )
        }

        tmp_data <- self$model$data
        tmp_data[[info$obs_weight_name]] <- NULL
        for (wn in info$weight_names) tmp_data[[wn]] <- NULL

        # Initial unconstrained parameters
        init_unc_list <- self$model$unconstrained_vector_to_list(self$par_unc)
        init_unc_list[[info$obs_weight_name]] <- rep(1, n)
        for (b in 1:info$num_bars) {
          init_unc_list[[info$weight_names[b]]] <- rep(1, info$num_groups_list[[b]])
        }

        random_names <- names(self$model$par_list)[vapply(self$model$par_list, function(p) isTRUE(p$random), logical(1))]
        
        f_ad_score <- function(y_unc_list) {
          para <- to_constrained(y_unc_list, tmp_par_list)
          -self$model$log_prob(dat = tmp_data, par = para)
        }

        obj_score <- RTMB::MakeADFun(
          func = f_ad_score,
          parameters = init_unc_list,
          random = random_names,
          map = self$map,
          intern = TRUE,
          silent = TRUE
        )
        message("AD object built. Calculating Hessian...")
        H <- obj_score$he(obj_score$par)
        message("Hessian calculated.")
        
        fe_pars <- intersect(c("Intercept", "Intercept_c", "b"), names(self$model$par_list))
        k_p <- 0
        for (nm in fe_pars) k_p <- k_p + self$model$par_list[[nm]]$unc_length
        
        n_weights <- n + sum(unlist(info$num_groups_list))
        n_total <- length(obj_score$par)
        n_theta <- n_total - n_weights
        theta_indices <- 1:n_theta
        
        meat <- matrix(0, n_theta, n_theta)
        
        current_pos <- n_theta + 1
        obs_weight_indices <- current_pos:(current_pos + n - 1)
        current_pos <- current_pos + n
        re_weight_indices <- list()
        for (b in 1:info$num_bars) {
          n_lev <- info$num_groups_list[[b]]
          re_weight_indices[[b]] <- current_pos:(current_pos + n_lev - 1)
          current_pos <- current_pos + n_lev
        }
        
        for (i in 1:m) {
          s_i <- rep(0, n_theta)
          obs_idx <- which(obs_to_cluster == i)
          if (length(obs_idx) > 0) {
            s_i <- s_i + rowSums(H[theta_indices, obs_weight_indices[obs_idx], drop = FALSE])
          }
          for (b in 1:info$num_bars) {
            lev_idx <- which(re_level_to_cluster == i)
            if (length(lev_idx) > 0) {
              s_i <- s_i + rowSums(H[theta_indices, re_weight_indices[[b]][lev_idx], drop = FALSE])
            }
          }
          meat <- meat + s_i %*% t(s_i)
        }
        
        H_theta <- H[theta_indices, theta_indices, drop = FALSE]
        H_inv <- tryCatch(solve(H_theta), error = function(e) MASS::ginv(H_theta))
        
        # 4. Final Calculations and Object Updates
        if (type == "CR1") {
          if (m <= 1L) {
            stop("CR1 robust SE requires at least two clusters.", call. = FALSE)
          }
          if (n <= k_p) {
            warning("CR1 finite-sample correction skipped because n <= p.", call. = FALSE)
            adj <- m / (m - 1)
          } else {
            adj <- (m / (m - 1)) * ((n - 1) / (n - k_p))
          }
        } else {
          adj <- 1
        }
        
        V_unc <- adj * H_inv %*% meat %*% H_inv
        
        # Correctly identify fixed-effect indices from AD object names
        par_names <- names(obj_score$par)
        fe_idx <- which(par_names %in% c("Intercept", "Intercept_c", "b"))
        
        # Standardized labels for summary alignment
        fe_labels <- unlist(lapply(fe_pars, function(nm) {
          generate_flat_names(
            base_name = nm,
            dims = self$model$par_list[[nm]]$dim,
            names_def = self$model$par_names[[nm]]
          )
        }), use.names = FALSE)
        
        if (length(fe_idx) != length(fe_labels)) {
          stop("Failed to align fixed-effect indices with displayed coefficient names.", call. = FALSE)
        }
        
        V_beta <- V_unc[fe_idx, fe_idx, drop = FALSE]
        rownames(V_beta) <- colnames(V_beta) <- fe_labels
        
        # Update object metadata
        self$vcov <- V_beta
        self$se_method <- "cluster-robust"
        self$ci_method <- "cluster-robust"
        self$cluster <- cluster
        self$robust_type <- type
        
        if (update) private$.update_fixed_effects_with_vcov(V_beta)
        
        message(sprintf("%s standard errors calculated.", tools::toTitleCase(self$se_method)))
        return(self)
      }
    )
  )

#' @export
summary.Classic_Fit <- function(object, ...) {
  object$summary(...)
}

#' @export
print.summary_Classic_Fit <- function(x, ...) {
  digits <- if (!is.null(x$digits)) x$digits else 5
  if (!is.null(x$type) && x$type == "table") { cat("\nContingency Table Analysis\n")
  } else if (!is.null(x$type) && x$type == "loglinear") { cat("\nLog-Linear Model Analysis\n")
  } else { cat("\nCall:\n"); type_label <- if (is.null(x$type)) "Generic Model" else x$type; cat(paste("Classical estimation via", type_label), "\n") }
  if (!is.null(x$se_method) && x$se_method != "wald") {
    method_label <- switch(
      x$se_method,
      robust = if (!is.null(x$cluster)) {
        paste("Robust (Cluster:", x$cluster, ")")
      } else {
        paste0("Robust (", x$robust_type %||% "HC3", ")")
      },
      `cluster-robust` = paste0(
        "Cluster-robust",
        if (!is.null(x$cluster)) paste0(" (cluster = ", x$cluster, ")") else "",
        if (!is.null(x$robust_type)) paste0(", ", x$robust_type) else ""
      ),
      bootstrap = paste0("Bootstrap (", x$bootstrap, " samples)"),
      "Standard"
    )
    cat(sprintf("Standard Errors: %s\n", method_label))
  }
  if (!is.null(x$logLik) && !is.na(x$logLik)) cat(sprintf("\nLog-Likelihood: %.3f, AIC: %.3f, BIC: %.3f\n", as.numeric(x$logLik), x$AIC, x$BIC))
  if (!is.null(x$coefficients)) {
    if (!is.null(x$type) && x$type == "table") {
      cat("\n---\n"); obs_tab <- x$extra$tab; N_tot <- sum(obs_tab); R_dim <- nrow(obs_tab); C_dim <- ncol(obs_tab)
      df_mu <- x$coefficients[grepl("^mu(\\[|$)", rownames(x$coefficients)), , drop = FALSE]
      if (nrow(df_mu) == R_dim * C_dim) {
        p_est <- df_mu$Estimate / N_tot; p_mat <- matrix(p_est, nrow = R_dim, ncol = C_dim); p_row <- rowSums(p_mat); p_col <- colSums(p_mat)
        E_mat <- matrix(0, nrow = R_dim, ncol = C_dim); for (i in 1:R_dim) for (j in 1:C_dim) E_mat[i, j] <- p_row[i] * p_col[j] * N_tot
        cat("Cell Probabilities (p) and Confidence Intervals:\n"); df_p <- df_mu
         # Correct order: Calculate then Round
         df_p$Estimate <- round(df_mu$Estimate / N_tot, digits)
         df_p$`Std. Error` <- round(df_mu$`Std. Error` / N_tot, digits)
         if ("Lower 95%" %in% names(df_p)) { 
           df_p$`Lower 95%` <- round(df_mu$`Lower 95%` / N_tot, digits)
           df_p$`Upper 95%` <- round(df_mu$`Upper 95%` / N_tot, digits)
         }
         rownames(df_p) <- gsub("^mu\\[", "p[", rownames(df_p)); cols_to_show <- intersect(c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "z value", "t value", "Pr"), names(df_p))
         print(df_p[, cols_to_show, drop = FALSE], quote = FALSE, right = TRUE)
         cat("\nExpected Counts (Independence) and Pearson Residuals:\n"); E_vec <- as.vector(E_mat); obs_vec <- as.vector(obs_tab); residuals <- (obs_vec - E_vec) / pmax(sqrt(E_vec), 1e-8)
         df_resid <- data.frame(Expected = round(E_vec, digits), Residual = round(residuals, digits), row.names = gsub("^mu\\[", "E\\[", rownames(df_mu)))
        print(df_resid, quote = FALSE, right = TRUE)
      } else {
        # Fallback for classic mode or models without mu parameters
        cat("\nTest Results:\n")
        print(x$coefficients, quote = FALSE, right = TRUE)
      }
    } else if (!is.null(x$type) && x$type == "loglinear") {
      cat("\nLog-linear Parameters and Confidence Intervals:\n"); df_params <- x$coefficients[!grepl("^mu\\[", rownames(x$coefficients)), , drop = FALSE]
      print(df_params[, intersect(c("Estimate", "Std. Error", "Lower 95%", "Upper 95%", "z value", "t value", "Pr"), names(df_params)), drop = FALSE], quote = FALSE, right = TRUE)
      df_mu <- x$coefficients[grepl("^mu\\[", rownames(x$coefficients)), , drop = FALSE]
      if (nrow(df_mu) > 0 && !is.null(x$extra$obs_Y)) {
        cat("\n---\nExpected Counts (Fitted) and Pearson Residuals:\n"); E_vec <- df_mu$Estimate; obs_vec <- x$extra$obs_Y; residuals <- (obs_vec - E_vec) / pmax(sqrt(E_vec), 1e-8)
        df_resid <- data.frame(Observed = obs_vec, Expected = round(E_vec, digits), Residual = round(residuals, digits), row.names = rownames(df_mu))
        max_r <- if (!is.null(x$max_rows)) x$max_rows else 20
        if (nrow(df_resid) > max_r) { print(df_resid[1:max_r, , drop=FALSE], quote=FALSE, right=TRUE); cat(sprintf("... (omitted %d rows)\n", nrow(df_resid) - max_r))
        } else { print(df_resid, quote = FALSE, right = TRUE) }
      }
    } else {
      cat("\nPoint Estimates and Confidence Intervals:\n"); df_to_print <- x$coefficients
      if (isTRUE(x$truncated)) df_to_print <- df_to_print[1:x$max_rows, , drop = FALSE]
      print(df_to_print, quote = FALSE, right = TRUE)
      if (isTRUE(x$truncated)) cat(sprintf("... (omitted %d parameters; use summary(max_rows = ...) to show more)\n", x$total_rows - x$max_rows))
    }
    if (!is.null(x$lm_info)) {
      cat("\n---\n"); info <- x$lm_info
      if (!is.null(info$dispersion)) { cat(sprintf("Dispersion parameter for %s family taken to be %s\n", x$family, format(round(info$dispersion, digits), nsmall = digits)))
        cat(sprintf("Null deviance: %s on %d degrees of freedom\n", format(round(info$null_deviance, 2), nsmall = 2), info$df_null))
        cat(sprintf("Residual deviance: %s on %d degrees of freedom\n", format(round(info$deviance, 2), nsmall = 2), info$df_residual))
      } else if (!is.null(info$sigma)) { cat(sprintf("Residual standard error: %s on %d degrees of freedom\n", format(round(info$sigma, digits), nsmall = digits), info$df_residual))
        cat(sprintf("Multiple R-squared: %s, Adjusted R-squared: %s\n", format(round(info$r_squared, 4), nsmall = 4), format(round(info$adj_r_squared, 4), nsmall = 4)))
        if (!is.null(info$fstatistic)) { f <- info$fstatistic; p_f <- stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
          cat(sprintf("F-statistic: %s on %d and %d DF, p-value: %s\n", format(round(f[1], 2), nsmall = 2), f[2], f[3], if (p_f < 0.001) "< .001" else format(round(p_f, 4), nsmall = 4)))
        }
      }
    }
  } else if (!is.null(x$fit_summary)) { print(x$fit_summary) }
  invisible(x)
}

#' @export
print.anova_rtmb_table <- function(x, ...) { cat("\nContingency Table Analysis Tests\n- \n"); if (!is.null(x$chisq)) print(x$chisq, row.names = TRUE, right = TRUE)
  if (!is.null(x$fisher) && !inherits(x$fisher, "try-error")) { cat("\nFisher's Exact Test for Count Data\np-value = ", sprintf("%.5f", x$fisher$p.value), "\n") }
  invisible(x)
}

#' @export
anova.Classic_Fit <- function(object, ...) {
  dots <- list(...)
  is_fit <- vapply(dots, inherits, logical(1), what = "Classic_Fit")

  if (length(dots) > 0L && any(is_fit)) {
    if (!all(is_fit)) {
      stop("When comparing models, all arguments after `object` must be Classic_Fit objects.", call. = FALSE)
    }

    call <- match.call(expand.dots = FALSE)
    dot_exprs <- as.list(call$...)
    labels <- c(deparse1(call$object), vapply(dot_exprs, deparse1, character(1)))
    fits <- c(list(object), dots)
    if (all(vapply(fits, .is_classic_lm_fit, logical(1)))) {
      return(.anova_classic_lm_f(fits, labels = labels))
    }
    return(.anova_classic_lrt(fits, labels = labels))
  }

  object$anova(...)
}

.is_classic_lm_fit <- function(x) {
  inherits(x, "Classic_Fit") &&
    !is.null(x$model) &&
    identical(x$model$type, "lm") &&
    is.finite(x$rss) &&
    is.finite(x$df_residual)
}

.anova_classic_lm_f <- function(fits, labels = NULL) {
  if (length(fits) < 2L) {
    stop("At least two Classic_Fit objects are required for model comparison.", call. = FALSE)
  }
  if (is.null(labels) || length(labels) != length(fits)) {
    labels <- paste0("Model ", seq_along(fits))
  }

  ll <- lapply(fits, stats::logLik)
  loglik <- vapply(ll, as.numeric, numeric(1))
  df <- vapply(ll, function(x) as.numeric(attr(x, "df"))[1], numeric(1))
  nobs <- vapply(ll, function(x) as.numeric(attr(x, "nobs"))[1], numeric(1))
  rss <- vapply(fits, function(x) as.numeric(x$rss)[1], numeric(1))
  res_df <- vapply(fits, function(x) as.numeric(x$df_residual)[1], numeric(1))

  if (length(unique(nobs[!is.na(nobs)])) > 1L) {
    warning("Models appear to have different numbers of observations; F tests may not be valid.", call. = FALSE)
  }

  ord <- order(df, rss, decreasing = c(FALSE, TRUE))
  labels <- labels[ord]
  fits <- fits[ord]
  loglik <- loglik[ord]
  df <- df[ord]
  rss <- rss[ord]
  res_df <- res_df[ord]

  df_diff <- c(NA_real_, diff(df))
  rss_diff <- c(NA_real_, -diff(rss))
  f_val <- rep(NA_real_, length(fits))
  p_val <- rep(NA_real_, length(fits))
  scale <- rss[length(rss)] / res_df[length(res_df)]
  df2 <- res_df[length(res_df)]
  ok <- !is.na(df_diff) & df_diff > 0 & !is.na(rss_diff) & rss_diff >= 0 &
    is.finite(scale) & scale > 0 & is.finite(df2) & df2 > 0
  f_val[ok] <- (rss_diff[ok] / df_diff[ok]) / scale
  p_val[ok] <- stats::pf(f_val[ok], df1 = df_diff[ok], df2 = df2, lower.tail = FALSE)

  res <- data.frame(
    Model = labels,
    Df = df,
    Res.Df = res_df,
    RSS = rss,
    `Df Diff` = df_diff,
    `Sum of Sq` = rss_diff,
    F = f_val,
    `Pr(>F)` = p_val,
    logLik = loglik,
    AIC = vapply(fits, stats::AIC, numeric(1)),
    BIC = vapply(fits, stats::BIC, numeric(1)),
    check.names = FALSE
  )
  rownames(res) <- seq_len(nrow(res))
  class(res) <- c("anova_classic_lm_f", "anova", "data.frame")
  attr(res, "heading") <- "Analysis of Variance Table (F tests)"
  res
}

.anova_classic_lrt <- function(fits, labels = NULL) {
  if (length(fits) < 2L) {
    stop("At least two Classic_Fit objects are required for model comparison.", call. = FALSE)
  }

  if (is.null(labels) || length(labels) != length(fits)) {
    labels <- paste0("Model ", seq_along(fits))
  }

  ll <- lapply(fits, stats::logLik)
  loglik <- vapply(ll, as.numeric, numeric(1))
  df <- vapply(ll, function(x) {
    val <- attr(x, "df")
    if (is.null(val) || length(val) == 0L) NA_real_ else as.numeric(val)[1]
  }, numeric(1))
  nobs <- vapply(ll, function(x) {
    val <- attr(x, "nobs")
    if (is.null(val) || length(val) == 0L) NA_real_ else as.numeric(val)[1]
  }, numeric(1))

  if (anyNA(loglik)) {
    stop("Likelihood-ratio comparison requires all models to have a finite logLik().", call. = FALSE)
  }
  if (anyNA(df)) {
    stop("Likelihood-ratio comparison requires all logLik() values to have a valid df attribute.", call. = FALSE)
  }

  known_nobs <- nobs[!is.na(nobs)]
  if (length(unique(known_nobs)) > 1L) {
    warning("Models appear to have different numbers of observations; likelihood-ratio tests may not be valid.", call. = FALSE)
  }

  ord <- order(df, loglik)
  labels <- labels[ord]
  loglik <- loglik[ord]
  df <- df[ord]
  fits <- fits[ord]

  aic <- vapply(fits, stats::AIC, numeric(1))
  bic <- vapply(fits, stats::BIC, numeric(1))
  chisq <- c(NA_real_, 2 * diff(loglik))
  df_diff <- c(NA_real_, diff(df))
  p_val <- rep(NA_real_, length(fits))
  ok <- !is.na(chisq) & !is.na(df_diff) & df_diff > 0
  p_val[ok] <- stats::pchisq(chisq[ok], df = df_diff[ok], lower.tail = FALSE)

  res <- data.frame(
    Model = labels,
    Df = df,
    logLik = loglik,
    AIC = aic,
    BIC = bic,
    Chisq = chisq,
    `Chi Df` = df_diff,
    `Pr(>Chisq)` = p_val,
    check.names = FALSE
  )
  rownames(res) <- seq_len(nrow(res))
  class(res) <- c("anova_classic_lrt", "anova", "data.frame")
  attr(res, "heading") <- "Likelihood Ratio Tests"
  res
}

#' @export
print.anova_classic_lrt <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  heading <- attr(x, "heading")
  if (!is.null(heading)) cat(heading, "\n", sep = "")
  x_print <- x
  class(x_print) <- "data.frame"
  print(x_print, digits = digits, na.print = "", row.names = FALSE, ...)
  invisible(x)
}

#' @export
print.anova_classic_lm_f <- print.anova_classic_lrt

#' @export
print.Classic_Fit <- function(x, ...) x$print(...)

#' Least Squares Means (Marginal Means)
#'
#' @description
#' Calculate least squares means (also known as marginal means or predicted means)
#' for categorical factors in a fitted model. It also supports pairwise comparisons
#' and simple main effects.
#'
#' @param object A fitted model object (e.g., `Classic_Fit`).
#' @param specs A character vector of factor names to calculate means for.
#' @param ... Additional arguments passed to the method.
#' @export
lsmeans <- function(object, specs, ...) UseMethod("lsmeans")

#' @export
lsmeans.Classic_Fit <- function(object, specs, ...) object$lsmeans(specs, ...)

#' @export
print.rtmb_lsmeans <- function(x, ...) {
  fmt_df <- function(v) sapply(round(v, 1), function(z) sprintf("%g", z))
  fmt_5 <- function(v) sprintf("%.5f", v)
  
  disp <- as.data.frame(x)
  
  cols_5 <- intersect(names(disp), c("estimate", "Std. Error", "Lower 95%", "Upper 95%"))
  for (col in cols_5) disp[[col]] <- fmt_5(disp[[col]])
  
  if ("df" %in% names(disp)) disp$df <- fmt_df(disp$df)
  
  print(disp, ...)
  invisible(x)
}

#' @export
print.rtmb_lsmeans_grouped <- function(x, ...) {
  adj <- attr(x, "adjustment"); prot <- attr(x, "protect"); cat("\nPost-hoc Analysis: Simple Main Effects and Pairwise Comparisons\n")
  cat(sprintf("P-value adjustment: %s (applied within groups)\n", adj)); if (prot) cat("Hierarchical testing (Protected): Pairwise results only shown if SME is significant.\n")
  fmt_df <- function(v) sapply(round(v, 1), function(z) sprintf("%g", z)); fmt_5 <- function(v) sprintf("%.5f", v); fmt_p <- function(v) sub("^0", "", sprintf("%.5f", v))
  for (grp in names(x)) { cat("\n", rep("-", nchar(grp) + 10), "\nGroup: ", grp, "\n", rep("-", nchar(grp) + 10), "\n", sep="")
    if (!is.null(x[[grp]]$sme)) { cat("\nSimple Main Effect (Wald F-test):\n"); sme_disp <- x[[grp]]$sme; sme_disp$`F value` <- fmt_5(sme_disp$`F value`); sme_disp$den_df <- fmt_df(sme_disp$den_df); sme_disp$Pr <- fmt_p(sme_disp$Pr); print(sme_disp, row.names = FALSE) }
    cat("\nPairwise Comparisons:\n"); pw_disp <- x[[grp]]$pairwise; pw_disp$estimate <- fmt_5(pw_disp$estimate); pw_disp$`Std. Error` <- fmt_5(pw_disp$`Std. Error`); pw_disp$df <- fmt_df(pw_disp$df); pw_disp$`t value` <- fmt_5(pw_disp$`t value`); pw_disp$Pr <- fmt_p(pw_disp$Pr); print(pw_disp); cat("\n")
  }
  invisible(x)
}

#' @export
logLik.Classic_Fit <- function(object, ...) object$logLik()

#' @export
AIC.Classic_Fit <- function(object, ..., k = 2) object$AIC()

#' @export
BIC.Classic_Fit <- function(object, ...) object$BIC()

#' @export
print.anova_rtmb <- function(x, ...) {
  # Use the standard anova print for the table
  # stats:::print.anova will handle the 'heading' attribute
  # Temporarily remove class to avoid infinite recursion
  x_copy <- x
  class(x_copy) <- c("anova", "data.frame")
  print(x_copy, ...)
  
  # Print LM info if available
  info <- attr(x, "lm_info")
  if (!is.null(info)) {
    cat("---\n")
    cat(sprintf("Residual standard error: %.4f on %.1f degrees of freedom\n", info$sigma, info$df_resid))
    cat(sprintf("Multiple R-squared: %.4f,  Adjusted R-squared: %.4f\n", info$r_squared, info$adj_r_squared))
    if (!is.null(info$fstatistic)) {
       cat(sprintf("F-statistic: %.4f on %d and %.1f DF,  p-value: %s\n", 
                   info$fstatistic["value"], 
                   as.integer(info$fstatistic["numdf"]), 
                   info$fstatistic["dendf"],
                   format.pval(info$p_value)))
    }
  }
  invisible(x)
}
