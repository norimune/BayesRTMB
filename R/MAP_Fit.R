#' MAP fit object
#'
#' An R6 class storing optimization results from maximum a posteriori
#' (MAP) estimation.
#
#' @param model The `RTMB_Model` object used for estimation.
#' @param par_vec Parameter vector on the unconstrained scale.
#' @param par Parameter list on the constrained scale.
#' @param objective RTMB objective function object.
#' @param log_ml Log marginal likelihood or related model criterion, if available.
#' @param convergence Optimizer convergence code.
#' @param sd_rep Standard deviation report object.
#' @param df_fixed Summary table for fixed-effect parameters.
#' @param random_effects Random effect estimates, if available.
#' @param df_transform Summary table for transformed parameter estimates, if available.
#' @param df_generate Summary table for generated quantity estimates, if available.
#' @param opt_history A vector of optimize objective history.
#' @param transform List of transformed parameters maintaining their original dimensions.
#' @param generate List of generated quantities maintaining their original dimensions.
#'
#' @field model The `RTMB_Model` object used for estimation.
#' @field par_vec Parameter vector on the unconstrained scale.
#' @field par Parameter list on the constrained scale.
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
#'
MAP_Fit <- R6::R6Class(
  classname = "map_fit",

  public = list(
    # --- フィールド ---
    model          = NULL,
    par_vec        = NULL,
    par            = NULL,
    objective      = NULL,
    log_ml         = NULL,
    convergence    = NULL,
    sd_rep         = NULL,
    df_fixed       = NULL,
    random_effects = NULL,
    df_transform   = NULL, # 変更
    df_generate    = NULL, # 変更
    opt_history    = NULL,
    transform      = NULL, # 追加
    generate       = NULL, # 追加

    #' @description Create a new `MAP_Fit` object.
    #' @param model The `RTMB_Model` object used for estimation.
    #' @param par_vec Parameter vector on the unconstrained scale.
    #' @param par Parameter list on the constrained scale.
    #' @param objective RTMB objective function object.
    #' @param log_ml Log marginal likelihood or related model criterion.
    #' @param convergence Optimizer convergence code.
    #' @param sd_rep Standard deviation report object.
    #' @param df_fixed Summary table for fixed-effect parameters.
    #' @param random_effects Random effect estimates.
    #' @param df_transform Summary table for transformed parameter estimates.
    #' @param df_generate Summary table for generated quantity estimates.
    #' @param opt_history A vector of optimize objective history.
    #' @param transform List of transformed parameters maintaining their original dimensions.
    #' @param generate List of generated quantities maintaining their original dimensions.
    initialize = function(model,par_vec, par, objective, log_ml, convergence, sd_rep, df_fixed,
                          random_effects, df_transform = NULL, df_generate = NULL, opt_history = NULL,
                          transform = NULL, generate = NULL) {
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

      # 複素数型（complex）の混入を防ぐための実数化処理
      if (!is.null(self$par_vec)) self$par_vec <- Re(self$par_vec)
      if (!is.null(self$par)) self$par <- lapply(self$par, Re)
      if (!is.null(self$random_effects)) self$random_effects <- lapply(self$random_effects, Re)
      if (!is.null(self$transform)) self$transform <- lapply(self$transform, Re)
      if (!is.null(self$generate)) self$generate <- lapply(self$generate, Re)
    },

    #' @description Summarize MAP estimates.
    #' @param pars Character vector specifying the names of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param max_rows Maximum number of rows to print in summaries. Default is 10.
    #' @param digits Number of digits to print.
    #' @return A summary object, typically a data frame or list.
    summary = function(pars = NULL, max_rows = 10, digits = 5) {
      cat("\nCall:\nMAP Estimation via RTMB\n")
      cat(sprintf("\nNegative Log-Posterior: %.2f\n", self$objective))

      if (!is.null(self$log_ml) && !is.na(self$log_ml)) {
        cat(sprintf("Approx. Log Marginal Likelihood (Laplace): %.2f\n", self$log_ml))
      } else {
        cat("Approx. Log Marginal Likelihood (Laplace): NA\n")
      }

      if (!is.null(self$random_effects)) {
        cat("Note: Random effects are stored in $random_effects\n")
      }

      all_dfs <- list()
      if (!is.null(self$df_fixed) && nrow(self$df_fixed) > 0) all_dfs$fixed <- self$df_fixed
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

      # --- lp と model$view による優先並び替え ---
      if (nrow(df_combined) > 0) {
        var_names <- rownames(df_combined)
        base_names <- gsub("\\[.*\\]$", "", var_names)

        # lp を常に最優先にする
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

      cat("\nPoint Estimates and 95% Wald CI:\n")

      num_cols <- sapply(df_combined, is.numeric)
      df_combined[num_cols] <- lapply(df_combined[num_cols], function(x) {
        x[abs(x) < 1e-12 & !is.na(x)] <- 0
        return(x)
      })

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
    #' @param code An `rtmb_code({ ... })` or `{ ... }` block containing the logic
    #' to be calculated using the MAP estimate.
    #' @return The `MAP_Fit` object itself (invisibly).
    #' Results are added or updated in the `generate` list.
    generated_quantities = function(code) {
      # 1. Capture and extract AST
      raw_code <- substitute(code)
      if (is.name(raw_code)) {
        evaluated <- tryCatch(eval(raw_code, envir = parent.frame()), error = function(e) NULL)
        if (is.language(evaluated) || is.call(evaluated)) code <- evaluated
      }

      gen_ast <- if (is.call(code) && identical(code[[1]], as.name("rtmb_code"))) code$generate else code

      # 2. Prepare function and environment
      gen_fn <- eval(bquote(transform_code(.(gen_ast))))
      environment(gen_fn) <- parent.env(globalenv())

      # 3. Create context (par + transform + generate)
      p_list <- self$par
      if (!is.null(self$transform)) p_list <- c(p_list, self$transform)
      if (!is.null(self$generate)) p_list <- c(p_list, self$generate)

      # 4. Execute
      res <- gen_fn(self$model$data, p_list)

      # 5. Store results
      if (is.null(self$generate)) self$generate <- list()
      for (n in names(res)) {
        self$generate[[n]] <- res[[n]]
      }

      cat("Generated quantities updated.\n")
      invisible(self)
    },
    #' @description Performs orthogonal Procrustes rotation on the MAP estimate.
    #' @param target Character string specifying the name of the matrix parameter
    #' to be used as the rotation reference.
    #' @param reference Matrix to rotate towards. If NULL, the target itself is
    #' used (no-op unless linked is provided).
    #' @param linked A character vector of other parameter names to be rotated
    #' in the same direction as the target.
    #' @return The `MAP_Fit` object itself (invisibly).
    #' Rotated values are saved to `generate` with a `_rot` suffix.
    rotate = function(target, reference = NULL, linked = NULL) {
      if (is.null(reference)) {
        warning("No reference provided for Procrustes rotation.")
        return(invisible(self))
      }

      # build AST for rotation
      exprs <- list()
      exprs[[length(exprs) + 1]] <- bquote(svd_res <- svd(t(.(as.name(target))) %*% .(reference)))
      exprs[[length(exprs) + 1]] <- quote(Q <- svd_res$u %*% t(svd_res$v))

      target_rot <- paste0(target, "_rot")
      exprs[[length(exprs) + 1]] <- bquote(.(as.name(target_rot)) <- .(as.name(target)) %*% Q)

      ret_list <- list()
      ret_list[[target_rot]] <- as.name(target_rot)

      if (!is.null(linked)) {
        for (l_var in linked) {
          l_rot <- paste0(l_var, "_rot")
          exprs[[length(exprs) + 1]] <- bquote(.(as.name(l_rot)) <- .(as.name(l_var)) %*% Q)
          ret_list[[l_rot]] <- as.name(l_rot)
        }
      }

      exprs[[length(exprs) + 1]] <- bquote(return(.(as.call(c(list(as.name("list")), ret_list)))))
      self$generated_quantities(as.call(c(list(as.name("{")), exprs)))
      invisible(self)
    },
    #' @description Rotates factor loadings and optionally rotates associated parameters.
    #' @param target Character string specifying the factor loadings matrix to
    #' base the rotation on. Defaults to "loadings".
    #' @param linked Character vector of parameters to which the same rotation
    #' matrix should be applied.
    #' @param scores Character vector of parameters to which the inverse-transpose
    #' of the rotation matrix should be applied.
    #' @param rotate Character string specifying the rotation method (e.g., "promax").
    #' @param ... Additional arguments passed to the rotation function.
    #' @return The `MAP_Fit` object itself (invisibly).
    #' Rotated values are saved to `generate` with a method-specific suffix.
    fa_rotate = function(target = "loadings", linked = NULL, scores = NULL, rotate = "promax", ...) {
      cat(sprintf("Applying %s rotation to %s...\n", rotate, target))

      # Identify function
      if (exists(rotate, mode = "function")) {
        fn_call <- as.name(rotate)
      } else {
        fn_call <- call("::", as.name("GPArotation"), as.name(rotate))
      }

      # Test rotation type to see if it's oblique (has Phi)
      t_val <- if (!is.null(self$par[[target]])) self$par[[target]] else self$transform[[target]]
      test_rot <- eval(as.call(list(fn_call, t_val)))
      is_matrix_rot <- is.matrix(test_rot)

      # Build AST
      exprs <- list()
      rot_name <- paste0(target, "_", rotate)
      exprs[[length(exprs) + 1]] <- bquote(rot_obj <- .(fn_call)(.(as.name(target))))

      ret_list <- list()
      if (is_matrix_rot) {
        exprs[[length(exprs) + 1]] <- quote(rot_mat <- unclass(rot_obj))
      } else {
        exprs[[length(exprs) + 1]] <- quote(rot_mat <- unclass(rot_obj$loadings))
        if (!is.null(test_rot$Phi)) ret_list[["fa_cor"]] <- quote(rot_obj$Phi)
      }
      ret_list[[rot_name]] <- as.name("rot_mat")

      # Handle linked and scores
      if (!is_matrix_rot && (!is.null(linked) || !is.null(scores))) {
        exprs[[length(exprs) + 1]] <- quote(rot_Th <- if (!is.null(rot_obj$Th)) rot_obj$Th else rot_obj$rotmat)

        if (!is.null(linked)) {
          for (l_var in linked) {
            l_rot <- paste0(l_var, "_", rotate)
            exprs[[length(exprs) + 1]] <- bquote(.(as.name(l_rot)) <- .(as.name(l_var)) %*% rot_Th)
            ret_list[[l_rot]] <- as.name(l_rot)
          }
        }
        if (!is.null(scores)) {
          exprs[[length(exprs) + 1]] <- quote(rot_Th_inv <- solve(t(rot_Th)))
          for (s_var in scores) {
            s_rot <- paste0(s_var, "_", rotate)
            exprs[[length(exprs) + 1]] <- bquote(.(as.name(s_rot)) <- .(as.name(s_var)) %*% rot_Th_inv)
            ret_list[[s_rot]] <- as.name(s_rot)
          }
        }
      }

      exprs[[length(exprs) + 1]] <- bquote(return(.(as.call(c(list(as.name("list")), ret_list)))))
      self$generated_quantities(as.call(c(list(as.name("{")), exprs)))
      invisible(self)
    }
  )
)
