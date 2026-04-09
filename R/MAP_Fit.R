#' MAP fit object
#'
#' An R6 class storing optimization results from maximum a posteriori
#' (MAP) estimation.
#'
#' @param par_vec Parameter vector on the unconstrained scale.
#' @param par Parameter list on the constrained scale.
#' @param objective RTMB objective function object.
#' @param log_ml Log marginal likelihood or related model criterion,
#'   if available.
#' @param convergence Optimizer convergence code.
#' @param sd_rep Standard deviation report object.
#' @param df_fixed Summary table for fixed-effect parameters.
#' @param random_effects Random effect estimates, if available.
#' @param tran_est List of transformed parameter estimates, if available.
#' @param gq_est List of generated quantity estimates, if available.
#' @param pars Character vector specifying the names of parameters to summarize.
#' @param max_rows Maximum number of rows to print in summaries.
#' @param digits Number of digits to print.
#' @param ... Additional arguments passed to print methods.
#'
#' @field par_vec Parameter vector on the unconstrained scale.
#' @field par Parameter list on the constrained scale.
#' @field objective RTMB objective function object.
#' @field log_ml Log marginal likelihood or related model criterion.
#' @field convergence Optimizer convergence code.
#' @field sd_rep Standard deviation report object.
#' @field df_fixed Summary table for fixed-effect parameters.
#' @field random_effects Random effect estimates.
#' @field df_tran List of transformed parameter estimates.
#' @field df_gq List of generated quantity estimates.
#' @field opt_history A vector of optimize objective history
#'
MAP_Fit <- R6::R6Class(
  classname = "map_fit",

  public = list(
    # --- フィールド ---
    par_vec        = NULL,
    par            = NULL,
    objective      = NULL,
    log_ml         = NULL,
    convergence    = NULL,
    sd_rep         = NULL,
    df_fixed       = NULL,
    random_effects = NULL,
    df_tran        = NULL,
    df_gq          = NULL,
    opt_history    = NULL,

    #' @description Create a new `MAP_Fit` object.
    #' @param par_vec Parameter vector on the unconstrained scale.
    #' @param par Parameter list on the constrained scale.
    #' @param objective RTMB objective function object.
    #' @param log_ml Log marginal likelihood or related model criterion.
    #' @param convergence Optimizer convergence code.
    #' @param sd_rep Standard deviation report object.
    #' @param df_fixed Summary table for fixed-effect parameters.
    #' @param random_effects Random effect estimates.
    #' @param df_tran List of transformed parameter estimates.
    #' @param df_gq List of generated quantity estimates.
    #' @param opt_history A vector of optimize objective history
    initialize = function(par_vec, par, objective, log_ml, convergence, sd_rep, df_fixed,
                          random_effects, df_tran = NULL, df_gq = NULL, opt_history = NULL) {
      self$par_vec <- par_vec
      self$par <- par
      self$objective <- objective
      self$log_ml <- log_ml
      self$convergence <- convergence
      self$sd_rep <- sd_rep
      self$df_fixed <- df_fixed
      self$random_effects <- random_effects
      self$df_tran <- df_tran
      self$df_gq <- df_gq
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

      # 1. すべてのデータフレームを縦に結合する
      all_dfs <- list()
      if (!is.null(self$df_fixed) && nrow(self$df_fixed) > 0) all_dfs$fixed <- self$df_fixed
      if (!is.null(self$df_tran) && nrow(self$df_tran) > 0) all_dfs$tran <- self$df_tran
      if (!is.null(self$df_gq) && nrow(self$df_gq) > 0) all_dfs$gq <- self$df_gq

      if (length(all_dfs) == 0) {
        cat("\nNo parameters to display.\n")
        return(invisible(NULL))
      }

      df_combined <- do.call(rbind, unname(all_dfs))

      # 2. 変数名での絞り込み
      if (!is.null(pars)) {
        base_names <- gsub("\\[.*\\]$", "", rownames(df_combined))
        keep_idx <- rownames(df_combined) %in% pars | base_names %in% pars
        df_combined <- df_combined[keep_idx, , drop = FALSE]
      }

      if (nrow(df_combined) == 0) {
        cat("\nNo matching parameters found.\n")
        return(invisible(df_combined))
      }

      cat("\nPoint Estimates and 95% Wald CI:\n")

      # 3. 表示用のデータフレーム作成
      out_df <- data.frame(variable = rownames(df_combined), df_combined, check.names = FALSE, stringsAsFactors = FALSE)
      rownames(out_df) <- NULL

      # 4. ここで全体に対して max_rows を適用する
      if (!is.null(max_rows) && nrow(out_df) > max_rows) {
        out_df <- head(out_df, max_rows)
      }

      # 5. MCMC/VB用のprint関数に回して表を整形
      class(out_df) <- c("summary_BayesRTMB", "data.frame")
      print(out_df, digits = digits)

      cat("\n")
      invisible(df_combined)
    },

    #' @description Print a brief summary of the fitted object.
    #' @param pars Character vector specifying the names of parameters to summarize.
    #' @param max_rows Maximum number of rows to print in summaries.
    #' @param digits Number of digits to print.
    #' @return The object itself, invisibly.
    print = function(pars = NULL, max_rows = 10, digits = 2, ...) {
      self$summary(pars = pars, max_rows = max_rows, digits = digits, ...)
      invisible(self)
    }
  )
)
