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
#' @field tran_est List of transformed parameter estimates.
#' @field gq_est List of generated quantity estimates.
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
    tran_est       = NULL, # 追加
    gq_est         = NULL, # 追加

    #' @description Create a new `MAP_Fit` object.
    #' @param par_vec Parameter vector on the unconstrained scale.
    #' @param par Parameter list on the constrained scale.
    #' @param objective RTMB objective function object.
    #' @param log_ml Log marginal likelihood or related model criterion.
    #' @param convergence Optimizer convergence code.
    #' @param sd_rep Standard deviation report object.
    #' @param df_fixed Summary table for fixed-effect parameters.
    #' @param random_effects Random effect estimates.
    #' @param tran_est List of transformed parameter estimates.
    #' @param gq_est List of generated quantity estimates.
    initialize = function(par_vec, par, objective, log_ml, convergence, sd_rep, df_fixed, random_effects, tran_est = NULL, gq_est = NULL) {
      self$par_vec <- par_vec
      self$par <- par
      self$objective <- objective
      self$log_ml <- log_ml
      self$convergence <- convergence
      self$sd_rep <- sd_rep
      self$df_fixed <- df_fixed
      self$random_effects <- random_effects
      self$tran_est <- tran_est
      self$gq_est <- gq_est
    },

    #' @description Summarize MAP estimates.
    #' @param pars Character vector specifying the names of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param max_rows Maximum number of rows to print in summaries.
    #' @param digits Number of digits to print.
    #' @return A summary object, typically a data frame or list.
    summary = function(pars = NULL, max_rows = 20, digits = 2) {
      cat("\nCall:\nMAP Estimation via RTMB\n\n")

      cat("Fixed Effects Coefficients (Constrained Scale):\n")
      if (!is.null(self$df_fixed)) {
        df <- self$df_fixed

        # --- 変数名のフィルタリング処理を追加 ---
        if (!is.null(pars)) {
          # "[1,1]" などのインデックス部分を削除してベース名を取得
          base_names <- gsub("\\[.*\\]$", "", rownames(df))
          # 完全一致、またはベース名での一致を判定
          keep_idx <- rownames(df) %in% pars | base_names %in% pars
          df <- df[keep_idx, , drop = FALSE]
        }

        if (nrow(df) == 0) {
          cat("No matching fixed effects found.\n")
        } else {
          # 行数の制限処理
          if (!is.null(max_rows) && nrow(df) > max_rows) {
            base::print(head(df, max_rows), digits = digits)
            #cat(sprintf("... (残り %d 行を省略。max_rows で表示数を変更できます)\n", nrow(df) - max_rows))
          } else {
            base::print(df, digits = digits)
          }
        }
      } else {
        cat("No fixed effects to display.\n")
      }

      cat("---\nNegative Log-Posterior:", self$objective, "\n")

      if (!is.null(self$log_ml) && !is.na(self$log_ml)) {
        cat("Approx. Log Marginal Likelihood (Laplace):", self$log_ml, "\n")
      } else {
        cat("Approx. Log Marginal Likelihood (Laplace): NA (Calculation failed or not applicable)\n")
      }

      if (!is.null(self$random_effects)) {
        cat("Note: Random effects are stored in $random_effects\n")
      }

      print_list_items <- function(est_list, category_name) {
        # 表示対象の変数を絞り込む
        target_names <- names(est_list)
        if (!is.null(pars)) {
          target_names <- intersect(target_names, pars)
        }

        if (length(target_names) > 0) {
          cat(sprintf("\n%s (Point Estimates):\n", category_name))
          for (name in target_names) {
            cat(sprintf("$%s\n", name))
            base::print(est_list[[name]], digits = digits)
          }
        }
      }

      if (!is.null(self$tran_est)) {
        print_list_items(self$tran_est, "Transformed Parameters")
      }

      if (!is.null(self$gq_est)) {
        print_list_items(self$gq_est, "Generated Quantities")
      }

      cat("\n")
      invisible(self$df_fixed) # 戻り値は元のdf_fixed全体のままとするか、フィルタ済みのdfにするかは用途に応じて調整してください
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
