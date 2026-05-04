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

    #' @description Create a new `Classic_Model` object.
    #' @param type Character string specifying the model type (e.g., "lm").
    #' @param formula The formula used for the model.
    #' @param data The data used for estimation.
    #' @param family The distribution family.
    #' @param view Parameter names to prioritize in summary display.
    initialize = function(type, formula, data, family = "gaussian", view = NULL) {
      self$type <- type
      self$formula <- formula
      self$data <- data
      self$family <- family
      self$view <- view
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
