#' VB fit object
#'
#' An R6 class storing posterior samples and related information
#' from Automatic Differentiation Variational Inference (ADVI).
#'
#' @field model An `RTMB_Model` object used for estimation.
#' @field fit A 3D array of posterior draws for fixed model parameters.
#' @field random_fit A 3D array of posterior draws for random effects, if available.
#' @field transform_fit A 3D array of posterior draws for transformed parameters, if available.
#' @field generate_fit A 3D array of posterior draws for generated quantities, if available.
#' @field transform_dims A list storing dimension information for transformed parameters.
#' @field generate_dims A list storing dimension information for generated quantities.
#' @field elbo_history A list of numeric vectors storing the Evidence Lower Bound (ELBO) history during optimization for each chain.
#' @field laplace Logical; whether Laplace approximation was used to marginalize random effects.
#' @field posterior_mean A named numeric vector of posterior mean estimates.
#' @field ELBO A numeric vector of final ELBO values for each chain.
#' @field rel_obj_vals A numeric vector of final relative objective tolerance values for each chain.
#' @field best_chain Integer; the index of the chain with the maximum ELBO.
#' @field mu_history Matrix of the parameter trajectory from the final window.
#'
VB_Fit <- R6::R6Class(
  classname = "advi_fit",
  inherit = RTMB_Fit_Base,

  public = list(
    # --- Fields ---
    model          = NULL, # Reference to RTMB_Model instance
    fit            = NULL,
    random_fit     = NULL,
    transform_fit  = NULL, # Store transformed quantities
    generate_fit   = NULL, # Store generated quantities
    transform_dims = NULL, # Store dimension info for transformed quantities
    generate_dims  = NULL, # Store dimension info for generated quantities
    elbo_history   = NULL, # Store ELBO history
    laplace        = NULL,
    posterior_mean = NULL,
    ELBO           = NULL,
    rel_obj_vals   = NULL,
    best_chain     = NULL,
    mu_history     = NULL,

    # 1. Constructor
    #' @description Get point estimate for a target parameter (internal use).
    #' @param target Target parameter name.
    #' @return Matrix or array of point estimate.
    get_point_estimate = function(target) {
      target_draws <- self$draws(pars = target, inc_transform = TRUE, inc_generate = TRUE)
      if (dim(target_draws)[3] == 0) stop("Parameter not found: ", target)

      lp_draws <- self$draws(pars = "lp", inc_transform = FALSE, inc_generate = FALSE)
      max_idx <- which(lp_draws == max(lp_draws, na.rm = TRUE), arr.ind = TRUE)
      best_iter <- max_idx[1, 1]
      best_chain <- max_idx[1, 2]

      t_info <- self$model$par_list[[target]]
      if (!is.null(t_info)) {
        t_dim <- t_info$dim
      } else if (!is.null(self$transform_dims[[target]])) {
        t_dim <- self$transform_dims[[target]]
      } else if (!is.null(self$generate_dims[[target]])) {
        t_dim <- self$generate_dims[[target]]
      } else {
        t_dim <- dim(target_draws)[3]
      }

      target_map_flat <- target_draws[best_iter, best_chain, ]
      target_map <- target_map_flat
      if (length(t_dim) > 1) dim(target_map) <- t_dim

      return(target_map)
    },

    #' @description Create a new `VB_Fit` object.
    #' @param model An `RTMB_Model` object.
    #' @param fit A 3D array of parameter draws.
    #' @param random_fit A 3D array of random effect draws.
    #' @param elbo_history A list of numeric vectors of ELBO values for each chain.
    #' @param laplace Logical; indicates if Laplace approximation was used.
    #' @param posterior_mean A named numeric vector of posterior means.
    #' @param ELBO A numeric vector of final ELBO values for each chain.
    #' @param rel_obj_vals A numeric vector of final relative objective tolerance values for each chain.
    #' @param best_chain Integer; the index of the chain with the maximum ELBO.
    #' @param mu_history Matrix of the parameter trajectory from the final window.
    initialize = function(model, fit, random_fit, elbo_history, laplace, posterior_mean, ELBO, rel_obj_vals, best_chain, mu_history) {
      self$model <- model
      self$fit <- fit
      self$random_fit <- random_fit
      self$elbo_history <- elbo_history
      self$laplace <- laplace
      self$posterior_mean <- posterior_mean
      self$ELBO <- ELBO
      self$rel_obj_vals <- rel_obj_vals
      self$best_chain <- best_chain
      self$mu_history <- mu_history
      self$transform_fit <- NULL
      self$transform_dims <- list()
      self$generate_fit <- NULL
      self$generate_dims <- list()
      # Ensure S3 dispatch works for base class methods
      class(self) <- c(class(self), "RTMB_Fit_Base")
    },

    #' @description Print a brief summary of the fitted object.
    #' @param ... Additional arguments passed to the `summary` method.
    #' @return The object itself, invisibly.
    print = function(...) {
      print(self$summary(...))
      invisible(self)
    },

    #' @description Extract posterior draws for selected parameters.
    #' @param pars Character or numeric vector specifying the names or indices of parameters to extract. If NULL, all available parameters are extracted.
    #' @param chains Numeric vector specifying the chains to extract. If NULL, draws from all chains are returned.
    #' @param best_chains Integer; number of best chains to retain based on ELBO. If supplied, chains with the highest ELBO are retained.
    #' @param inc_random Logical; whether to include random effects in the output. Default is FALSE.
    #' @param inc_transform Logical; whether to include transformed parameters in the output. Default is TRUE.
    #' @param inc_generate Logical; whether to include generated quantities in the output. Default is TRUE.
    #' @param best_only Logical; whether to extract only from the chain with the maximum ELBO. Default is FALSE unless explicitly requested.
    #' @return A 3D array of posterior draws `[iterations, chains, parameters]`.
    draws = function(pars = NULL, chains = NULL, best_chains = NULL,
                     inc_random = FALSE, inc_transform = TRUE, inc_generate = TRUE,
                     best_only = FALSE) {
      total_chains <- dim(self$fit)[2]
      available_chains <- seq_len(total_chains)

      if (!is.null(chains)) {
        available_chains <- intersect(available_chains, chains)
        if (length(available_chains) == 0) {
          stop("The specified chains were not found.", call. = FALSE)
        }
      }

      if (!is.null(best_chains)) {
        elbo_vals <- self$ELBO[available_chains]
        n_best <- min(best_chains, length(available_chains))
        ordered_idx <- order(elbo_vals, decreasing = TRUE, na.last = NA)
        available_chains <- available_chains[ordered_idx[seq_len(n_best)]]
      }

      available_chains <- sort(unique(available_chains))
      out_array <- self$fit[, available_chains, , drop = FALSE]

      if (inc_random && !is.null(self$random_fit)) {
        P1 <- dim(out_array)[3]; P2 <- dim(self$random_fit)[3]
        I <- dim(out_array)[1]; C <- dim(out_array)[2]
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$random_fit[, available_chains, , drop = FALSE]
        dimnames(new_out) <- list(
          iteration = dimnames(out_array)[[1]],
          chain = dimnames(out_array)[[2]],
          variable = c(dimnames(out_array)[[3]], dimnames(self$random_fit)[[3]])
        )
        out_array <- new_out
      }

      if (inc_transform && !is.null(self$transform_fit)) {
        P1 <- dim(out_array)[3]; P2 <- dim(self$transform_fit)[3]
        I <- dim(out_array)[1]; C <- dim(out_array)[2]
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$transform_fit[, available_chains, , drop = FALSE]
        dimnames(new_out) <- list(
          iteration = dimnames(out_array)[[1]],
          chain = dimnames(out_array)[[2]],
          variable = c(dimnames(out_array)[[3]], dimnames(self$transform_fit)[[3]])
        )
        out_array <- new_out
      }

      if (inc_generate && !is.null(self$generate_fit)) {
        P1 <- dim(out_array)[3]; P2 <- dim(self$generate_fit)[3]
        I <- dim(out_array)[1]; C <- dim(out_array)[2]
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$generate_fit[, available_chains, , drop = FALSE]
        dimnames(new_out) <- list(
          iteration = dimnames(out_array)[[1]],
          chain = dimnames(out_array)[[2]],
          variable = c(dimnames(out_array)[[3]], dimnames(self$generate_fit)[[3]])
        )
        out_array <- new_out
      }

      P <- dim(out_array)[3]
      param_names <- dimnames(out_array)[[3]]
      if (is.null(param_names)) param_names <- paste0("V", 1:P)
      target_idx <- 1:P

      if (!is.null(pars)) {
        if (is.numeric(pars)) {
          valid_idx <- pars[pars >= 1 & pars <= P]
          if (length(valid_idx) == 0) {
            stop("The index specified in 'pars' was not found.", call. = FALSE)
          }
          target_idx <- valid_idx

        } else if (is.character(pars)) {
          base_names <- gsub("\\[.*\\]$", "", param_names)
          matched <- which(param_names %in% pars | base_names %in% pars)
          if (length(matched) == 0) {
            stop("The variable name specified in 'pars' was not found.", call. = FALSE)
          }
          target_idx <- matched

        } else {
          stop("'pars' must be either numeric or character.", call. = FALSE)
        }
      }

      if (best_only) {
        if (!(self$best_chain %in% available_chains)) {
          stop("The best chain is not included in the selected chains.", call. = FALSE)
        }
        best_pos <- which(available_chains == self$best_chain)[1]
        return(out_array[, best_pos, target_idx, drop = FALSE])
      }

      return(out_array[, , target_idx, drop = FALSE])
    },

    #' @description Summarize posterior draws.
    #' @param pars Character or numeric vector specifying the names or indices of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param max_rows Integer; maximum number of rows to print in the summary table. Default is 10.
    #' @param digits Integer; number of decimal places to print. Default is 2.
    #' @param inc_random Logical; whether to include random effects in the summary. Default is FALSE.
    #' @param inc_transform Logical; whether to include transformed parameters in the summary. Default is TRUE.
    #' @param inc_generate Logical; whether to include generated quantities in the summary. Default is TRUE.
    #' @return A data frame containing the summarized posterior statistics.
    summary = function(pars = NULL, max_rows = 10, digits = 2,
                       inc_random = FALSE, inc_transform = TRUE, inc_generate = TRUE) {

      draws_best <- self$draws(pars = pars,
                               inc_random = inc_random,
                               inc_transform = inc_transform,
                               inc_generate = inc_generate,
                               best_only = TRUE)

      draws_all <- self$draws(pars = pars,
                              inc_random = inc_random,
                              inc_transform = inc_transform,
                              inc_generate = inc_generate,
                              best_only = FALSE)

      P <- dim(draws_best)[3]
      param_names <- dimnames(draws_best)[[3]]
      target_idx <- seq_len(P)

      if (length(target_idx) > 0) {
        current_names <- param_names[target_idx]
        base_names <- gsub("\\[.*\\]$", "", current_names)
        target_views <- c("lp")
        if (!is.null(self$model$view)) {
          target_views <- c(target_views, self$model$view)
        }

        priority_sub_idx <- integer(0)
        for (v in target_views) {
          match_idx <- which(current_names == v | base_names == v)
          priority_sub_idx <- c(priority_sub_idx, match_idx)
        }
        priority_sub_idx <- unique(priority_sub_idx)
        other_sub_idx <- setdiff(seq_along(current_names), priority_sub_idx)
        target_idx <- target_idx[c(priority_sub_idx, other_sub_idx)]
      }

      if (!is.null(max_rows)) {
        limit <- min(length(target_idx), as.integer(max_rows))
        target_idx <- target_idx[seq_len(limit)]
      }

      max_rel_obj <- if (!is.null(self$rel_obj_vals) && length(self$rel_obj_vals) >= self$best_chain) {
        self$rel_obj_vals[self$best_chain]
      } else {
        NA_real_
      }

      res_list_sum <- vector("list", length(target_idx))

      for (i in seq_along(target_idx)) {
        p <- target_idx[i]
        vec_best <- as.vector(draws_best[, 1, p])
        valid_vec <- vec_best[is.finite(vec_best)]
        mat_all <- as.matrix(draws_all[, , p, drop = FALSE])
        var_name <- param_names[p]

        if (length(valid_vec) == 0) {
          res_list_sum[[i]] <- data.frame(
            variable = var_name,
            mean = NA, sd = NA, map = NA,
            q2.5 = NA, q97.5 = NA,
            `Max rel_obj` = max_rel_obj,
            rhat = NA,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          next
        }

        sd_val <- sd(valid_vec)
        if (is.na(sd_val) || sd_val < 1e-10) {
          map_val <- valid_vec[1]
          q95 <- c(valid_vec[1], valid_vec[1])
        } else {
          map_val <- map_est(valid_vec)
          q95 <- quantile95(valid_vec)
        }

        if (ncol(mat_all) >= 2 && all(apply(mat_all, 2, function(z) sum(is.finite(z)) > 1))) {
          rhat_val <- posterior::rhat(mat_all)
        } else {
          rhat_val <- NA_real_
        }

        res_list_sum[[i]] <- data.frame(
          variable = var_name,
          mean = mean(valid_vec),
          sd = sd_val,
          map = map_val,
          q2.5 = unname(q95[1]),
          q97.5 = unname(q95[2]),
          #`Max rel_obj` = max_rel_obj,
          #rhat = rhat_val,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }

      res_df <- do.call(rbind, res_list_sum)
      num_cols <- sapply(res_df, is.numeric)
      res_df[num_cols] <- lapply(res_df[num_cols], function(x) {
        x[abs(x) < 1e-12 & !is.na(x)] <- 0
        x
      })
      class(res_df) <- c("summary_BayesRTMB", "data.frame")
      attr(res_df, "digits") <- digits
      return(res_df)
    },

    #' @description Plot the ELBO history to diagnose convergence.
    #' @param tail_n Integer; the number of recent iterations to plot. If NULL, plots the entire history. Default is 2000.
    #' @param ests Character string `"best"`, numeric vector of estimate indices (e.g., `c(1, 3)`), or `NULL` to plot all. Default is `NULL`.
    #' @param type Character string; the type of plot. Default is "l" (lines).
    #' @param ... Additional arguments passed to the `plot` function.
    #' @return The object itself, invisibly.
    plot_elbo = function(tail_n = 1000, ests = NULL, type = "l", ...) {
      history_list <- self$elbo_history
      num_estimate <- length(history_list)
      if (num_estimate == 0) {
        cat("No ELBO history available.\n")
        return(invisible(self))
      }

      n_total <- length(history_list[[1]])
      plot_data <- do.call(cbind, history_list)
      colnames(plot_data) <- paste0("est", seq_len(num_estimate))

      start_iter <- 1
      if (!is.null(tail_n) && tail_n > 0 && tail_n < n_total) {
        plot_data <- plot_data[(n_total - tail_n + 1):n_total, , drop = FALSE]
        start_iter <- n_total - tail_n + 1
      }

      x_axis <- seq(start_iter, n_total)
      main_title <- if (!is.null(tail_n) && tail_n < n_total) {
        sprintf("ELBO History (Last %d iterations)", tail_n)
      } else {
        "ELBO History"
      }

      best_c <- self$best_chain
      max_rel_obj_val <- self$rel_obj_vals[best_c]
      sub_title <- if(!is.na(max_rel_obj_val)) sprintf("Best (est%d) Final rel_obj: %.5f", best_c, max_rel_obj_val) else ""

      # --- Chain selection for plotting ---
      if (is.null(ests)) {
        target_idx <- seq_len(num_estimate)
      } else if (identical(ests, "best")) {
        target_idx <- best_c
      } else if (is.numeric(ests)) {
        target_idx <- intersect(ests, seq_len(num_estimate))
        if (length(target_idx) == 0) {
          cat("Specified estimates are out of bounds.\n")
          return(invisible(self))
        }
      } else {
        cat("Invalid 'ests' argument. Use NULL, \"best\", or a numeric vector.\n")
        return(invisible(self))
      }

      plot_data <- plot_data[, target_idx, drop = FALSE]

      # --- Line color and thickness adjustment ---
      # Default is gray (black if only 1), highlight blue only if best chain is included
      cols <- rep("gray", length(target_idx))
      lwds <- rep(1, length(target_idx))

      best_in_plot <- which(target_idx == best_c)
      if (length(best_in_plot) > 0) {
        cols[best_in_plot] <- "dodgerblue"
        lwds[best_in_plot] <- 2
      } else if (length(target_idx) == 1) {
        cols <- "black"
      }

      matplot(x_axis, plot_data, type = type, lty = 1, col = cols, lwd = lwds,
              xlab = "Iteration", ylab = "ELBO",
              main = main_title, sub = sub_title, ...)

      # --- Draw legend ---
      if (length(target_idx) > 1) {
        if (length(best_in_plot) > 0) {
          legend("bottomright", legend = c(paste("Best (est", best_c, ")"), "Others"),
                 col = c("dodgerblue", "gray"), lwd = c(2, 1), lty = 1, bty = "n")
        } else {
          legend("bottomright", legend = colnames(plot_data),
                 col = cols, lwd = lwds, lty = 1, bty = "n")
        }
      } else {
        legend("bottomright", legend = colnames(plot_data),
               col = cols, lwd = lwds, lty = 1, bty = "n")
      }

      invisible(self)
    },

    #' @description Plot the parameter trajectory from the final optimization window.
    #' @param pars Character vector specifying the names of parameters to plot. If NULL, plots all available parameters.
    #' @param type Character string; the type of plot. Default is "l" (lines).
    #' @param ... Additional arguments passed to the `matplot` function.
    #' @return The object itself, invisibly.
    plot_trajectory = function(pars = NULL, type = "l", ...) {
      if (is.null(self$mu_history)) {
        cat("No trajectory data available (mu_history is NULL).\n")
        return(invisible(self))
      }

      plot_data <- self$mu_history

      if (!is.null(pars)) {
        valid_pars <- character(0)
        for (p in pars) {
          # 1. Search by exact match (e.g., "delta[1]" or "mat[1,2]")
          if (p %in% colnames(plot_data)) {
            valid_pars <- c(valid_pars, p)
          } else {
            # 2. Bulk search by base name (e.g., "delta" -> "delta[1]", "mat" -> "mat[1,1]")
            pattern <- paste0("^", p, "\\[.*\\]$")
            matched <- grep(pattern, colnames(plot_data), value = TRUE)
            valid_pars <- c(valid_pars, matched)
          }
        }
        valid_pars <- unique(valid_pars)

        if (length(valid_pars) == 0) {
          cat("None of the specified parameters were found in the trajectory data.\n")
          return(invisible(self))
        }
        plot_data <- plot_data[, valid_pars, drop = FALSE]
      }

      n_steps <- nrow(plot_data)
      x_axis <- seq_len(n_steps)

      # --- Automatic adjustment algorithm for vertical axis (ylim) ---
      y_min <- min(plot_data, na.rm = TRUE)
      y_max <- max(plot_data, na.rm = TRUE)
      y_mean <- mean(plot_data, na.rm = TRUE)
      y_range <- y_max - y_min

      args_list <- list(...)

      # Auto-adjust only if ylim is not explicitly specified by user
      if (!"ylim" %in% names(args_list)) {
        # Minimum vertical axis range to secure (max of absolute 0.1 or 5% of mean)
        min_range <- max(0.1, abs(y_mean) * 0.05)

        if (y_range < min_range) {
          # If fluctuation is too small, apply minimum range to make it look "flat"
          y_min <- y_mean - min_range / 2
          y_max <- y_mean + min_range / 2
        } else {
          # If fluctuation is sufficient, add 5% margin to top and bottom
          y_min <- y_min - y_range * 0.05
          y_max <- y_max + y_range * 0.05
        }
        args_list$ylim <- c(y_min, y_max)
      }

      main_title <- sprintf("Parameter Trajectory (Last %d steps)", n_steps)

      # --- Dynamically build arguments for matplot and call ---
      args_list$x <- x_axis
      args_list$y <- plot_data
      args_list$type <- type
      args_list$lty <- 1
      args_list$xlab <- "Iterations (in final window)"
      if (!"ylab" %in% names(args_list)) args_list$ylab <- "Parameter Value"
      if (!"main" %in% names(args_list)) args_list$main <- main_title

      do.call(matplot, args_list)

      if (ncol(plot_data) <= 10) {
        legend("topright", legend = colnames(plot_data),
               col = seq_len(ncol(plot_data)), lty = 1, cex = 0.8,
               bty = "n", inset = c(0.02, 0.02))
      }

      invisible(self)
    },

    #' @description Compute transformed parameters from posterior draws.
    #' @param tran_fn An optional user-supplied function that takes data and parameter lists to return transformed quantities.
    #' @return The `VB_Fit` object itself, invisibly.
    transformed_draws = function(tran_fn = NULL) {
      all_draws <- self$draws(
        inc_random = TRUE,
        inc_transform = FALSE,
        inc_generate = FALSE,
        best_only = FALSE
      )

      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2]

      if (is.null(tran_fn)) tran_fn <- self$model$transform

      wrapper_tran_fn <- function(dat, param) {
        res <- list()

        if (!is.null(tran_fn)) {
          user_res <- tran_fn(dat, param)
          if (is.null(user_res)) user_res <- list()
          res <- c(res, user_res)
        }
        return(res)
      }

      test_para <- constrained_vector_to_list(all_draws[1, 1, -1], self$model$par_list)
      test_tran <- wrapper_tran_fn(self$model$data, test_para)

      if (length(test_tran) == 0) return(invisible(self))

      tran_names <- character(0)
      self$transform_dims <- list()

      for (name in names(test_tran)) {
        val <- test_tran[[name]]
        len <- length(val)
        dim_val <- dim(val)
        if (is.null(dim_val)) dim_val <- len
        self$transform_dims[[name]] <- dim_val

        names_def <- self$model$par_names[[name]]
        flat_nms <- generate_flat_names(name, dim_val, names_def)
        tran_names <- c(tran_names, flat_nms)
      }

      tran_array <- array(NA, dim = c(iter, chains, length(tran_names)))
      dimnames(tran_array) <- list(
        iteration = NULL,
        chain = paste0("est", seq_len(chains)),
        variable = tran_names
      )

      pb <- txtProgressBar(min = 0, max = iter * chains, style = 3)
      counter <- 0
      for (c in seq_len(chains)) {
        for (i in seq_len(iter)) {
          p_list <- constrained_vector_to_list(all_draws[i, c, -1], self$model$par_list)
          res <- wrapper_tran_fn(self$model$data, p_list)
          tran_array[i, c, ] <- unlist(res, use.names = FALSE)
          counter <- counter + 1
          setTxtProgressBar(pb, counter)
        }
      }
      close(pb)

      self$transform_fit <- tran_array
      return(invisible(self))
    },

    #' @description Compute generated quantities from posterior draws.
    #' @param code An `rtmb_code(\{ ... \})` or `\{ ... \}` block containing the logic
    #' to be calculated using posterior samples.
    #' @return The `VB_Fit` object itself (invisibly).
    #' Results are appended to the `generate_fit` field.
    generated_quantities = function(code) {
      raw_code <- substitute(code)

      if (is.name(raw_code)) {
        evaluated <- tryCatch(eval(raw_code, envir = parent.frame()), error = function(e) NULL)
        if (is.language(evaluated) || is.call(evaluated)) {
          code <- evaluated
          raw_code <- evaluated
        }
      }

      if (is.call(raw_code) && identical(raw_code[[1]], as.name("rtmb_code"))) {
        parsed_code <- eval(raw_code, envir = parent.frame())
        if (!"generate" %in% names(parsed_code)) {
          stop("There is no 'generate' block in rtmb_code().")
        }
        gen_ast <- parsed_code$generate
      } else if (is.call(raw_code) && identical(raw_code[[1]], as.name("{"))) {
        gen_ast <- raw_code
      } else if (is.call(code) && identical(code[[1]], as.name("{"))) {
        gen_ast <- code
      } else if (is.list(code) && "generate" %in% names(code)) {
        gen_ast <- code$generate
      } else {
        stop("'code' must be specified in the format rtmb_code(generate = { ... }) or { ... }.")
      }

      gen_fn <- eval(bquote(transform_code(.(gen_ast))))
      environment(gen_fn) <- parent.env(globalenv())

      # 3. Get all samples
      all_draws <- self$draws(inc_random = TRUE, inc_transform = TRUE, inc_generate = FALSE, best_only = FALSE)
      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2]

      # 4. Get dimension information
      test_p_list <- constrained_vector_to_list(all_draws[1, 1, -1], self$model$par_list)
      if (!is.null(self$model$transform)) {
        tran_res <- self$model$transform(self$model$data, test_p_list)
        if (is.list(tran_res)) test_p_list <- c(test_p_list, tran_res)
      }

      if (!is.null(self$generate_fit)) {
        existing_gq <- dimnames(self$generate_fit)[[3]]
        for (g_name in names(self$generate_dims)) {
          g_dim <- self$generate_dims[[g_name]]
          flat_nms <- generate_flat_names(g_name, g_dim, self$model$par_names[[g_name]])
          # Extract only those existing in generate_fit
          if (all(flat_nms %in% existing_gq)) {
            val <- self$generate_fit[1, 1, flat_nms]
            if (length(g_dim) > 1) dim(val) <- g_dim
            test_p_list[[g_name]] <- val
          }
        }
      }
      test_gq <- gen_fn(self$model$data, test_p_list)

      if (is.null(test_gq) || length(test_gq) == 0) return(invisible(self))

      gq_names <- character(0)
      for (name in names(test_gq)) {
        val <- test_gq[[name]]
        dim_val <- if (is.null(dim(val))) length(val) else dim(val)
        self$generate_dims[[name]] <- dim_val
        gq_names <- c(gq_names, generate_flat_names(name, dim_val, self$model$par_names[[name]]))
      }

      # 5. Array for storing results
      new_gq_array <- array(NA, dim = c(iter, chains, length(gq_names)))
      dimnames(new_gq_array) <- list(iteration = NULL, chain = paste0("est", seq_len(chains)), variable = gq_names)

      # 6. Execute for all estimates and samples
      pb <- txtProgressBar(min = 0, max = iter * chains, style = 3)
      counter <- 0
      for (c in seq_len(chains)) {
        for (i in seq_len(iter)) {
          p_list <- constrained_vector_to_list(all_draws[i, c, -1], self$model$par_list)
          if (!is.null(self$model$transform)) {
            tran_res <- self$model$transform(self$model$data, p_list)
            if (is.list(tran_res)) p_list <- c(p_list, tran_res)
          }
          if (!is.null(self$generate_fit)) {
            existing_gq <- dimnames(self$generate_fit)[[3]]
            for (g_name in names(self$generate_dims)) {
              g_dim <- self$generate_dims[[g_name]]
              flat_nms <- generate_flat_names(g_name, g_dim, self$model$par_names[[g_name]])
              # Extract only those existing in generate_fit
              if (all(flat_nms %in% existing_gq)) {
                val <- self$generate_fit[i, c, flat_nms]
                if (length(g_dim) > 1) dim(val) <- g_dim
                p_list[[g_name]] <- val
              }
            }
          }
          res <- gen_fn(self$model$data, p_list)
          new_gq_array[i, c, ] <- unlist(res, use.names = FALSE)
          counter <- counter + 1
          setTxtProgressBar(pb, counter)
        }
      }
      close(pb)

      # 7. Merge with existing results
      if (is.null(self$generate_fit)) {
        self$generate_fit <- new_gq_array
      } else {
        old_gq <- self$generate_fit
        merged_gq <- array(NA, dim = c(iter, chains, dim(old_gq)[3] + dim(new_gq_array)[3]))
        merged_gq[,,1:dim(old_gq)[3]] <- old_gq
        merged_gq[,,(dim(old_gq)[3]+1):dim(merged_gq)[3]] <- new_gq_array
        dimnames(merged_gq) <- list(NULL, dimnames(old_gq)[[2]], c(dimnames(old_gq)[[3]], gq_names))
        self$generate_fit <- merged_gq
      }

      cat("Generated quantities added to samples.\n")
      invisible(self)
    }
  )
)
