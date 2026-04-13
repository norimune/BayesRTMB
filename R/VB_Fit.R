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

  public = list(
    # --- フィールド ---
    model          = NULL, # RTMB_Model のインスタンスへの参照
    fit            = NULL,
    random_fit     = NULL,
    transform_fit       = NULL, # 変換量を保存
    generate_fit         = NULL, # 生成量を保存
    transform_dims      = NULL, # 変換量の次元情報を保存
    generate_dims        = NULL, # GQ変数の次元情報を保存
    elbo_history   = NULL, # ELBOの推移を保存
    laplace        = NULL,
    posterior_mean = NULL,
    ELBO           = NULL,
    rel_obj_vals   = NULL,
    best_chain     = NULL,
    mu_history     = NULL,

    # 1. コンストラクタ
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
    #' @param inc_random Logical; whether to include random effects in the output. Default is FALSE.
    #' @param inc_transform Logical; whether to include transformed parameters in the output. Default is TRUE.
    #' @param inc_generate Logical; whether to include generated quantities in the output. Default is TRUE.
    #' @param best_only Logical; whether to extract only from the chain with the maximum ELBO. Default is TRUE.
    #' @return A 3D array of posterior draws `[iterations, chains, parameters]`.
    draws = function(pars = NULL, inc_random = FALSE, inc_transform = TRUE, inc_generate = TRUE, best_only = TRUE) {
      out_array <- self$fit

      if (inc_random && !is.null(self$random_fit)) {
        P1 <- dim(out_array)[3]; P2 <- dim(self$random_fit)[3]
        I <- dim(out_array)[1]; C <- dim(out_array)[2]
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$random_fit
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
        new_out[,,(P1+1):(P1+P2)] <- self$transform_fit
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
        new_out[,,(P1+1):(P1+P2)] <- self$generate_fit
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
            stop("`pars` に指定されたインデックスが見つかりません。", call. = FALSE)
          }
          target_idx <- valid_idx

        } else if (is.character(pars)) {
          base_names <- gsub("\\[.*\\]$", "", param_names)
          matched <- which(param_names %in% pars | base_names %in% pars)
          if (length(matched) == 0) {
            stop("`pars` に指定された変数名が見つかりません。", call. = FALSE)
          }
          target_idx <- matched

        } else {
          stop("`pars` は numeric か character で指定してください。", call. = FALSE)
        }
      }

      if (best_only) {
        return(out_array[, self$best_chain, target_idx, drop = FALSE])
      } else {
        return(out_array[, , target_idx, drop = FALSE])
      }
    },

    #' @description Summarize posterior draws. (Note: Rhat and ESS are not computed for ADVI).
    #' @param pars Character or numeric vector specifying the names or indices of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param max_rows Integer; maximum number of rows to print in the summary table. Default is 10.
    #' @param digits Integer; number of decimal places to print. Default is 2.
    #' @param inc_random Logical; whether to include random effects in the summary. Default is FALSE.
    #' @param inc_transform Logical; whether to include transformed parameters in the summary. Default is TRUE.
    #' @param inc_generate Logical; whether to include generated quantities in the summary. Default is TRUE.
    #' @return A data frame containing the summarized posterior statistics.
    summary = function(pars = NULL, max_rows = 10, digits = 2,
                       inc_random = FALSE, inc_transform = TRUE, inc_generate = TRUE) {

      draws_array <- self$draws(pars = pars,
                                inc_random = inc_random,
                                inc_transform = inc_transform,
                                inc_generate = inc_generate,
                                best_only = TRUE)

      P <- dim(draws_array)[3]
      param_names <- dimnames(draws_array)[[3]]

      target_idx <- 1:P

      # --- lp と model$view による優先並び替え ---
      if (length(target_idx) > 0) {
        current_names <- param_names[target_idx]
        base_names <- gsub("\\[.*\\]$", "", current_names)

        # lp を常に最優先にする
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
        target_idx <- target_idx[1:limit]
      }

      res_list_sum <- vector("list", length(target_idx))

      for (i in seq_along(target_idx)) {
        p <- target_idx[i]

        mat_best <- as.matrix(draws_array[, 1, p])
        vec_best <- as.vector(mat_best)
        valid_vec <- vec_best[is.finite(vec_best)]

        var_name <- param_names[p]

        # 軌跡データから変動幅を計算
        if (!is.null(self$mu_history) && var_name %in% colnames(self$mu_history)) {
          hist_vec <- self$mu_history[, var_name]
          n_hist <- length(hist_vec)
          half_n <- floor(n_hist / 2)

          mean_all   <- mean(hist_vec)
          mean_first <- mean(hist_vec[1:half_n])
          mean_last  <- mean(hist_vec[(half_n + 1):n_hist])

          scale_factor <- abs(mean_all) + 1e-4
          drift_val <- abs(mean_last - mean_first) / scale_factor
          rel_sd_val <- sd(hist_vec) / scale_factor

        } else {
          drift_val <- NA_real_
          rel_sd_val <- NA_real_
        }

        if (length(valid_vec) == 0) {
          res_list_sum[[i]] <- data.frame(
            variable = var_name, mean = NA, sd = NA, map = NA,
            q2.5 = NA, q97.5 = NA, opt_diff = NA,
            stringsAsFactors = FALSE, check.names = FALSE
          )
          next
        }

        sd_val <- sd(valid_vec)
        if (is.na(sd_val) || sd_val < 1e-10) {
          map_val <- valid_vec[1]; q95 <- c(valid_vec[1], valid_vec[1])
          diff_sd_ratio <- NA
        } else {
          map_val   <- map_est(valid_vec)
          q95       <- quantile95(valid_vec)
        }

        res_list_sum[[i]] <- data.frame(
          variable = param_names[p],
          mean     = mean(valid_vec),
          sd       = sd_val,
          map      = map_val,
          q2.5     = unname(q95[1]),
          q97.5    = unname(q95[2]),
          #drift    = drift_val,   # NAの場合はそのままNAになるためif文は不要です
          #rel_sd   = rel_sd_val,
          stringsAsFactors = FALSE
        )
      }
      res_df <- do.call(rbind, res_list_sum)
      class(res_df) <- c("summary_BayesRTMB", "data.frame")
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

      # --- 描画するチェインの選択処理 ---
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

      # --- 線の色と太さの調整 ---
      # 基本はグレー（1本のみの場合は黒）で、ベストチェインが含まれる場合のみ青くハイライトする
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

      # --- 凡例の描画 ---
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
          # 1. 完全一致による検索 (例: "delta[1]" や "mat[1,2]")
          if (p %in% colnames(plot_data)) {
            valid_pars <- c(valid_pars, p)
          } else {
            # 2. ベース名による一括検索 (例: "delta" -> "delta[1]", "mat" -> "mat[1,1]")
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

      # --- 縦軸(ylim)の自動調整アルゴリズム ---
      y_min <- min(plot_data, na.rm = TRUE)
      y_max <- max(plot_data, na.rm = TRUE)
      y_mean <- mean(plot_data, na.rm = TRUE)
      y_range <- y_max - y_min

      args_list <- list(...)

      # ユーザーが明示的に ylim を指定していない場合のみ自動調整
      if (!"ylim" %in% names(args_list)) {
        # 最小限確保する縦軸の幅 (絶対値0.1 または 平均値の5% のどちらか大きい方)
        min_range <- max(0.1, abs(y_mean) * 0.05)

        if (y_range < min_range) {
          # 変動が小さすぎる場合は、最小幅を適用して「平坦」に見せる
          y_min <- y_mean - min_range / 2
          y_max <- y_mean + min_range / 2
        } else {
          # 変動が十分ある場合は、上下に5%の余白を持たせる
          y_min <- y_min - y_range * 0.05
          y_max <- y_max + y_range * 0.05
        }
        args_list$ylim <- c(y_min, y_max)
      }

      main_title <- sprintf("Parameter Trajectory (Last %d steps)", n_steps)

      # --- matplotの引数を動的に構築して呼び出し ---
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
    #' @param code An `rtmb_code({ ... })` or `{ ... }` block containing the logic
    #' to be calculated using posterior samples.
    #' @return The `VB_Fit` object itself (invisibly).
    #' Results are appended to the `generate_fit` field.
    generated_quantities = function(code) {
      # 1. 引数のキャプチャとASTの抽出 (MCMC_Fitと同様)
      raw_code <- substitute(code)
      if (is.name(raw_code)) {
        evaluated <- tryCatch(eval(raw_code, envir = parent.frame()), error = function(e) NULL)
        if (is.language(evaluated) || is.call(evaluated)) code <- evaluated
      }

      if (is.call(code) && identical(code[[1]], as.name("rtmb_code"))) {
        gen_ast <- code$generate
      } else {
        gen_ast <- code
      }

      # 2. 関数化
      gen_fn <- eval(bquote(transform_code(.(gen_ast))))
      environment(gen_fn) <- parent.env(globalenv())

      # 3. 全サンプルの取得
      all_draws <- self$draws(inc_random = TRUE, inc_transform = TRUE, inc_generate = FALSE, best_only = FALSE)
      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2]

      # 4. 次元情報の取得
      test_p_list <- constrained_vector_to_list(all_draws[1, 1, -1], self$model$par_list)
      if (!is.null(self$model$transform)) {
        tran_res <- self$model$transform(self$model$data, test_p_list)
        if (is.list(tran_res)) test_p_list <- c(test_p_list, tran_res)
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

      # 5. 結果格納用の配列
      new_gq_array <- array(NA, dim = c(iter, chains, length(gq_names)))
      dimnames(new_gq_array) <- list(iteration = NULL, chain = paste0("est", seq_len(chains)), variable = gq_names)

      # 6. 全エスティメイト・全サンプルに対して実行
      pb <- txtProgressBar(min = 0, max = iter * chains, style = 3)
      counter <- 0
      for (c in seq_len(chains)) {
        for (i in seq_len(iter)) {
          p_list <- constrained_vector_to_list(all_draws[i, c, -1], self$model$par_list)
          if (!is.null(self$model$transform)) {
            tran_res <- self$model$transform(self$model$data, p_list)
            if (is.list(tran_res)) p_list <- c(p_list, tran_res)
          }
          res <- gen_fn(self$model$data, p_list)
          new_gq_array[i, c, ] <- unlist(res, use.names = FALSE)
          counter <- counter + 1
          setTxtProgressBar(pb, counter)
        }
      }
      close(pb)

      # 7. 既存の結果とマージ
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
    },

    #' @description Apply internal rotation to sampled parameters.
    #' @param target Character string specifying the target variable to base the rotation on.
    #' @param method Character string specifying the rotation method.
    #' @param type Character string specifying the rotation type.
    #' @param linked_straight Character vector of variable names to be rotated in the same direction.
    #' @param linked_inverse Character vector of variable names to be rotated in the inverse direction.
    #' @param overwrite Logical; whether to overwrite the stored draws. If FALSE, adds to generated quantities. Default is FALSE.
    #' @param suffix Character string to append to the rotated variable names when overwrite is FALSE.
    #' @param ... Additional arguments passed to the rotation function.
    #' @return The updated object invisibly.
    internal_rotate = function(target, method = "procrustes", type = "orthogonal",
                               linked_straight = NULL, linked_inverse = NULL,
                               overwrite = FALSE, suffix = "rot", ...) {

      f_arr <- self$fit
      r_arr <- self$random_fit
      t_arr <- self$transform_fit
      g_arr <- self$generate_fit

      v_names_f <- dimnames(f_arr)[[3]]
      v_names_r <- if (!is.null(r_arr)) dimnames(r_arr)[[3]] else character(0)
      v_names_t <- if (!is.null(t_arr)) dimnames(t_arr)[[3]] else character(0)
      v_names_g <- if (!is.null(g_arr)) dimnames(g_arr)[[3]] else character(0)

      get_var_info <- function(vname) {
        pattern <- paste0("^", vname, "\\[")
        idx_f <- grep(pattern, v_names_f)
        if (length(idx_f) > 0) return(list(loc = "fixed", idx = idx_f, dim = self$model$par_list[[vname]]$dim))
        idx_r <- grep(pattern, v_names_r)
        if (length(idx_r) > 0) return(list(loc = "random", idx = idx_r, dim = self$model$par_list[[vname]]$dim))
        idx_t <- grep(pattern, v_names_t)
        if (length(idx_t) > 0) return(list(loc = "transform", idx = idx_t, dim = self$transform_dims[[vname]]))
        idx_g <- grep(pattern, v_names_g)
        if (length(idx_g) > 0) return(list(loc = "generate", idx = idx_g, dim = self$generate_dims[[vname]]))

        if (vname %in% v_names_f) return(list(loc = "fixed", idx = which(v_names_f == vname), dim = 1))
        if (vname %in% v_names_r) return(list(loc = "random", idx = which(v_names_r == vname), dim = 1))
        if (vname %in% v_names_t) return(list(loc = "transform", idx = which(v_names_t == vname), dim = 1))
        if (vname %in% v_names_g) return(list(loc = "generate", idx = which(v_names_g == vname), dim = 1))

        stop(paste0("Rotation failed: Variable '", vname, "' not found."))
      }

      t_info <- get_var_info(target)
      if (length(t_info$dim) != 2) stop(paste0("Target variable '", target, "' must be a matrix."))
      R_t <- t_info$dim[1]; C_t <- t_info$dim[2]

      # MAP推定値 (最大lp) を取得して基準とする
      lp_mat <- f_arr[, , 1]
      if (is.null(dim(lp_mat))) { # VB_Fit用のフォールバック
        best_iter <- which.max(lp_mat)
        best_chain <- 1
      } else {
        max_idx <- which(lp_mat == max(lp_mat, na.rm = TRUE), arr.ind = TRUE)
        best_iter <- max_idx[1, 1]
        best_chain <- max_idx[1, 2]
      }

      if (t_info$loc == "fixed") Y_vec <- f_arr[best_iter, best_chain, t_info$idx]
      else if (t_info$loc == "random") Y_vec <- r_arr[best_iter, best_chain, t_info$idx]
      else if (t_info$loc == "transform") Y_vec <- t_arr[best_iter, best_chain, t_info$idx]
      else Y_vec <- g_arr[best_iter, best_chain, t_info$idx]

      X_map <- matrix(Y_vec, nrow = R_t, ncol = C_t)

      # --- 変換行列の決定 ---
      if (method == "procrustes") {
        Th_straight <- diag(1, C_t)
        Th_inverse <- diag(1, C_t)
      } else {
        if (!requireNamespace("GPArotation", quietly = TRUE)) stop("GPArotation is required.")
        rot_fn <- tryCatch(
          match.fun(method),
          error = function(e) {
            if (exists(method, where = asNamespace("GPArotation"), mode = "function")) return(getFromNamespace(method, "GPArotation"))
            stop("Rotation method not found.")
          }
        )
        map_rot <- rot_fn(X_map, ...)
        if (type == "orthogonal") {
          Th_straight <- map_rot$Th
          Th_inverse <- map_rot$Th
        } else if (type == "oblique") {
          Th_straight <- solve(t(map_rot$Th))
          Th_inverse <- map_rot$Th
        }
      }

      iter_total <- dim(f_arr)[1]
      chains <- dim(f_arr)[2]

      if (!isTRUE(overwrite)) {
        all_vars <- c(target, linked_straight, linked_inverse)
        new_names <- character(0)
        new_dims <- list()
        new_idx_map <- list()
        current_p <- 1

        for (v in all_vars) {
          info <- get_var_info(v)
          # 変更点: "_rot" 固定ではなく "_suffix" を使う
          v_rot <- paste0(v, "_", suffix)
          new_dims[[v_rot]] <- info$dim

          if (info$loc == "fixed") orig_names <- dimnames(f_arr)[[3]][info$idx]
          else if (info$loc == "random") orig_names <- dimnames(r_arr)[[3]][info$idx]
          else if (info$loc == "transform") orig_names <- dimnames(t_arr)[[3]][info$idx]
          else orig_names <- dimnames(g_arr)[[3]][info$idx]

          new_flat_nms <- gsub(paste0("^", v, "(?=\\[|$)"), v_rot, orig_names, perl = TRUE)
          new_names <- c(new_names, new_flat_nms)
          new_idx_map[[v]] <- current_p:(current_p + length(new_flat_nms) - 1)
          current_p <- current_p + length(new_flat_nms)
        }
        rot_arr <- array(NA, dim = c(iter_total, chains, length(new_names)))
        dimnames(rot_arr) <- list(iteration = NULL, chain = dimnames(f_arr)[[2]], variable = new_names)
      }

      for (c in 1:chains) {
        for (i in 1:iter_total) {

          N_draws <- iter_total * chains

          # ターゲット変数の抽出 (N_draws x (R_t * C_t) の行列にする)
          if (t_info$loc == "fixed") X_all <- matrix(f_arr[,,t_info$idx], nrow = N_draws)
          else if (t_info$loc == "random") X_all <- matrix(r_arr[,,t_info$idx], nrow = N_draws)
          else if (t_info$loc == "transform") X_all <- matrix(t_arr[,,t_info$idx], nrow = N_draws)
          else X_all <- matrix(g_arr[,,t_info$idx], nrow = N_draws)

          # 行ごとに回転処理を適用
          rotated_list <- lapply(1:N_draws, function(idx) {
            X <- matrix(X_all[idx, ], nrow = R_t, ncol = C_t)
            svd_out <- svd(t(X) %*% X_map)
            R_proc <- svd_out$u %*% t(svd_out$v)
            X_rot <- as.numeric((X %*% R_proc) %*% Th_straight)

            # 戻り値として、回転済みターゲット変数と R_proc を返す
            list(X_rot = X_rot, R_proc = R_proc)
          })

          # 結果を元の配列に戻す
          for (idx in 1:N_draws) {
            # 1次元インデックスから i(iter) と c(chain) を逆算
            i <- (idx - 1) %% iter_total + 1
            c <- (idx - 1) %/% iter_total + 1

            res_rot <- rotated_list[[idx]]
            R_proc <- res_rot$R_proc

            if (isTRUE(overwrite)) {
              if (t_info$loc == "fixed") f_arr[i, c, t_info$idx] <- res_rot$X_rot
              else if (t_info$loc == "random") r_arr[i, c, t_info$idx] <- res_rot$X_rot
              else if (t_info$loc == "transform") t_arr[i, c, t_info$idx] <- res_rot$X_rot
              else g_arr[i, c, t_info$idx] <- res_rot$X_rot
            } else {
              rot_arr[i, c, new_idx_map[[target]]] <- res_rot$X_rot
            }

            # リンクされた変数の処理 (R_proc を再利用)
            process_linked <- function(linked_vars, Th_mat) {
              if (is.null(linked_vars)) return()
              for (lvar in linked_vars) {
                l_info <- get_var_info(lvar)
                if (length(l_info$dim) != 2 || l_info$dim[2] != C_t) next

                if (l_info$loc == "fixed") Z_vec <- f_arr[i, c, l_info$idx]
                else if (l_info$loc == "random") Z_vec <- r_arr[i, c, l_info$idx]
                else if (l_info$loc == "transform") Z_vec <- t_arr[i, c, l_info$idx]
                else Z_vec <- g_arr[i, c, l_info$idx]

                Z <- matrix(Z_vec, nrow = l_info$dim[1], ncol = C_t)
                Z_rot <- as.numeric((Z %*% R_proc) %*% Th_mat)

                if (isTRUE(overwrite)) {
                  if (l_info$loc == "fixed") f_arr[i, c, l_info$idx] <- Z_rot
                  else if (l_info$loc == "random") r_arr[i, c, l_info$idx] <- Z_rot
                  else if (l_info$loc == "transform") t_arr[i, c, l_info$idx] <- Z_rot
                  else g_arr[i, c, l_info$idx] <- Z_rot
                } else {
                  rot_arr[i, c, new_idx_map[[lvar]]] <- Z_rot
                }
              }
            }

            process_linked(linked_straight, Th_straight)
            process_linked(linked_inverse, Th_inverse)
          }
        }
      }

      if (isTRUE(overwrite)) {
        self$fit <- f_arr
        self$random_fit <- r_arr
        self$transform_fit <- t_arr
        self$generate_fit <- g_arr

        fixed_mean_new <- apply(self$fit[, , -1, drop = FALSE], 3, mean)
        new_posterior_mean <- self$posterior_mean
        new_posterior_mean[names(fixed_mean_new)] <- fixed_mean_new
        if (!is.null(self$random_fit)) {
          random_mean_new <- apply(self$random_fit, 3, mean)
          new_posterior_mean[names(random_mean_new)] <- random_mean_new
        }
        self$posterior_mean <- new_posterior_mean
      } else {
        if (is.null(self$generate_fit)) {
          self$generate_fit <- rot_arr
        } else {
          old_gq <- self$generate_fit
          I <- dim(old_gq)[1]; C <- dim(old_gq)[2]
          P1 <- dim(old_gq)[3]; P2 <- dim(rot_arr)[3]
          new_gq <- array(NA, dim = c(I, C, P1 + P2))
          new_gq[,,1:P1] <- old_gq
          new_gq[,,(P1+1):(P1+P2)] <- rot_arr
          dimnames(new_gq) <- list(dimnames(old_gq)[[1]], dimnames(old_gq)[[2]], c(dimnames(old_gq)[[3]], dimnames(rot_arr)[[3]]))
          self$generate_fit <- new_gq
        }
        for (v in names(new_dims)) {
          self$generate_dims[[v]] <- new_dims[[v]]
        }
      }

      return(invisible(self))
    },

    #' @description Performs orthogonal Procrustes rotation on posterior samples
    #' based on a specified reference (the MAP value at the point of maximum ELBO).
    #' @param target Character string specifying the name of the matrix parameter
    #' to be used as the rotation reference.
    #' @param linked A character vector of other parameter names to be rotated
    #' in the same direction as the target.
    #' @return The `VB_Fit` object itself (invisibly).
    #' Rotated values are saved to `generate_fit` with a `_rot` suffix added
    #' to the variable names.
    rotate = function(target, linked = NULL) {
      cat("Applying orthogonal Procrustes rotation (Saving to generate as _rot)...\n")

      t_info <- self$model$par_list[[target]]
      t_dim <- if (is.null(t_info)) self$transform_dims[[target]] else t_info$dim

      # 1. ターゲットとLPの事後サンプルを取得 (全エスティメイト)
      target_draws <- self$draws(pars = target, inc_transform = TRUE, inc_generate = TRUE, best_only = FALSE)
      lp_draws <- self$draws(pars = "lp", best_only = FALSE)

      # 2. 全エスティメイトの中で最大ELBO/LPを持つ箇所を特定 (基準値の抽出)
      max_idx <- which(lp_draws == max(lp_draws, na.rm = TRUE), arr.ind = TRUE)
      target_map <- target_draws[max_idx[1,1], max_idx[1,2], ]
      dim(target_map) <- t_dim

      # 3. 基準値をモデルのデータに一時保存
      ref_name <- paste0(target, "_ref")
      self$model$data[[ref_name]] <- target_map

      # 4. ASTの構築
      exprs <- list()
      exprs[[length(exprs) + 1]] <- bquote(svd_res <- svd(t(.(as.name(target))) %*% .(as.name(ref_name))))
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

      # 5. 実行
      self$generated_quantities(as.call(c(list(as.name("{")), exprs)))
      invisible(self)
    },

    #' @description Rotates factor loadings and optionally rotates associated parameters.
    #' @param target Character string specifying the factor loadings matrix to
    #' base the rotation on. Defaults to "loadings".
    #' @param linked Character vector of parameters (e.g., item parameters)
    #' to which the same rotation matrix should be applied.
    #' @param scores Character vector of parameters (e.g., factor scores)
    #' to which the inverse-transpose of the rotation matrix should be applied.
    #' @param rotate Character string specifying the rotation method (e.g., "promax",
    #' "varimax", "oblimin"). Supports `GPArotation` functions or "promax".
    #' @param ... Additional arguments passed to the rotation function.
    #' @return The `VB_Fit` object itself (invisibly).
    #' Rotated values are saved to `generate_fit` with a suffix indicating
    #' the rotation method.
    fa_rotate = function(target, linked = NULL, scores = NULL, rotate = "promax") {
      cat(sprintf("Applying %s rotation to %s (Saving to generate as _%s)...\n", rotate, target, rotate))

      # 回転関数の特定
      if (exists(rotate, mode = "function")) {
        fn_call <- as.name(rotate)
      } else {
        fn_call <- call("::", as.name("GPArotation"), as.name(rotate))
      }

      # 回転タイプとPhiの有無をテスト
      all_draws <- self$draws(pars = target, inc_transform = TRUE, best_only = TRUE)
      dummy_L <- all_draws[1, 1, ]
      dim(dummy_L) <- if (is.null(self$model$par_list[[target]])) self$transform_dims[[target]] else self$model$par_list[[target]]$dim

      test_rot <- eval(as.call(list(fn_call, dummy_L)))
      is_matrix_rot <- is.matrix(test_rot)

      # ASTの構築
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

      # 連動変数の処理 (linked=同方向, scores=逆転置)
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

      # 実行
      self$generated_quantities(as.call(c(list(as.name("{")), exprs)))
      invisible(self)
    }
  )
)
