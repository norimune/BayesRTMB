#' VB fit object
#'
#' An R6 class storing posterior samples and related information
#' from Automatic Differentiation Variational Inference (ADVI).
#'
#' @field model An `RTMB_Model` object used for estimation.
#' @field fit A 3D array of posterior draws for fixed model parameters.
#' @field random_fit A 3D array of posterior draws for random effects, if available.
#' @field tran_fit A 3D array of posterior draws for transformed parameters, if available.
#' @field gq_fit A 3D array of posterior draws for generated quantities, if available.
#' @field tran_dims A list storing dimension information for transformed parameters.
#' @field gq_dims A list storing dimension information for generated quantities.
#' @field elbo_history A numeric vector storing the Evidence Lower Bound (ELBO) history during optimization.
#' @field laplace Logical; whether Laplace approximation was used to marginalize random effects.
#' @field posterior_mean A named numeric vector of posterior mean estimates.
#'
#' @export
VB_Fit <- R6::R6Class(
  classname = "advi_fit",

  public = list(
    # --- フィールド ---
    model          = NULL, # RTMB_Model のインスタンスへの参照
    fit            = NULL,
    random_fit     = NULL,
    tran_fit       = NULL, # 変換量を保存
    gq_fit         = NULL, # 生成量を保存
    tran_dims      = NULL, # 変換量の次元情報を保存
    gq_dims        = NULL, # GQ変数の次元情報を保存
    elbo_history   = NULL, # ELBOの推移を保存
    laplace        = NULL,
    posterior_mean = NULL,

    # 1. コンストラクタ
    #' @description Create a new `VB_Fit` object.
    #' @param model An `RTMB_Model` object.
    #' @param fit A 3D array of parameter draws.
    #' @param random_fit A 3D array of random effect draws.
    #' @param elbo_history A numeric vector of ELBO values.
    #' @param laplace Logical; indicates if Laplace approximation was used.
    #' @param posterior_mean A named numeric vector of posterior means.
    initialize = function(model, fit, random_fit, elbo_history, laplace, posterior_mean) {
      self$model <- model
      self$fit <- fit
      self$random_fit <- random_fit
      self$elbo_history <- elbo_history
      self$laplace <- laplace
      self$posterior_mean <- posterior_mean
      self$tran_fit <- NULL
      self$tran_dims <- list()
      self$gq_fit <- NULL
      self$gq_dims <- list()
    },

    #' @description Print a brief summary of the fitted object.
    #' @param ... Additional arguments passed to the `summary` method.
    #' @return The object itself, invisibly.
    print = function(...) {
      out <- self$summary(...)
      base::print(out)
      invisible(self)
    },

    #' @description Extract posterior draws for selected parameters.
    #' @param pars Character or numeric vector specifying the names or indices of parameters to extract. If NULL, all available parameters are extracted.
    #' @param inc_random Logical; whether to include random effects in the output. Default is FALSE.
    #' @param inc_tran Logical; whether to include transformed parameters in the output. Default is TRUE.
    #' @param inc_gq Logical; whether to include generated quantities in the output. Default is TRUE.
    #' @return A 3D array of posterior draws `[iterations, chains, parameters]`.
    draws = function(pars = NULL, inc_random = FALSE, inc_tran = TRUE, inc_gq = TRUE) {
      out_array <- self$fit

      if (inc_random && !is.null(self$random_fit)) {
        P1 <- dim(out_array)[3]
        P2 <- dim(self$random_fit)[3]
        I <- dim(out_array)[1]
        C <- 1 # ADVIはチェイン1固定
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

      if (inc_tran && !is.null(self$tran_fit)) {
        P1 <- dim(out_array)[3]
        P2 <- dim(self$tran_fit)[3]
        I <- dim(out_array)[1]
        C <- 1
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$tran_fit
        dimnames(new_out) <- list(
          iteration = dimnames(out_array)[[1]],
          chain = dimnames(out_array)[[2]],
          variable = c(dimnames(out_array)[[3]], dimnames(self$tran_fit)[[3]])
        )
        out_array <- new_out
      }

      if (inc_gq && !is.null(self$gq_fit)) {
        P1 <- dim(out_array)[3]
        P2 <- dim(self$gq_fit)[3]
        I <- dim(out_array)[1]
        C <- 1
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$gq_fit
        dimnames(new_out) <- list(
          iteration = dimnames(out_array)[[1]],
          chain = dimnames(out_array)[[2]],
          variable = c(dimnames(out_array)[[3]], dimnames(self$gq_fit)[[3]])
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

      return(out_array[, 1, target_idx, drop = FALSE])
    },

    #' @description Summarize posterior draws. (Note: Rhat and ESS are not computed for ADVI).
    #' @param pars Character or numeric vector specifying the names or indices of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param max_rows Integer; maximum number of rows to print in the summary table. Default is 10.
    #' @param digits Integer; number of decimal places to print. Default is 2.
    #' @param inc_random Logical; whether to include random effects in the summary. Default is FALSE.
    #' @param inc_tran Logical; whether to include transformed parameters in the summary. Default is TRUE.
    #' @param inc_gq Logical; whether to include generated quantities in the summary. Default is TRUE.
    #' @return A data frame containing the summarized posterior statistics.
    summary = function(pars = NULL, max_rows = 10, digits = 2,
                       inc_random = FALSE, inc_tran = TRUE, inc_gq = TRUE){
      draws_array <- self$draws(pars = pars,
                                inc_random = inc_random,
                                inc_tran = inc_tran,
                                inc_gq = inc_gq)

      P <- dim(draws_array)[3]
      param_names <- dimnames(draws_array)[[3]]

      target_idx <- 1:P
      if (!is.null(max_rows)) {
        limit <- min(length(target_idx), as.integer(max_rows))
        target_idx <- target_idx[1:limit]
      }

      res_list_sum <- vector("list", length(target_idx))

      for (i in seq_along(target_idx)) {
        p <- target_idx[i]
        mat_p <- as.matrix(draws_array[, 1, p]) # チェインは1固定
        vec_p <- as.vector(mat_p)
        valid_vec <- vec_p[is.finite(vec_p)]

        if (length(valid_vec) == 0) {
          res_list_sum[[i]] <-
            data.frame(
              variable = param_names[p],
              mean = NA,
              sd = NA,
              map = NA,
              q2.5 = NA,
              q97.5 = NA,
              stringsAsFactors = FALSE
            )
          next
        }

        sd_val <- sd(valid_vec)
        if (is.na(sd_val) || sd_val < 1e-10) {
          map_val <- valid_vec[1]; q95 <- c(valid_vec[1], valid_vec[1])
        } else {
          # MCMC_Fit等で定義されているmap_est, quantile95関数を利用する前提
          map_val   <- map_est(valid_vec)
          q95       <- quantile95(valid_vec)
        }

        res_list_sum[[i]] <- data.frame(
          variable = param_names[p],
          mean     = round(mean(valid_vec), digits),
          sd       = round(sd_val, digits),
          map      = round(map_val, digits),
          q2.5     = round(unname(q95[1]), digits),
          q97.5    = round(unname(q95[2]), digits),
          row.names = NULL,
          stringsAsFactors = FALSE
        )
      }
      return(do.call(rbind, res_list_sum))
    },

    #' @description Compute transformed parameters from posterior draws.
    #' @param tran_fn An optional user-supplied function that takes data and parameter lists to return transformed quantities.
    #' @return The `VB_Fit` object itself, invisibly.
    transformed_draws = function(tran_fn = NULL) {
      all_draws <- self$draws(
        inc_random = TRUE,
        inc_tran = FALSE,
        inc_gq = FALSE
      )

      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2] # 基本的に1

      wrapper_tran_fn <- function(dat, param) {
        res <- list()
        for (name in names(self$model$par_list)) {
          if (self$model$par_list[[name]]$type == "CF_corr") {
            mat_name <- if (grepl("^CF_", name)) sub("^CF_", "", name) else paste0(name, "_corr")
            res[[mat_name]] <- param[[name]] %*% t(param[[name]])
          }
        }
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
      self$tran_dims <- list()

      for (name in names(test_tran)) {
        val <- test_tran[[name]]
        len <- length(val)
        dim_val <- dim(val)
        if (is.null(dim_val)) dim_val <- len
        self$tran_dims[[name]] <- dim_val

        if (len == 1) {
          tran_names <- c(tran_names, name)
        } else {
          if (length(dim_val) == 1) {
            tran_names <- c(tran_names, paste0(name, "[", seq_len(len), "]"))
          } else {
            grid <- expand.grid(lapply(dim_val, seq_len))
            indices <- apply(grid, 1, paste, collapse = ",")
            tran_names <- c(tran_names, paste0(name, "[", indices, "]"))
          }
        }
      }

      tran_array <- array(NA, dim = c(iter, chains, length(tran_names)))
      dimnames(tran_array) <- list(
        iteration = NULL,
        chain = paste0("chain", seq_len(chains)),
        variable = tran_names
      )

      cat("Calculating transformed parameters...\n")
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

      self$tran_fit <- tran_array
      return(invisible(self))
    },

    #' @description Compute generated quantities from posterior draws.
    #' @param gq_fn An optional user-supplied function that takes data and parameter lists to return generated quantities.
    #' @return The `VB_Fit` object itself, invisibly.
    generated_quantities = function(gq_fn = NULL) {
      all_draws <- self$draws(
        inc_random = TRUE,
        inc_tran = FALSE,
        inc_gq = FALSE
      )
      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2]

      if (is.null(gq_fn)) return(invisible(self))

      test_para <- constrained_vector_to_list(all_draws[1, 1, -1], self$model$par_list)
      test_gq <- gq_fn(self$model$data, test_para)

      if (is.null(test_gq) || length(test_gq) == 0) return(invisible(self))

      gq_names <- character(0)
      self$gq_dims <- list()

      for (name in names(test_gq)) {
        val <- test_gq[[name]]
        len <- length(val)
        dim_val <- dim(val)
        if (is.null(dim_val)) dim_val <- len
        self$gq_dims[[name]] <- dim_val

        if (len == 1) {
          gq_names <- c(gq_names, name)
        } else {
          if (length(dim_val) == 1) {
            gq_names <- c(gq_names, paste0(name, "[", seq_len(len), "]"))
          } else {
            grid <- expand.grid(lapply(dim_val, seq_len))
            indices <- apply(grid, 1, paste, collapse = ",")
            gq_names <- c(gq_names, paste0(name, "[", indices, "]"))
          }
        }
      }

      gq_array <- array(NA, dim = c(iter, chains, length(gq_names)))
      dimnames(gq_array) <- list(
        iteration = NULL,
        chain = paste0("chain", seq_len(chains)),
        variable = gq_names
      )

      cat("Calculating generated quantities...\n")
      pb <- txtProgressBar(min = 0, max = iter * chains, style = 3)
      counter <- 0

      for (c in seq_len(chains)) {
        for (i in seq_len(iter)) {
          p_list <- constrained_vector_to_list(all_draws[i, c, -1], self$model$par_list)
          res <- gq_fn(self$model$data, p_list)
          gq_array[i, c, ] <- unlist(res, use.names = FALSE)

          counter <- counter + 1
          setTxtProgressBar(pb, counter)
        }
      }
      close(pb)

      self$gq_fit <- gq_array
      return(invisible(self))
    },
    #' @description Apply internal rotation to sampled parameters.
    #' @param target Character string specifying the target variable to base the rotation on.
    #' @param method Character string specifying the rotation method (e.g., "procrustes", "promax"). Default is "procrustes".
    #' @param type Character string specifying the rotation type ("orthogonal" or "oblique"). Default is "orthogonal".
    #' @param linked_straight Character vector of variable names to be rotated in the same direction. Default is NULL.
    #' @param linked_inverse Character vector of variable names to be rotated in the inverse direction. Default is NULL.
    #' @param overwrite Logical; whether to overwrite the stored draws in the current object. Default is NULL.
    #' @param ... Additional arguments passed to the rotation function.
    #' @return Rotated draws or updated object.
    internal_rotate = function(target, method = "procrustes", type = "orthogonal",
                               linked_straight = NULL, linked_inverse = NULL,
                               overwrite = NULL, ...) {

      # 上書きしない場合はクローンを作成
      obj <- if (isTRUE(overwrite)) self else self$clone(deep = TRUE)

      f_arr <- obj$fit
      r_arr <- obj$random_fit
      t_arr <- obj$tran_fit
      g_arr <- obj$gq_fit

      v_names_f <- dimnames(f_arr)[[3]]
      v_names_r <- if (!is.null(r_arr)) dimnames(r_arr)[[3]] else character(0)
      v_names_t <- if (!is.null(t_arr)) dimnames(t_arr)[[3]] else character(0)
      v_names_g <- if (!is.null(g_arr)) dimnames(g_arr)[[3]] else character(0)

      # 配列統合の恩恵: fit, random_fit, tran_fit, gq_fit のどこにある変数か自動判定
      get_var_info <- function(vname) {
        pattern <- paste0("^", vname, "\\[")
        idx_f <- grep(pattern, v_names_f)
        if (length(idx_f) > 0)
          return(list(
            loc = "fixed",
            idx = idx_f,
            dim = obj$model$par_list[[vname]]$dim
          ))

        idx_r <- grep(pattern, v_names_r)
        if (length(idx_r) > 0)
          return(list(
            loc = "random",
            idx = idx_r,
            dim = obj$model$par_list[[vname]]$dim
          ))

        idx_t <- grep(pattern, v_names_t)
        if (length(idx_t) > 0)
          return(list(
            loc = "tran",
            idx = idx_t,
            dim = obj$tran_dims[[vname]]
          ))

        idx_g <- grep(pattern, v_names_g)
        if (length(idx_g) > 0)
          return(list(
            loc = "gq",
            idx = idx_g,
            dim = obj$gq_dims[[vname]]
          ))

        stop(paste0("Rotation failed: Variable '", vname, "' not found."))
      }

      t_info <- get_var_info(target)
      if (length(t_info$dim) != 2)
        stop(paste0("Target variable '", target, "' must be a matrix."))
      R_t <- t_info$dim[1]; C_t <- t_info$dim[2]

      # MAP推定値 (最大lp) を取得して基準とする
      lp_mat <- f_arr[, , 1]
      max_idx <- which(lp_mat == max(lp_mat, na.rm = TRUE), arr.ind = TRUE)
      best_iter <- max_idx[1, 1]
      best_chain <- max_idx[1, 2]

      if (t_info$loc == "fixed")
        Y_vec <- f_arr[best_iter, best_chain, t_info$idx]
      else if (t_info$loc == "random")
        Y_vec <- r_arr[best_iter, best_chain, t_info$idx]
      else if (t_info$loc == "tran")
        Y_vec <- t_arr[best_iter, best_chain, t_info$idx]
      else
        Y_vec <- g_arr[best_iter, best_chain, t_info$idx]

      X_map <- matrix(Y_vec, nrow = R_t, ncol = C_t)

      # --- 変換行列の決定 ---
      if (method == "procrustes") {
        Th_straight <- diag(1, C_t)
        Th_inverse <- diag(1, C_t)
      } else {
        if (!requireNamespace("GPArotation", quietly = TRUE))
          stop("GPArotation is required.")
        rot_fn <- tryCatch(
          match.fun(method),
          error = function(e) {
            if (exists(method, where = asNamespace("GPArotation"), mode = "function"))
              return(getFromNamespace(method, "GPArotation"))
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

      for (c in 1:chains) {
        for (i in 1:iter_total) {

          if (t_info$loc == "fixed") X_vec <- f_arr[i, c, t_info$idx]
          else if (t_info$loc == "random") X_vec <- r_arr[i, c, t_info$idx]
          else if (t_info$loc == "tran") X_vec <- t_arr[i, c, t_info$idx]
          else X_vec <- g_arr[i, c, t_info$idx]

          X <- matrix(X_vec, nrow = R_t, ncol = C_t)
          svd_out <- svd(t(X) %*% X_map)
          R_proc <- svd_out$u %*% t(svd_out$v)
          X_rot <- (X %*% R_proc) %*% Th_straight

          if (t_info$loc == "fixed") f_arr[i, c, t_info$idx] <- as.numeric(X_rot)
          else if (t_info$loc == "random") r_arr[i, c, t_info$idx] <- as.numeric(X_rot)
          else if (t_info$loc == "tran") t_arr[i, c, t_info$idx] <- as.numeric(X_rot)
          else g_arr[i, c, t_info$idx] <- as.numeric(X_rot)

          # Linked straight
          if (!is.null(linked_straight)) {
            for (lvar in linked_straight) {
              l_info <- get_var_info(lvar)
              if (length(l_info$dim) != 2 || l_info$dim[2] != C_t) next
              if (l_info$loc == "fixed") {
                Z <- matrix(f_arr[i, c, l_info$idx], nrow = l_info$dim[1], ncol = C_t)
                f_arr[i, c, l_info$idx] <- as.numeric((Z %*% R_proc) %*% Th_straight)
              } else if (l_info$loc == "random") {
                Z <- matrix(r_arr[i, c, l_info$idx], nrow = l_info$dim[1], ncol = C_t)
                r_arr[i, c, l_info$idx] <- as.numeric((Z %*% R_proc) %*% Th_straight)
              } else if (l_info$loc == "tran") {
                Z <- matrix(t_arr[i, c, l_info$idx], nrow = l_info$dim[1], ncol = C_t)
                t_arr[i, c, l_info$idx] <- as.numeric((Z %*% R_proc) %*% Th_straight)
              } else {
                Z <- matrix(g_arr[i, c, l_info$idx], nrow = l_info$dim[1], ncol = C_t)
                g_arr[i, c, l_info$idx] <- as.numeric((Z %*% R_proc) %*% Th_straight)
              }
            }
          }

          # Linked inverse
          if (!is.null(linked_inverse)) {
            for (lvar in linked_inverse) {
              l_info <- get_var_info(lvar)
              if (length(l_info$dim) != 2 || l_info$dim[2] != C_t) next
              if (l_info$loc == "fixed") {
                Z <- matrix(f_arr[i, c, l_info$idx], nrow = l_info$dim[1], ncol = C_t)
                f_arr[i, c, l_info$idx] <- as.numeric((Z %*% R_proc) %*% Th_inverse)
              } else if (l_info$loc == "random") {
                Z <- matrix(r_arr[i, c, l_info$idx], nrow = l_info$dim[1], ncol = C_t)
                r_arr[i, c, l_info$idx] <- as.numeric((Z %*% R_proc) %*% Th_inverse)
              } else if (l_info$loc == "tran") {
                Z <- matrix(t_arr[i, c, l_info$idx], nrow = l_info$dim[1], ncol = C_t)
                t_arr[i, c, l_info$idx] <- as.numeric((Z %*% R_proc) %*% Th_inverse)
              } else {
                Z <- matrix(g_arr[i, c, l_info$idx], nrow = l_info$dim[1], ncol = C_t)
                g_arr[i, c, l_info$idx] <- as.numeric((Z %*% R_proc) %*% Th_inverse)
              }
            }
          }
        }
      }

      obj$fit <- f_arr
      obj$random_fit <- r_arr
      obj$tran_fit <- t_arr
      obj$gq_fit <- g_arr

      # posterior_meanの再計算
      fixed_mean_new <- apply(obj$fit[, , -1, drop = FALSE], 3, mean)
      new_posterior_mean <- obj$posterior_mean
      new_posterior_mean[names(fixed_mean_new)] <- fixed_mean_new
      if (!is.null(obj$random_fit)) {
        random_mean_new <- apply(obj$random_fit, 3, mean)
        new_posterior_mean[names(random_mean_new)] <- random_mean_new
      }
      obj$posterior_mean <- new_posterior_mean

      if (isTRUE(overwrite)) return(invisible(self)) else return(obj)
    },

    #' #' @description Rotate sampled parameters using orthogonal Procrustes.
    #' #' @param target Character string specifying the target variable.
    #' #' @param linked Character vector of variable names to be rotated. Default is NULL.
    #' #' @param overwrite Logical; whether to overwrite the stored draws. Default is TRUE.
    #' #' @param ... Additional arguments.
    #' #' @return Rotated draws or updated object.
    #' rotate = function(target,
    #'                   linked = NULL,
    #'                   overwrite = TRUE,
    #'                   ...) {
    #'   cat("Applying orthogonal Procrustes rotation to VB samples...\n")
    #'   self$internal_rotate(
    #'     target = target,
    #'     method = "procrustes",
    #'     type = "orthogonal",
    #'     linked_straight = linked,
    #'     linked_inverse = NULL,
    #'     overwrite = overwrite,
    #'     ...
    #'   )
    #' },

    #' @description Rotate factor loadings and optional factor scores.
    #' @param loadings Character string specifying the factor loadings variable.
    #' @param scores Character vector specifying the factor scores variable. Default is NULL.
    #' @param method Character string specifying the rotation method. Default is "promax".
    #' @param type Character string specifying the rotation type ("orthogonal" or "oblique"). Default is "oblique".
    #' @param linked_loadings Character vector of linked loading variables. Default is NULL.
    #' @param overwrite Logical; whether to overwrite the stored draws. Default is TRUE.
    #' @param ... Additional arguments.
    #' @return Rotated draws or updated object.
    fa_rotate = function(loadings,
                         scores = NULL,
                         method = "promax",
                         type = "oblique",
                         linked_loadings = NULL,
                         overwrite = TRUE,
                         ...) {
      cat(sprintf("Applying %s rotation to VB samples...\n", method))
      self$internal_rotate(
        target = loadings,
        method = method,
        type = type,
        linked_straight = linked_loadings,
        linked_inverse = scores,
        overwrite = overwrite,
        ...
      )
    }
  )
)
