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
    }
  )
)
