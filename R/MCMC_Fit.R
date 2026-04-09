#' MCMC fit object
#'
#' An R6 class storing posterior samples and related information
#' from MCMC estimation.
#'
#' @param model An `RTMB_Model` object used for estimation.
#' @param fit Posterior draws for model parameters.
#' @param random_fit Posterior draws for random effects, if available.
#' @param tran_fit Posterior draws for transformed parameters, if available.
#' @param tran_dims Dimension information for transformed parameters.
#' @param gq_fit Posterior draws for generated quantities, if available.
#' @param gq_dims Dimension information for generated quantities.
#' @param eps Step size used by the sampler.
#' @param accept Acceptance statistics from sampling.
#' @param treedepth Tree depth used in HMC/NUTS sampling.
#' @param laplace Logical; whether Laplace approximation was used.
#' @param posterior_mean Posterior mean estimates.
#' @param pars Names of parameters to extract or summarize.
#' @param chains Chains to include.
#' @param inc_random Logical; whether to include random effects.
#' @param inc_generate Logical; whether to include generated quantities.
#' @param max_rows Maximum number of rows to print in summaries.
#' @param digits Number of digits to print.
#' @param method Method name.
#' @param seed Random seed.
#' @param tran_fn A function for transformed parameters.
#' @param gq_fn A function for generated quantities.
#' @param target Target variable(s) for rotation or relabeling.
#' @param type Rotation type.
#' @param linked_straight Linked variables transformed in the same direction.
#' @param linked_inverse Linked variables transformed in the opposite direction.
#' @param overwrite Logical; whether to overwrite stored draws.
#' @param linked Linked variables.
#' @param loadings Factor loading variables.
#' @param scores Factor score variables.
#' @param linked_loadings Linked loading variables.
#' @param ... Additional arguments.
#'
#' @field model An `RTMB_Model` object used for estimation.
#' @field fit Posterior draws for model parameters.
#' @field random_fit Posterior draws for random effects.
#' @field transform_fit Posterior draws for transformed parameters.
#' @field transform_dims Dimension information for transformed parameters.
#' @field generate_fit Posterior draws for generated quantities.
#' @field generate_dims Dimension information for generated quantities.
#' @field eps Step size used by the sampler.
#' @field accept Acceptance statistics from sampling.
#' @field treedepth Tree depth used in HMC/NUTS sampling.
#' @field laplace Logical; whether Laplace approximation was used.
#' @field posterior_mean Posterior mean estimates.
#'
MCMC_Fit <- R6::R6Class(
  classname = "mcmc_fit",

  public = list(
    # --- フィールド ---
    model          = NULL, # RTMB_Model のインスタンスへの参照
    fit            = NULL,
    random_fit     = NULL,
    transform_fit  = NULL, # 変換量を保存
    generate_fit   = NULL, # 生成量を保存
    transform_dims = NULL, # 変換量の次元情報を保存
    generate_dims  = NULL, # GQ変数の次元情報を保存
    eps            = NULL,
    accept         = NULL,
    treedepth      = NULL,
    laplace        = NULL,
    posterior_mean = NULL,

    # 1. コンストラクタ
    #' @description Create a new `MCMC_Fit` object.
    initialize = function(model, fit, random_fit, eps, accept, treedepth, laplace, posterior_mean) {
      self$model <- model
      self$fit <- fit
      self$random_fit <- random_fit
      self$eps <- eps
      self$accept <- accept
      self$treedepth <- treedepth
      self$laplace <- laplace
      self$posterior_mean <- posterior_mean
      self$transform_fit <- NULL
      self$transform_dims <- list()
      self$generate_fit <- NULL
      self$generate_dims <- list()
    },

    #' @description Print a brief summary of the fitted object.
    #' @return The object itself, invisibly.
    print = function(...) {
      print(self$summary(...))
      invisible(self)
    },
    #' @description Extract posterior draws for selected parameters.
    #' @param pars Character or numeric vector specifying the names or indices of parameters to extract. If NULL, all available parameters are extracted.
    #' @param chains Numeric vector specifying the chains to extract. If NULL, draws from all chains are returned.
    #' @param inc_random Logical; whether to include random effects in the output. Default is FALSE.
    #' @param inc_tranform Logical; whether to include transformed parameters in the output. Default is TRUE.
    #' @param inc_generate Logical; whether to include generated quantities in the output. Default is TRUE.
    #' @return Posterior draws.
    draws = function(pars = NULL, chains = NULL,
                     inc_random = FALSE, inc_tranform = TRUE, inc_generate = TRUE) {
      out_array <- self$fit

      if (inc_random && !is.null(self$random_fit)) {
        P1 <- dim(out_array)[3]
        P2 <- dim(self$random_fit)[3]
        I <- dim(out_array)[1]
        C <- dim(out_array)[2]
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
        P1 <- dim(out_array)[3];P2 <- dim(self$generate_fit)[3]
        I <- dim(out_array)[1];C <- dim(out_array)[2]
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

      if (is.null(chains)) {
        return(out_array[, , target_idx, drop = FALSE])
      } else {
        return(out_array[, chains, target_idx, drop = FALSE])
      }
    },

    #' @description Summarize posterior draws.
    #' @param pars Character or numeric vector specifying the names or indices of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param max_rows Integer; maximum number of rows to print in the summary table. Default is 10.
    #' @param digits Integer; number of decimal places to print. Default is 2.
    #' @param inc_random Logical; whether to include random effects in the summary. Default is FALSE.
    #' @param inc_transform Logical; whether to include transformed parameters in the summary. Default is TRUE.
    #' @param inc_generate Logical; whether to include generated quantities in the summary. Default is TRUE.
    #' @return A summary object.
    summary = function(pars = NULL, max_rows = 10, digits = 2,
                       inc_random = FALSE, inc_transform = TRUE, inc_generate = TRUE){
      draws_array <- self$draws(pars = pars,
                                chains = NULL,
                                inc_random = inc_random,
                                inc_tranform = inc_transform,
                                inc_generate = inc_generate)

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
        mat_p <- as.matrix(draws_array[, , p])
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
              ess_bulk = NA,
              ess_tail = NA,
              rhat = NA,
              stringsAsFactors = FALSE
            )
          next
        }

        sd_val <- sd(valid_vec)
        if (is.na(sd_val) || sd_val < 1e-10) {
          map_val <- valid_vec[1]; q95 <- c(valid_vec[1], valid_vec[1])
          rhat_val <- NA; ebulk_val <- NA; etail_val <- NA
        } else {
          map_val   <- map_est(valid_vec)
          q95       <- quantile95(valid_vec)
          rhat_val  <- r_hat(mat_p)
          ebulk_val <- ess_bulk(mat_p)
          etail_val <- ess_tail95(mat_p)
        }

        res_list_sum[[i]] <- data.frame(
          variable = param_names[p],
          mean     = mean(valid_vec),
          sd       = sd_val,
          map      = map_val,
          q2.5     = unname(q95[1]),
          q97.5    = unname(q95[2]),
          ess_bulk = ebulk_val,
          ess_tail = etail_val,
          rhat     = rhat_val,
          stringsAsFactors = FALSE
        )
      }
      res_df <- do.call(rbind, res_list_sum)
      class(res_df) <- c("summary_BayesRTMB", "data.frame")
      return(res_df)
    },

    #' @description Transform posterior draws to the unconstrained scale.
    #' @return Posterior draws on the unconstrained scale.
    unconstrain_draws = function() {
      I <- dim(self$fit)[1]
      C <- dim(self$fit)[2]
      orig_pl <- self$model$par_list
      random_flags <- sapply(orig_pl, function(x) isTRUE(x$random))

      target_par_list <- if (self$laplace && any(random_flags)) orig_pl[!random_flags] else orig_pl
      draws_ob <- self$fit[, , -1, drop = FALSE]

      Q <- sum(sapply(target_par_list, function(x) x$unc_length))
      draws_uc <- array(NA, dim = c(I, C, Q))

      for (i in 1:I) {
        for (j in 1:C) {
          con_vec <- draws_ob[i, j, ]
          con_list <- constrained_vector_to_list(con_vec, target_par_list)
          unc_list <- to_unconstrained(con_list, target_par_list)
          draws_uc[i, j, ] <- unlist(unc_list, use.names = FALSE)
        }
      }
      return(matrix(draws_uc, nrow = I * C, ncol = Q))
    },

    #' @description Evaluate log-probability values.
    #' @return Numeric vector of log-probability values.
    log_prob = function() {
      ad_setup <- self$model$build_ad_obj(init = self$posterior_mean, laplace = self$laplace, jacobian_target = "all")
      ad_obj <- ad_setup$ad_obj

      fn <- function(q) {
        val <- tryCatch({ -ad_obj$fn(q) }, error = function(e) -Inf)
        if (is.na(val) || is.nan(val)) return(-Inf)
        return(val)
      }
      return(fn)
    },

    #' @description Estimate the marginal likelihood by bridge sampling.
    #' @return Bridge sampling result.
    bridgesampling = function(method = "warp3", seed = NULL) {

      if (!is.null(seed)) set.seed(seed)

      draws_uc <- self$unconstrain_draws()
      log_prob_fn <- self$log_prob()

      N_total <- nrow(draws_uc)
      Q <- ncol(draws_uc)

      index1 <- seq(1, N_total - 1, by = 2)
      index2 <- seq(2, N_total, by = 2)
      z_post <- draws_uc[index1, , drop = FALSE]
      z_fit  <- draws_uc[index2, , drop = FALSE]

      M1 <- nrow(z_post)
      M2 <- nrow(z_fit)
      M  <- M1 + M2

      S1 <- M1/M
      S2 <- M2/M

      meanz <- apply(z_fit, 2, mean)
      covz <- cov(z_fit)

      if (method == "normal") {
        z_propose <- MASS::mvrnorm(M2, meanz, covz)
        log_propose <- function(z) mvtnorm::dmvnorm(z, meanz, covz, log = TRUE)

        log_L1 <- apply(z_post, 1, log_prob_fn) - log_propose(z_post)
        log_L2 <- apply(z_propose, 1, log_prob_fn) - log_propose(z_propose)

        valid_log_L1 <- log_L1[is.finite(log_L1)]
        if (length(valid_log_L1) == 0) stop("有効な事後確率サンプルがありません。")
        log_Lm <- median(valid_log_L1)

        L1 <- exp(log_L1 - log_Lm)
        L2 <- exp(log_L2 - log_Lm)
        L1[is.na(L1) | is.nan(L1)] <- 0
        L2[is.na(L2) | is.nan(L2)] <- 0

        ml_t <- exp(log_prob_fn(meanz) - log_Lm)
        if (is.na(ml_t) || is.nan(ml_t) || ml_t == 0) ml_t <- 1e-10

      } else if (method == "warp3") {
        L_chol <- t(chol(covz + diag(0.0001, Q)))
        log_det_L <- sum(log(diag(L_chol)))

        eta_post <- t(solve(L_chol, t(sweep(z_post, 2, meanz, "-"))))
        eta_prop <- matrix(rnorm(M2 * Q), nrow = M2, ncol = Q)

        log_propose_std <- function(eta) {
          sum(dnorm(eta, mean = 0, sd = 1, log = TRUE))
        }

        log_target_warp3 <- function(eta) {
          z_plus  <- meanz + as.vector(L_chol %*% eta)
          z_minus <- meanz - as.vector(L_chol %*% eta)

          lp_plus  <- log_prob_fn(z_plus)
          lp_minus <- log_prob_fn(z_minus)

          max_lp <- max(lp_plus, lp_minus)
          if (is.infinite(max_lp)) return(-Inf)

          log_mix <-
            log(0.5) + max_lp + log(exp(lp_plus - max_lp) + exp(lp_minus - max_lp))
          return(log_mix + log_det_L)
        }

        log_L1 <- apply(eta_post, 1, log_target_warp3) - apply(eta_post, 1, log_propose_std)
        log_L2 <- apply(eta_prop, 1, log_target_warp3) - apply(eta_prop, 1, log_propose_std)

        valid_log_L1 <- log_L1[is.finite(log_L1)]
        if (length(valid_log_L1) == 0) stop("有効な事後確率サンプルがありません。")
        log_Lm <- median(valid_log_L1)

        L1 <- exp(log_L1 - log_Lm)
        L2 <- exp(log_L2 - log_Lm)
        L1[is.na(L1) | is.nan(L1)] <- 0
        L2[is.na(L2) | is.nan(L2)] <- 0

        ml_t <- exp(log_target_warp3(rep(0, Q)) - log_Lm)
        if (is.na(ml_t) || is.nan(ml_t) || ml_t == 0) ml_t <- 1e-10

      } else {
        stop("method引数は 'normal' または 'warp3' を指定してください。")
      }

      Trial <- 100
      for(t in 1:Trial){
        ml_t_old <- ml_t
        bunbo <- 0
        bunshi <- 0
        for(m in 1:M1) bunbo <- bunbo + 1/(S1*L1[m] + S2*ml_t)
        for(m in 1:M2) bunshi <- bunshi + L2[m] / (S1*L2[m] + S2*ml_t)

        ml_t <- (bunshi/M2) / (bunbo/M1)
        cat(paste0("iteration: ", t, "\n"))

        if (is.na(ml_t) || is.nan(ml_t) || ml_t <= 0) {
          warning("計算中に ml_t が 0 または NA になりました。")
          break
        }
        if(abs(log(ml_t) - log(ml_t_old)) < 0.000001){
          break
        }
      }

      logml.bs <- as.numeric(log(ml_t) + log_Lm)
      return(logml.bs)
    },

    #' @description Compute transformed parameters from posterior draws.
    #' @return Transformed parameter draws.
    transformed_draws = function(tran_fn = NULL) {
      all_draws <- self$draws(
        inc_random = TRUE,
        inc_transform = FALSE,
        inc_generate = FALSE
      )

      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2]

      wrapper_tran_fn <- function(dat, param) {
        res <- list()

        # # CF_corr から相関行列を自動生成
        # for (name in names(self$model$par_list)) {
        #   if (self$model$par_list[[name]]$type == "CF_corr") {
        #     mat_name <- if (grepl("^CF_", name)) sub("^CF_", "", name) else paste0(name, "_corr")
        #     res[[mat_name]] <- param[[name]] %*% t(param[[name]])
        #   }
        # }

        # ユーザー定義 transformed_parameters を追加
        if (!is.null(tran_fn)) {
          user_res <- tran_fn(dat, param)

          if (is.null(user_res)) {
            user_res <- list()
          } else if (!is.list(user_res)) {
            stop("`transformed_parameters` は list を返す必要があります。", call. = FALSE)
          }

          # 名前の重複チェック
          dup_names <- intersect(names(res), names(user_res))
          if (length(dup_names) > 0) {
            stop(
              sprintf(
                "transformed parameters 名が重複しています: %s",
                paste(dup_names, collapse = ", ")
              ),
              call. = FALSE
            )
          }

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

        names_def <- self$model$par_names[[name]]
        flat_nms <- generate_flat_names(name, dim_val, names_def)
        tran_names <- c(tran_names, flat_nms)
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
    #' @return Generated quantities draws.
    generated_quantities = function(gq_fn = NULL) {
      all_draws <- self$draws(
        inc_random = TRUE,
        inc_transform = FALSE,
        inc_generate = FALSE
      )
      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2]

      if (is.null(gq_fn)) return(invisible(self))

      test_para <- constrained_vector_to_list(all_draws[1, 1, -1], self$model$par_list)

      if (!is.null(self$model$transform)) {
        tran_res <- self$model$transform(self$model$data, test_para)
        if (is.list(tran_res)) {
          test_para <- c(test_para, tran_res)
        }
      }

      test_gq <- gq_fn(self$model$data, test_para)

      if (is.null(test_gq) || length(test_gq) == 0) return(invisible(self))
      if (!is.list(test_gq)) {
        stop("`generate` は list を返す必要があります。", call. = FALSE)
      }

      gq_names <- character(0)
      self$gq_dims <- list()

      for (name in names(test_gq)) {
        val <- test_gq[[name]]
        len <- length(val)
        dim_val <- dim(val)
        if (is.null(dim_val)) dim_val <- len
        self$gq_dims[[name]] <- dim_val

        names_def <- self$model$par_names[[name]]
        flat_nms <- generate_flat_names(name, dim_val, names_def)
        gq_names <- c(gq_names, flat_nms)
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

          if (!is.null(self$model$transform)) {
            tran_res <- self$model$transform(self$model$data, p_list)
            if (is.list(tran_res)) {
              p_list <- c(p_list, tran_res)
            }
          }

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

          if (t_info$loc == "fixed") X_vec <- f_arr[i, c, t_info$idx]
          else if (t_info$loc == "random") X_vec <- r_arr[i, c, t_info$idx]
          else if (t_info$loc == "transform") X_vec <- t_arr[i, c, t_info$idx]
          else X_vec <- g_arr[i, c, t_info$idx]

          X <- matrix(X_vec, nrow = R_t, ncol = C_t)
          svd_out <- svd(t(X) %*% X_map)
          R_proc <- svd_out$u %*% t(svd_out$v)
          X_rot <- as.numeric((X %*% R_proc) %*% Th_straight)

          if (isTRUE(overwrite)) {
            if (t_info$loc == "fixed") f_arr[i, c, t_info$idx] <- X_rot
            else if (t_info$loc == "random") r_arr[i, c, t_info$idx] <- X_rot
            else if (t_info$loc == "transform") t_arr[i, c, t_info$idx] <- X_rot
            else g_arr[i, c, t_info$idx] <- X_rot
          } else {
            rot_arr[i, c, new_idx_map[[target]]] <- X_rot
          }

          if (!is.null(linked_straight)) {
            for (lvar in linked_straight) {
              l_info <- get_var_info(lvar)
              if (length(l_info$dim) != 2 || l_info$dim[2] != C_t) next

              if (l_info$loc == "fixed") Z_vec <- f_arr[i, c, l_info$idx]
              else if (l_info$loc == "random") Z_vec <- r_arr[i, c, l_info$idx]
              else if (l_info$loc == "transform") Z_vec <- t_arr[i, c, l_info$idx]
              else Z_vec <- g_arr[i, c, l_info$idx]

              Z <- matrix(Z_vec, nrow = l_info$dim[1], ncol = C_t)
              Z_rot <- as.numeric((Z %*% R_proc) %*% Th_straight)

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

          if (!is.null(linked_inverse)) {
            for (lvar in linked_inverse) {
              l_info <- get_var_info(lvar)
              if (length(l_info$dim) != 2 || l_info$dim[2] != C_t) next

              if (l_info$loc == "fixed") Z_vec <- f_arr[i, c, l_info$idx]
              else if (l_info$loc == "random") Z_vec <- r_arr[i, c, l_info$idx]
              else if (l_info$loc == "transform") Z_vec <- t_arr[i, c, l_info$idx]
              else Z_vec <- g_arr[i, c, l_info$idx]

              Z <- matrix(Z_vec, nrow = l_info$dim[1], ncol = C_t)
              Z_rot <- as.numeric((Z %*% R_proc) %*% Th_inverse)

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

    #' @description Rotate sampled parameters.
    #' @param target Character string specifying the target variable to base the rotation on.
    #' @param linked Character vector of variable names to be rotated in the same direction.
    #' @param overwrite Logical; whether to overwrite the stored draws. If FALSE, adds to generated quantities. Default is FALSE.
    #' @param suffix Character string to append to the rotated variable names when overwrite is FALSE. Default is "rot".
    #' @param ... Additional arguments passed to the rotation function.
    #' @return The updated object invisibly.
    rotate = function(target,
                      linked = NULL,
                      overwrite = FALSE,
                      suffix = "rot",
                      ...) {
      cat(if(isTRUE(overwrite)) "Applying orthogonal Procrustes rotation (Overwriting)...\n" else sprintf("Applying orthogonal Procrustes rotation (Saving to generate as _%s)...\n", suffix))
      self$internal_rotate(
        target = target,
        method = "procrustes",
        type = "orthogonal",
        linked_straight = linked,
        linked_inverse = NULL,
        overwrite = overwrite,
        suffix = suffix,
        ...
      )
    },

    #' @description Rotate factor loadings and optional factor scores.
    #' @param loadings Character string specifying the factor loadings variable.
    #' @param scores Character vector specifying the factor scores variable.
    #' @param method Character string specifying the rotation method.
    #' @param type Character string specifying the rotation type.
    #' @param linked_loadings Character vector of linked loading variables.
    #' @param overwrite Logical; whether to overwrite the stored draws. If FALSE, adds to generated quantities. Default is FALSE.
    #' @param suffix Character string to append to the rotated variable names when overwrite is FALSE. Default is the method name.
    #' @param ... Additional arguments passed to the rotation function.
    #' @return The updated object invisibly.
    fa_rotate = function(loadings,
                         scores = NULL,
                         method = "promax",
                         type = "oblique",
                         linked_loadings = NULL,
                         overwrite = FALSE,
                         suffix = method,
                         ...) {
      cat(if(isTRUE(overwrite)) sprintf("Applying %s rotation (Overwriting)...\n", method) else sprintf("Applying %s rotation (Saving to generate as _%s)...\n", method, suffix))
      self$internal_rotate(
        target = loadings,
        method = method,
        type = type,
        linked_straight = linked_loadings,
        linked_inverse = scores,
        overwrite = overwrite,
        suffix = suffix,
        ...
      )
    },

    #' @description Resolve label switching in posterior draws.
    #' @param target Character string specifying the target variable to base the relabeling on.
    #' @param linked Character vector of variable names to be relabeled in the same order as the target. Default is NULL.
    #' @param overwrite Logical; whether to overwrite the stored draws in the current object. Default is TRUE.
    #' @param scalar_fns A named list of functions to apply to scalar variables for relabeling. Default is an empty list.
    #' @return Relabeled draws or updated object.
    resolve_switching = function(target, linked = NULL, overwrite = TRUE, scalar_fns = list()) {
      cat(sprintf("Resolving label switching based on '%s'...\n", target))

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

      # 変数の場所と次元を特定するヘルパー
      get_var_info <- function(vname) {
        pattern <- paste0("^", vname, "\\[")

        idx_f <- grep(pattern, v_names_f)
        if (length(idx_f) > 0) {
          return(list(
            loc = "fixed",
            idx = idx_f,
            dim = obj$model$par_list[[vname]]$dim
          ))
        }

        idx_r <- grep(pattern, v_names_r)
        if (length(idx_r) > 0) {
          return(list(
            loc = "random",
            idx = idx_r,
            dim = obj$model$par_list[[vname]]$dim
          ))
        }

        idx_t <- grep(pattern, v_names_t)
        if (length(idx_t) > 0) {
          return(list(
            loc = "tran",
            idx = idx_t,
            dim = obj$tran_dims[[vname]]
          ))
        }

        idx_g <- grep(pattern, v_names_g)
        if (length(idx_g) > 0) {
          return(list(
            loc = "gq",
            idx = idx_g,
            dim = obj$gq_dims[[vname]]
          ))
        }

        # スカラー変数も許す
        idx_f0 <- which(v_names_f == vname)
        if (length(idx_f0) == 1) {
          return(list(loc = "fixed_scalar", idx = idx_f0, dim = 1))
        }

        idx_r0 <- which(v_names_r == vname)
        if (length(idx_r0) == 1) {
          return(list(loc = "random_scalar", idx = idx_r0, dim = 1))
        }

        idx_t0 <- which(v_names_t == vname)
        if (length(idx_t0) == 1) {
          return(list(loc = "tran_scalar", idx = idx_t0, dim = 1))
        }

        idx_g0 <- which(v_names_g == vname)
        if (length(idx_g0) == 1) {
          return(list(loc = "gq_scalar", idx = idx_g0, dim = 1))
        }

        stop(paste0("Label switching failed: Variable '", vname, "' not found."))
      }

      # ターゲット取得
      t_info <- get_var_info(target)

      # 並べ替え対象はベクトルのみ許可
      if (grepl("scalar", t_info$loc)) {
        stop(paste0("Target variable '", target, "' must be a vector, not a scalar."))
      }

      if (length(t_info$dim) != 1) {
        stop(paste0("Target variable '", target, "' must be a vector."))
      }

      K <- t_info$dim[1]
      if (length(t_info$idx) != K) {
        stop(paste0("Target variable '", target, "' has inconsistent dimension information."))
      }

      linked_info_list <- list()
      if (!is.null(linked)) {
        for (lvar in linked) {
          l_info <- get_var_info(lvar)

          if (grepl("scalar", l_info$loc)) {
            linked_info_list[[lvar]] <- l_info
          } else {
            # ベクトルのみ連動可
            if (length(l_info$dim) != 1) {
              warning(sprintf(
                "Variable '%s' is not a vector. Skipping.",
                lvar
              ))
            } else if (l_info$dim[1] != K) {
              warning(sprintf(
                "Length of '%s' (%d) != target '%s' (%d). Skipping.",
                lvar, l_info$dim[1], target, K
              ))
            } else {
              linked_info_list[[lvar]] <- l_info
            }
          }
        }
      }

      iter_total   <- dim(f_arr)[1]
      chains_total <- dim(f_arr)[2]

      get_values <- function(arr, i, c, idx) {
        arr[i, c, idx]
      }

      set_values <- function(arr, i, c, idx, value) {
        arr[i, c, idx] <- value
        arr
      }

      for (c in seq_len(chains_total)) {
        for (i in seq_len(iter_total)) {

          # ターゲット値取得
          if (t_info$loc == "fixed") {
            t_vals <- get_values(f_arr, i, c, t_info$idx)
          } else if (t_info$loc == "random") {
            t_vals <- get_values(r_arr, i, c, t_info$idx)
          } else if (t_info$loc == "tran") {
            t_vals <- get_values(t_arr, i, c, t_info$idx)
          } else if (t_info$loc == "gq") {
            t_vals <- get_values(g_arr, i, c, t_info$idx)
          } else {
            stop("Target must be a vector variable.")
          }

          ord <- order(t_vals)

          if (any(ord != seq_len(K))) {
            # target 自体を並べ替え
            if (t_info$loc == "fixed") {
              f_arr <- set_values(f_arr, i, c, t_info$idx, t_vals[ord])
            } else if (t_info$loc == "random") {
              r_arr <- set_values(r_arr, i, c, t_info$idx, t_vals[ord])
            } else if (t_info$loc == "tran") {
              t_arr <- set_values(t_arr, i, c, t_info$idx, t_vals[ord])
            } else if (t_info$loc == "gq") {
              g_arr <- set_values(g_arr, i, c, t_info$idx, t_vals[ord])
            }

            # linked を同じ順序で並べ替え
            for (lvar in names(linked_info_list)) {
              l_info <- linked_info_list[[lvar]]

              if (l_info$loc == "fixed") {
                f_arr <- set_values(f_arr, i, c, l_info$idx, f_arr[i, c, l_info$idx][ord])

              } else if (l_info$loc == "random") {
                r_arr <- set_values(r_arr, i, c, l_info$idx, r_arr[i, c, l_info$idx][ord])

              } else if (l_info$loc == "tran") {
                t_arr <- set_values(t_arr, i, c, l_info$idx, t_arr[i, c, l_info$idx][ord])

              } else if (l_info$loc == "gq") {
                g_arr <- set_values(g_arr, i, c, l_info$idx, g_arr[i, c, l_info$idx][ord])

              } else if (K == 2) {
                # 2群ラベルスイッチ用にスカラー補助変数の処理を適用
                s_fn <- function(x) 1 - x
                if (!is.null(scalar_fns) && lvar %in% names(scalar_fns)) {
                  s_fn <- scalar_fns[[lvar]]
                }

                if (l_info$loc == "fixed_scalar") {
                  f_arr[i, c, l_info$idx] <- s_fn(f_arr[i, c, l_info$idx])
                } else if (l_info$loc == "random_scalar") {
                  r_arr[i, c, l_info$idx] <- s_fn(r_arr[i, c, l_info$idx])
                } else if (l_info$loc == "tran_scalar") {
                  t_arr[i, c, l_info$idx] <- s_fn(t_arr[i, c, l_info$idx])
                } else if (l_info$loc == "gq_scalar") {
                  g_arr[i, c, l_info$idx] <- s_fn(g_arr[i, c, l_info$idx])
                }
              }
            }
          }
        }
      }

      obj$fit <- f_arr
      obj$random_fit <- r_arr
      obj$tran_fit <- t_arr
      obj$gq_fit <- g_arr

      # posterior_mean の再計算
      fixed_mean_new <- apply(obj$fit[, , -1, drop = FALSE], 3, mean)
      new_posterior_mean <- obj$posterior_mean
      new_posterior_mean[names(fixed_mean_new)] <- fixed_mean_new

      if (!is.null(obj$random_fit)) {
        random_mean_new <- apply(obj$random_fit, 3, mean)
        new_posterior_mean[names(random_mean_new)] <- random_mean_new
      }

      obj$posterior_mean <- new_posterior_mean

      if (isTRUE(overwrite)) invisible(self) else obj
    }
  )
)
