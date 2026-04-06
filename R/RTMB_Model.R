#' RTMB model object
#'
#' An R6 class representing a Bayesian model built with RTMB.
#' This class stores model components and provides methods for
#' building the automatic differentiation object, optimizing the
#' posterior, and drawing posterior samples.
#'
#' @param data A list of observed data used in the model.
#' @param par_list A list defining model parameters and their dimensions.
#' @param log_prob A user-supplied function that returns the log-probability.
#' @param transform An optional function for transformed parameters.
#' @param generate An optional function for generated quantities.
#' @param init Optional initial values for parameters.
#' @param laplace Logical; whether to use Laplace approximation for
#'   random effects.
#' @param include_jacobian Logical; whether to include Jacobian
#'   adjustments in the objective function.
#' @param control A list of control settings passed to the optimizer.
#' @param sampling A character string specifying the sampler to use.
#' @param warmup Number of warmup iterations.
#' @param chains Number of MCMC chains.
#' @param thin Thinning interval.
#' @param seed Random seed.
#' @param delta Target acceptance rate for HMC/NUTS.
#' @param max_treedepth Maximum tree depth for HMC/NUTS.
#' @param parallel Logical; whether to run chains in parallel.
#' @param ... Additional arguments.
#'
#' @field data A list of observed data.
#' @field par_list A list defining model parameters.
#' @field log_prob A user-supplied log-probability function.
#' @field transform An optional transformed parameters function.
#' @field generate An optional generated quantities function.
#' @field pl_full Full parameter information used internally.
#' @import RTMB
#' @export
RTMB_Model <- R6::R6Class(
  classname = "RTMB_Model",

  public = list(
    data       = NULL,
    par_list   = NULL,
    log_prob   = NULL,
    transform  = NULL,
    generate   = NULL,
    pl_full    = NULL,

    # 1. コンストラクタ
    #' @description Create a new `RTMB_Model` object.
    initialize = function(data, par_list, log_prob, transform = NULL, generate = NULL) {
      self$data <- data
      self$par_list <- lapply(par_list, function(x) {
        if (typeof(x) %in% c("language", "symbol")) {
          eval(x, envir = self$data, enclos = parent.frame())
        } else {
          x
        }
      })
      self$log_prob   <- log_prob
      self$transform  <- transform
      self$generate   <- generate

      self$pl_full <- parse_parameters(self$par_list)
      names(self$pl_full$init) <- self$pl_full$names

      init_vec <- generate_random_init(self$pl_full, self$par_list, range = 2)
      test_para <- constrained_vector_to_list(init_vec, self$par_list)

      test_val <- tryCatch({
        if (!is.null(self$transform)) {
          test_tran <- self$transform(self$data, test_para)
          test_para <- c(test_para, test_tran)
        }
        self$log_prob(self$data, test_para)
      }, error = function(e) {
        stop("R関数(log_prob)のテスト実行でエラーが発生しました。\n[エラー]: "
             , e$message, call. = FALSE)
      })

      if (!is.numeric(test_val) || length(test_val) != 1) {
        stop("log_prob は長さ1の数値（スカラー）を返す必要があります。")
      }

      if (!is.null(self$transform)) {
        test_tran <- tryCatch({
          self$transform(self$data, test_para)
        }, error = function(e) {
          stop("R関数(transform)のテスト実行でエラーが発生しました。\n[エラー]: ",
               e$message, call. = FALSE)
        })

        if (!is.null(test_tran)) {
          if (!is.list(test_tran)) {
            stop("transform は list を返す必要があります。", call. = FALSE)
          }
          if (is.null(names(test_tran)) || any(names(test_tran) == "")) {
            stop("transform は名前付き list を返す必要があります。", call. = FALSE)
          }
        }
      }

      cat("Checking RTMB setup...\n")
      test_ad <- self$build_ad_obj(init = init_vec, laplace = FALSE, jacobian_target = "all")
      test_gr <- tryCatch(test_ad$ad_obj$gr(test_ad$ad_obj$par), error = function(e) e)
      if (inherits(test_gr, "error")) {
        stop("MakeADFun の勾配計算でエラーが発生しました。\n[エラー]: ",
             test_gr$message, call. = FALSE)
      }
    },

    # 2. ADオブジェクト生成ファクトリ
    # ヤコビアンを追加するかどうかを引数で制御
    #' @description Build the RTMB automatic differentiation object.
    #' @param init Optional numeric vector or list of initial values for the parameters. Default is NULL.
    #' @param laplace Logical; whether to use Laplace approximation to marginalize random effects. Default is FALSE.
    #' @param jacobian_target Character string specifying which parameters to apply Jacobian adjustments to (e.g., "all", "random", or "none"). Default is "all".
    #' @return An RTMB objective object..
    build_ad_obj = function(init = NULL, laplace = FALSE, jacobian_target = "all") {
      random_effs <-
        names(self$par_list)[sapply(self$par_list, function(x) isTRUE(x$random))]
      use_random <- if (laplace && length(random_effs) > 0) random_effs else NULL

      current_init <- if (is.null(init))
        generate_random_init(self$pl_full, self$par_list, range = 2)
      else
        init

      current_init_list <- constrained_vector_to_list(current_init, self$par_list)
      init_unc_list <- to_unconstrained(current_init_list, self$par_list)

      pl_full_local    <- self$pl_full
      par_list_local   <- self$par_list

      log_prob_local   <- self$log_prob
      data_local       <- self$data

      f_ad <- function(y_unc_list) {
        para <- to_constrained(y_unc_list, par_list_local)
        if (!is.null(self$transform)) {
          tran_res <- self$transform(data_local, para)
          para <- c(para, tran_res)
        }
        lp <- log_prob_local(data_local, para)

        # ヤコビアンの対象を条件分岐
        if (jacobian_target == "all") {
          lj <- calc_log_jacobian(y_unc_list, par_list_local, only_random = FALSE)
          return(-(lp + lj))
        } else if (jacobian_target == "random") {
          lj <- calc_log_jacobian(y_unc_list, par_list_local, only_random = TRUE)
          return(-(lp + lj))
        } else {
          return(-lp)
        }
      }

      ad_obj <- tryCatch({
        RTMB::MakeADFun(
          func = f_ad,
          parameters = init_unc_list,
          random = use_random,
          silent = TRUE,
          #inner.control = list(
          #  smartsearch = FALSE, # 無駄な直線探索を省略
          #  tol10 = 1e-4,        # 収束判定の許容誤差を緩和（デフォルトより緩く）
          #  iter.max = 20        # 最大イテレーション数を制限
          #)
        )
      }, error = function(e) {
        stop("MakeADFun のセットアップに失敗しました。\n[エラー]: "
             , e$message, call. = FALSE)
      })

      return(list(
        ad_obj = ad_obj,
        pl_full = pl_full_local,
        use_random = use_random
      ))
    },


    # 3. 最尤法 / MAP推定 メソッド
    #' @description Optimize the posterior or marginal posterior.
    #' @return A fitted `MAP_Fit` object.
    optimize = function(laplace = TRUE, init = NULL, control = list()) {
      cat("Starting optimization...\n")

      # MAP推定ではヤコビアンは不要
      jac_target <- if (laplace) "random" else "none"
      ad_setup <- self$build_ad_obj(init = init, laplace = laplace, jacobian_target = jac_target)
      ad_obj <- ad_setup$ad_obj

      opt <- nlminb(
        start = ad_obj$par,
        objective = ad_obj$fn,
        gradient = ad_obj$gr,
        control = control
      )

      ad_obj$fn(opt$par)

      sd_rep <- tryCatch(RTMB::sdreport(ad_obj), error = function(e) NULL)
      if (!is.null(ad_obj$env$last.par.best)) {
        unc_est_vec <- ad_obj$env$last.par.best
      } else {
        unc_est_vec <- ad_obj$env$last.par
      }

      unc_se_vec  <- rep(NA, length(unc_est_vec))

      if (!is.null(sd_rep)) {
        idx_ran <- ad_obj$env$random
        smry_fix <- summary(sd_rep, select = "fixed")
        if (laplace && length(idx_ran) > 0) {
          smry_ran <- summary(sd_rep, select = "random")
          unc_se_vec[-idx_ran] <- smry_fix[, "Std. Error"]
          unc_se_vec[idx_ran]  <- smry_ran[, "Std. Error"]
        } else {
          unc_se_vec <- smry_fix[, "Std. Error"]
        }
      }

      unc_est_list <- unconstrained_vector_to_list(unc_est_vec, self$par_list)
      unc_se_list  <- unconstrained_vector_to_list(unc_se_vec, self$par_list)

      con_est_list <- to_constrained(unc_est_list, self$par_list)
      con_se_list  <- list()
      con_lower_list <- list()
      con_upper_list <- list()

      z_95 <- qnorm(0.975)
      eps_diff <- 1e-5

      for (name in names(self$par_list)) {
        p_info <- self$par_list[[name]]

        u_val <- unc_est_list[[name]]
        u_se  <- unc_se_list[[name]]
        c_val <- con_est_list[[name]]

        L_u <- p_info$unc_length
        L_c <- p_info$length

        transform_single <- function(u) {
          tmp_list <- unc_est_list
          tmp_list[[name]] <- u
          tmp_con <- to_constrained(tmp_list, self$par_list)
          return(as.numeric(tmp_con[[name]]))
        }

        J <- matrix(0, nrow = L_c, ncol = L_u)
        for (i in 1:L_u) {
          u_tmp <- u_val
          u_tmp[i] <- u_tmp[i] + eps_diff
          c_tmp <- transform_single(u_tmp)
          J[, i] <- (c_tmp - as.numeric(c_val)) / eps_diff
        }

        c_se <- numeric(L_c)
        for (j in 1:L_c) {
          c_se[j] <- sqrt(sum((J[j, ] * u_se)^2))
        }

        if (length(p_info$dim) > 1) dim(c_se) <- p_info$dim
        con_se_list[[name]] <- c_se

        if (p_info$bounds == "corr_matrix") {
          c_low <- rep(NA, L_c)
          c_up  <- rep(NA, L_c)
        } else {
          u_low <- u_val - z_95 * u_se
          u_up  <- u_val + z_95 * u_se
          c_low <- transform_single(u_low)
          c_up  <- transform_single(u_up)

          c_low_final <- pmin(c_low, c_up)
          c_up_final  <- pmax(c_low, c_up)
          c_low <- c_low_final
          c_up  <- c_up_final
        }

        if (length(p_info$dim) > 1) {
          dim(c_low) <- p_info$dim
          dim(c_up)  <- p_info$dim
        }
        con_lower_list[[name]] <- c_low
        con_upper_list[[name]] <- c_up
      }

      build_df <- function(target_random = FALSE) {
        names_vec <- c()
        est_vec <- c()
        se_vec <- c()
        low_vec <- c()
        up_vec <- c()

        for (name in names(self$par_list)) {
          p_info <- self$par_list[[name]]
          if (p_info$random == target_random) {

            # --- Omega (相関行列) を計算して表に追加 ---
            if (p_info$type == "CF_corr") {
              L_est <- con_est_list[[name]]
              Omega_est <- L_est %*% t(L_est)
              R <- nrow(Omega_est)
              display_name <- sub("^CF_", "", name)
              if (display_name == name) display_name <- paste0(name, "_corr")

              u_val <- unc_est_list[[name]]
              u_se  <- unc_se_list[[name]]

              get_omega_elem <- function(u_vec_name, row, col) {
                tmp_unc_list <- unc_est_list
                tmp_unc_list[[name]] <- u_vec_name
                tmp_con <- to_constrained(tmp_unc_list, self$par_list)
                Om_tmp <- tmp_con[[name]] %*% t(tmp_con[[name]])
                return(as.numeric(Om_tmp[row, col]))
              }

              for (r in 1:R) {
                for (c in 1:R) {
                  est_val <- Omega_est[r, c]
                  grad <- numeric(length(u_val))
                  for (k in seq_along(u_val)) {
                    u_step <- u_val; u_step[k] <- u_step[k] + 1e-5
                    grad[k] <- (get_omega_elem(u_step, r, c) - est_val) / 1e-5
                  }
                  omega_se <- sqrt(sum((grad * u_se)^2))

                  names_vec <- c(names_vec, paste0(display_name, "[", r, ",", c, "]"))
                  est_vec   <- c(est_vec, est_val)
                  se_vec    <- c(se_vec,  omega_se)

                  low_val <- est_val - 1.96 * omega_se
                  up_val  <- est_val + 1.96 * omega_se

                  if (r != c) {
                    low_val <- max(-1, min(1, low_val))
                    up_val  <- max(-1, min(1, up_val))
                  } else {
                    omega_se <- 0; low_val <- 1; up_val <- 1
                  }

                  low_vec   <- c(low_vec, low_val)
                  up_vec    <- c(up_vec,  up_val)
                }
              }
            }

            # --- 本来のパラメータ処理 ---
            len <- p_info$length
            f_names <- character(len)

            if (len == 1) {
              f_names <- name
            } else {
              if (length(p_info$dim) > 1) {
                grid <- expand.grid(lapply(p_info$dim, seq_len))
                indices <- apply(grid, 1, paste, collapse = ",")
                f_names <- paste0(name, "[", indices, "]")
              } else {
                f_names <- paste0(name, "[", 1:len, "]")
              }
            }

            names_vec <- c(names_vec, f_names)
            est_vec <- c(est_vec, as.numeric(con_est_list[[name]]))
            se_vec  <- c(se_vec,  as.numeric(con_se_list[[name]]))
            low_vec <- c(low_vec, as.numeric(con_lower_list[[name]]))
            up_vec  <- c(up_vec,  as.numeric(con_upper_list[[name]]))
          }
        }

        if (length(names_vec) == 0) return(NULL)

        df <- data.frame(
          Estimate     = est_vec,
          `Std. Error` = se_vec,
          `Lower 95%`  = low_vec,
          `Upper 95%`  = up_vec,
          row.names    = names_vec,
          check.names  = FALSE
        )
        return(df)
      }

      df_fixed <- build_df(target_random = FALSE)
      df_random <- if (laplace) build_df(target_random = TRUE) else NULL

      con_est_vec <- unlist(con_est_list, use.names = FALSE)

      log_ml <- NA
      if (!is.null(sd_rep) && !is.null(sd_rep$cov.fixed)) {
        D <- length(opt$par)
        cov_mat <- sd_rep$cov.fixed

        eig <- tryCatch(eigen(cov_mat, symmetric = TRUE), error = function(e) NULL)

        if (!is.null(eig)) {
          vals <- eig$values
          if (any(vals <= 0)) {
            vals[vals <= 0] <- 1e-10
          }
          log_det_cov <- sum(log(vals))
          log_ml <- - opt$objective + (D / 2) * log(2 * pi) + 0.5 * log_det_cov
        }
      }

      # --- MAP_Fit インスタンスを返す (後で定義します) ---
      res_obj <- MAP_Fit$new(
        par_vec = con_est_vec,
        par = con_est_list,
        objective = opt$objective,
        log_ml = log_ml,
        convergence = opt$convergence,
        sd_rep = sd_rep,
        df_fixed = df_fixed,
        random_effects = df_random
      )

      return(res_obj)
    },


    # 4. MCMC NUTS サンプリング メソッド
    #' @description Draw posterior samples from the model.
    #' @return A fitted `MCMC_Fit` object.
    sample = function(sampling=1000, warmup=1000, chains=4,
                      thin=1, seed=sample.int(1e6,1),
                      delta=0.8, max_treedepth = 10,
                      parallel = TRUE, laplace = FALSE,
                      init = NULL) {

      if (!is.null(init)) init <- as.numeric(init)
      set.seed(seed)
      orig_pl <- self$par_list

      random_flags <- sapply(orig_pl, function(x) isTRUE(x$random))

      if (laplace && any(random_flags)) {
        pl_fixed  <- parse_parameters(orig_pl[!random_flags])
        pl_random <- parse_parameters(orig_pl[random_flags])

        fixed_idx  <- which(self$pl_full$names %in% pl_fixed$names)
        random_idx <- which(self$pl_full$names %in% pl_random$names)
      } else {
        pl_fixed   <- self$pl_full
        pl_random  <- NULL
        fixed_idx  <- 1:length(self$pl_full$names)
        random_idx <- integer(0)
      }

      P_fixed <- length(pl_fixed$names)
      P_random <- if(!is.null(pl_random)) length(pl_random$names) else 0

      run_chain <- function(c, p_callback = NULL) {

        if (!is.null(init)) {
          init_list <- constrained_vector_to_list(init, self$par_list)
          unc_init_list <- to_unconstrained(init_list, self$par_list)
          unc_init_vec <- unlist(unc_init_list, use.names = FALSE)

          unc_init_vec <- unc_init_vec + rnorm(length(unc_init_vec), mean = 0, sd = 0.1)
          unc_init_list_new <- unconstrained_vector_to_list(unc_init_vec, self$par_list)
          init_full_list <- to_constrained(unc_init_list_new, self$par_list)
          init_full <- unlist(init_full_list, use.names = FALSE)
        } else {
          init_full <- generate_random_init(self$pl_full, self$par_list, range = 2)
        }

        # MCMCではヤコビアンが必須
        ad_setup <- self$build_ad_obj(init = init_full, laplace = laplace, jacobian_target = "all")
        ad_obj <- ad_setup$ad_obj

        #if (laplace) {
        #  # Laplace近似が有効な場合は IMH 法を実行
        #  res <- IMH_method(
        #    model = ad_obj,
        #    sampling = sampling,
        #    warmup = warmup,
        #    chain = c,
        #    update_progress = p_callback,
        #    df = 4
        #  )
        #} else {
          # 通常は NUTS を実行
          res <- NUTS_method(
            model = ad_obj,
            sampling = sampling,
            warmup = warmup,
            delta = delta,
            max_treedepth = max_treedepth,
            chain = c,
            update_progress = p_callback,
            laplace = laplace
          )
        #}

        P_all_true <- length(self$pl_full$names)
        iter <- sampling + warmup
        para_final <- array(NA, dim = c(iter, P_all_true))

        for (i in 1:iter) {
          x_in <- as.numeric(res$para_fixed[i, ])

          if (laplace && length(ad_obj$env$random) > 0) {
            ad_obj$fn(x_in)
            para_list <- ad_obj$env$parList()
          } else {
            para_list <- ad_obj$env$parList(x = x_in)
          }

          con_list <- to_constrained(para_list, self$par_list)
          para_final[i, ] <- unlist(con_list, use.names = FALSE)
        }
        res$para <- para_final
        res$para_full <- NULL
        return(res)
      }

      results_list <- list()
      if (parallel) {
        future::plan(future::multisession, workers = chains)
        cat(paste0("並列サンプリングを開始します (chains = ", chains, ")...\n"))
        iter <- sampling + warmup
        total_updates <- chains * floor(iter / 100)

        progressr::with_progress({
          p <- progressr::progressor(steps = total_updates)
          results_list <- future.apply::future_lapply(1:chains, function(c) {
            run_chain(c, p_callback = function() p())
          }, future.seed = TRUE,
          future.globals = c(
            "unconstrained_vector_to_list",
            "constrained_vector_to_list",
            "stz_basis",
            "to_constrained",
            "to_unconstrained",
            "calc_log_jacobian",
            "generate_random_init",
            "NUTS_method",
            "IMH_method",
            "create_NUTS_core",
            "lpdf",
            "math",
            "log_sum_exp",
            "inv_logit",
            "logit",
            "distance",
            "squared_distance",
            "softmax"
          ),
          future.packages = c("RTMB","BayesRTMB")
          )
        })
        future::plan(future::sequential)
      } else {
        cat(paste0("直列サンプリングを開始します (chains = ", chains, ")...\n"))
        results_list <- lapply(1:chains, function(c) {
          run_chain(c, p_callback = NULL)
        })
      }

      mcmc_index <- seq(from = (warmup+1), to = (warmup+sampling), by = thin)
      accept_mat <- array(NA, dim=c(length(mcmc_index), chains))
      td_mat <- array(NA, dim=c(length(mcmc_index), chains))
      eps_vec <- numeric(chains)

      fit <- array(NA, dim=c(length(mcmc_index), chains, P_fixed + 1))
      dimnames(fit) <- list(
        iteration = NULL,
        chain = paste0("chain", 1:chains),
        variable = c("lp", pl_fixed$names)
      )

      if (P_random > 0) {
        random_fit <- array(NA, dim=c(length(mcmc_index), chains, P_random))
        dimnames(random_fit) <- list(
          iteration = NULL,
          chain = paste0("chain", 1:chains),
          variable = pl_random$names
        )
      } else {
        random_fit <- NULL
      }

      for(c in 1:chains){
        res <- results_list[[c]]
        fit[, c, 1] <- res$lp[mcmc_index]
        for(j in 1:P_fixed) fit[, c, j+1] <- res$para[mcmc_index, fixed_idx[j]]
        if (P_random > 0) {
          for (j in 1:P_random)
            random_fit[, c, j] <- res$para[mcmc_index, random_idx[j]]
        }
        accept_mat[,c] <- res$accept[mcmc_index]
        td_mat[,c]     <- res$treedepth[mcmc_index]
        eps_vec[c]     <- res$eps
      }

      eps_chains <- eps_vec
      accept_chains <- apply(accept_mat, 2, mean)
      treedepth_chains <- apply(td_mat, 2, max)
      names(eps_chains) <-
        names(accept_chains) <-
        names(treedepth_chains) <- paste0("chain", 1:chains)

      posterior_mean <- numeric(length(self$pl_full$names))
      names(posterior_mean) <- self$pl_full$names
      fixed_mean <- apply(fit[, , -1, drop = FALSE], 3, mean)
      posterior_mean[names(fixed_mean)] <- fixed_mean
      if (!is.null(random_fit)) {
        random_mean <- apply(random_fit, 3, mean)
        posterior_mean[names(random_mean)] <- random_mean
      }

      # --- MCMC_Fit インスタンスを返す ---
      res_obj <- MCMC_Fit$new(
        model          = self,
        fit            = fit,
        random_fit     = random_fit,
        eps            = eps_chains,
        accept         = accept_chains,
        treedepth      = treedepth_chains,
        laplace        = laplace,
        posterior_mean = posterior_mean
      )

      has_tran <- !is.null(self$transform)
      has_generate <- !is.null(self$generate)
      has_cf_corr <- any(sapply(self$par_list, function(x) x$type == "CF_corr"))

      if (has_tran || has_cf_corr) {
        res_obj$transformed_draws(self$transform)
      }

      if (has_generate) {
        res_obj$generated_quantities(self$generate)
      }

      return(res_obj)
    },

    #' @description Run Automatic Differentiation Variational Inference (ADVI).
    #' @param iter Integer; maximum number of iterations for the SGD optimization. Default is 10000.
    #' @param min_iter Integer; minimum number of iterations to run before checking for convergence. Default is 1000.
    #' @param tol_rel_obj Numeric; relative tolerance for the ELBO change to determine convergence. Default is 0.001.
    #' @param window_size Integer; window size for median smoothing in the convergence check. Default is 100.
    #' @param num_samples Integer; number of posterior samples to generate from the fitted variational distribution. Default is 1000.
    #' @param chains Integer; number of output chains (primarily for formatting compatibility with MCMC outputs). Default is 1.
    #' @param alpha Numeric; learning rate for the Adam optimizer. Default is 0.01.
    #' @param laplace Logical; whether to use Laplace approximation to marginalize random effects. Default is TRUE.
    #' @param print_freq Integer; iterations interval for progress output. Set to 0 to disable. Default is 100.
    #' @param fullrank Logical; whether to use a full-rank approximation to capture posterior correlations. Default is FALSE.
    #' @param seed Integer; random seed for reproducibility.
    #' @param init Optional numeric vector or list for initial parameter values. Default is NULL.
    #' @return A fitted `VB_Fit` object containing posterior samples and diagnostic information.
    variational = function(iter = 10000, min_iter = 1000,
                           tol_rel_obj = 0.001, window_size = 100,
                           num_samples = 1000, chains = 1, alpha = 0.01,
                           laplace = TRUE, print_freq = 100,
                           fullrank = FALSE,
                           seed = sample.int(1e6, 1), init = NULL) {

      set.seed(seed)
      cat("Preparing ADVI...\n")

      if (!is.null(init)) {
        init_list <- constrained_vector_to_list(init, self$par_list)
        unc_init_list <- to_unconstrained(init_list, self$par_list)
        init_full <- unlist(unc_init_list, use.names = FALSE)
      } else {
        init_full <- generate_random_init(self$pl_full, self$par_list, range = 2)
      }

      # ADVIは無制約空間で分布を学習するためヤコビアンが必須
      ad_setup <- self$build_ad_obj(init = init_full, laplace = laplace, jacobian_target = "all")
      ad_obj <- ad_setup$ad_obj

      res <- ADVI_method(
        model = ad_obj,
        par_list = self$par_list,
        pl_full = self$pl_full,
        iter = iter,
        min_iter = min_iter,
        tol_rel_obj = tol_rel_obj,
        window_size = window_size,
        num_samples = num_samples,
        chains = chains,
        alpha = alpha,
        laplace = laplace,
        print_freq = print_freq
      )

      # MCMC_Fitインスタンス作成のためのダミー変数
      samples_per_chain <- dim(res$fit)[1]
      eps_chains <- rep(NA, chains)
      accept_chains <- rep(NA, chains)
      treedepth_chains <- rep(NA, chains)
      names(eps_chains) <- names(accept_chains) <- names(treedepth_chains) <- paste0("chain", 1:chains)

      posterior_mean <- numeric(length(self$pl_full$names))
      names(posterior_mean) <- self$pl_full$names
      fixed_mean <- apply(res$fit[, , -1, drop = FALSE], 3, mean)
      posterior_mean[names(fixed_mean)] <- fixed_mean

      if (!is.null(res$random_fit)) {
        random_mean <- apply(res$random_fit, 3, mean)
        posterior_mean[names(random_mean)] <- random_mean
      }

      # MCMC_Fit オブジェクトとして結果を返す
      # MCMC_Fit ではなく ADVI_Fit をインスタンス化する
      res_obj <- VB_Fit$new(
        model          = self,
        fit            = res$fit,
        random_fit     = res$random_fit,
        elbo_history   = res$elbo_history, # ADVI.Rの戻り値に追加が必要なら渡す
        laplace        = laplace,
        posterior_mean = posterior_mean
      )

      # 変換量・生成量関数の適用
      has_tran <- !is.null(self$transform)
      has_generate <- !is.null(self$generate)
      has_cf_corr <- any(sapply(self$par_list, function(x) x$type == "CF_corr"))

      if (has_tran || has_cf_corr) {
        res_obj$transformed_draws(self$transform)
      }

      if (has_generate) {
        res_obj$generated_quantities(self$generate)
      }

      return(res_obj)
    },


    # 5. モデルコードの表示メソッド
    #' @description Print model code or model structure.
    #' @return The object itself, invisibly.
    print_code = function() {
      param_lines <- c("par_list <- list(")
      names_pl <- names(self$par_list)

      for (i in seq_along(names_pl)) {
        name <- names_pl[i]
        p <- self$par_list[[name]]
        args <- c()
        if (length(p$dim) > 1) {
          args <- c(args, paste0("dim = c(", paste(p$dim, collapse = ", "), ")"))
        } else {
          args <- c(args, paste0("dim = ", p$dim))
        }

        default_type <-
          if (length(p$dim) == 1)
            "vector"
        else if (length(p$dim) == 2)
          "matrix"
        else
          "array"
        if (!is.null(p$type) &&
            p$type != default_type)
          args <- c(args, paste0("type = \"", p$type, "\""))
        if (!is.null(p$lower)) args <- c(args, paste0("lower = ", p$lower))
        if (!is.null(p$upper)) args <- c(args, paste0("upper = ", p$upper))
        if (isTRUE(p$random)) args <- c(args, "random = TRUE")

        line <- paste0("  ", name, " = Dim(", paste(args, collapse = ", "), ")")
        if (i < length(names_pl)) line <- paste0(line, ",")
        param_lines <- c(param_lines, line)
      }
      param_lines <- c(param_lines, ")", "")

      lp_code <- deparse(self$log_prob, control = "useSource")
      lp_code[1] <- paste("log_prob <-", lp_code[1])

      full_code <- c(param_lines, lp_code)
      cat(paste(full_code, collapse = "\n"), "\n")
      invisible(full_code)
    }
  )
)
