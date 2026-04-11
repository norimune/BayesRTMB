#' RTMB model object
#'
#' @description
#' An R6 class representing a Bayesian model built with RTMB.
#' This class stores model components and provides methods for
#' building the automatic differentiation object, optimizing the
#' posterior, and drawing posterior samples.
#'
#' @param data 観測データを含む名前付きリスト。
#' @param par_list 評価済みのパラメータ定義リスト。
#' @param log_prob 実行可能な尤度・事前分布の計算関数（エラーハンドリング済み）。
#' @param transform 変換パラメータの計算関数（省略可能）。
#' @param generate 生成量の計算関数（省略可能）。
#' @param par_names パラメータの次元に対応する変数名のリスト。
#' @param init 初期値のリストまたは数値ベクトル。
#' @param view summaryで優先表示する変数名の文字ベクトル。
#' @param code `rtmb_code()` で生成された、ユーザー記述の元のモデル構文木 (AST)。
#'
#' @field data 観測データを含む名前付きリスト。
#' @field par_list 評価済みのパラメータ定義リスト。
#' @field log_prob 実行可能な尤度・事前分布の計算関数。
#' @field transform 変換パラメータの計算関数。
#' @field generate 生成量の計算関数。
#' @field par_names パラメータの次元に対応する変数名のリスト。
#' @field pl_full 内部で使用される、制約や次元情報を展開した完全なパラメータリスト。
#' @field formula モデル生成に使用されたフォーミュラ（ラッパー関数経由の場合）。
#' @field raw_data モデル生成に使用された元のデータフレーム（ラッパー関数経由の場合）。
#' @field family 分布族の文字列（ラッパー関数経由の場合）。
#' @field init 初期値のリストまたは数値ベクトル。
#' @field view summaryで優先表示する変数名の文字ベクトル。
#' @field code `rtmb_code()` で生成された、ユーザー記述の元のモデル構文木 (AST)。
#'
#' @import RTMB
RTMB_Model <- R6::R6Class(
  classname = "RTMB_Model",

  public = list(
    data       = NULL,
    par_list   = NULL,
    log_prob   = NULL,
    transform  = NULL,
    generate   = NULL,
    par_names  = NULL,
    pl_full    = NULL,
    formula    = NULL,
    raw_data   = NULL,
    family     = NULL,
    init       = NULL,
    view       = NULL,
    code       = NULL,

    # 1. コンストラクタ
    #' @description Create a new `RTMB_Model` object.
    initialize = function(data, par_list, log_prob,
                          transform = NULL, generate = NULL, par_names = NULL,
                          init = NULL, view = NULL, code = NULL) {
      self$data <- data
      self$par_names <- par_names
      self$init <- init
      self$view <- view
      self$code <- code
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

      self$pl_full <- parse_parameters(self$par_list, self$par_names)
      names(self$pl_full$init) <- self$pl_full$names

      # ここを prepare_init を経由するように変更
      init_vec <- self$prepare_init(self$init)
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

    # 初期値の整形と補完を行うヘルパーメソッド (引数がない場合はself$initを使うように変更)
    #' @description Prepare and format initial values for the model parameters.
    #' @param init_arg Optional list or numeric vector of initial values. If NULL, defaults to `self$init` or random generation.
    #' @return A flat numeric vector of constrained initial values.
    prepare_init = function(init_arg) {
      target_init <- if (!is.null(init_arg)) init_arg else self$init

      if (is.null(target_init)) {
        return(generate_random_init(self$pl_full, self$par_list, range = 2))
      }

      # 名前付きリストが渡された場合（部分的な初期値指定）
      if (is.list(target_init)) {
        # 1. まず全体をランダムな値で初期化し、リスト形式に変換
        base_vec <- generate_random_init(self$pl_full, self$par_list, range = 2)
        base_list <- constrained_vector_to_list(base_vec, self$par_list)

        # 2. ユーザーが指定したパラメータだけを上書き
        for (name in names(target_init)) {
          if (name %in% names(base_list)) {
            base_list[[name]] <- as.numeric(target_init[[name]])
          } else {
            warning(sprintf("初期値として指定された '%s' はモデルに存在しません。", name), call. = FALSE)
          }
        }
        # 3. 再び完全なフラットベクトルに戻して返す
        return(unlist(base_list, use.names = FALSE))
      }

      # ベクトルが渡された場合（全体指定）
      if (is.numeric(target_init)) {
        return(target_init)
      }

      stop("init は名前付きリスト、または数値ベクトルで指定してください。")
    },

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

      current_init <- self$prepare_init(init)

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
    #' @param num_estimate Integer; Number of estimate
    #' @param laplace Logical; whether to use Laplace approximation. Default is TRUE.
    #' @param init Optional initial values for parameters.
    #' @param control A list of control settings passed to the optimizer.
    #' @param optimizer Character; The optimizer to use, either "optim" or "nlminb". Default is "optim".
    #' @param method Character; The method for "optim" (e.g. "BFGS", "L-BFGS-B"). Default is "BFGS".
    #' @return A fitted `MAP_Fit` object.
    optimize = function(laplace = TRUE, init = NULL, num_estimate = 1, control = list(),
                        optimizer = "nlminb", method = "BFGS") {
      cat("Starting optimization...\n")

      opt_results <- list()
      obj_vals <- numeric(num_estimate)
      conv_codes <- numeric(num_estimate)

      # MAP推定ではヤコビアンは不要
      jac_target <- if (laplace) "random" else "none"

      # --- 修正: ループ外でADオブジェクトを1度だけ構築 ---
      ad_setup <- self$build_ad_obj(init = init, laplace = laplace, jacobian_target = jac_target)
      base_ad_obj <- ad_setup$ad_obj

      for (i in 1:num_estimate) {
        if (num_estimate > 1)
          cat(sprintf("Optimization run %d/%d...\r", i, num_estimate))

        res <- tryCatch({
          # --- 修正: 2回目以降、かつユーザー指定の初期値がない場合はランダムに更新 ---
          if (i > 1 && is.null(init)) {
            current_init <- self$prepare_init(NULL) # ランダム生成
            init_unc_list <- to_unconstrained(constrained_vector_to_list(current_init, self$par_list), self$par_list)
            base_ad_obj$par <- unlist(init_unc_list, use.names = FALSE)
          }

          if (optimizer == "nlminb") {
            if (is.null(control$iter.max)) control$iter.max <- 5000
            if (is.null(control$eval.max)) control$eval.max <- 5000
            if (is.null(control$rel.tol)) control$rel.tol <- 1e-8

            opt <- nlminb(
              start     = base_ad_obj$par,
              objective = base_ad_obj$fn,
              gradient  = base_ad_obj$gr,
              control   = control
            )
          } else if (optimizer == "optim") {
            if (is.null(control$maxit)) control$maxit <- 5000

            opt <- optim(
              par     = base_ad_obj$par,
              fn      = base_ad_obj$fn,
              gr      = base_ad_obj$gr,
              method  = method,
              control = control
            )
            opt$objective <- opt$value
          } else {
            stop("optimizerは 'optim' または 'nlminb' を指定してください。")
          }

          list(opt = opt, ad_obj = base_ad_obj)
        }, error = function(e) {
          warning("Optimization error: ", e$message, call. = FALSE)
          NULL
        })

        if (!is.null(res)) {
          opt_results[[i]] <- res
          obj_vals[i]      <- res$opt$objective
          conv_codes[i]    <- res$opt$convergence
        } else {
          opt_results[[i]] <- NULL
          obj_vals[i]      <- NA
          conv_codes[i]    <- NA
        }
      }
      cat("\n")

      # 有効な結果の中から Objective が最小のものを選択
      valid_idx <- which(!is.na(obj_vals))
      if (length(valid_idx) == 0) {
        stop("すべての最適化試行が失敗しました。初期値やモデルの定義を確認してください。")
      }

      best_idx <- valid_idx[which.min(obj_vals[valid_idx])]
      best_res <- opt_results[[best_idx]]

      # 推定結果の診断表示
      if (num_estimate == 1) {
        # 1回だけのときは、結果がどうなったかだけを短く報告
        status <- if (conv_codes[1] == 0) "converged" else "Not Converged"
        cat(sprintf("Optimization %s. Final objective: %.2f\n", status, obj_vals[1]))
      } else {
        # 複数回のときは、診断テーブルを表示
        cat("\nOptimization Diagnostics per estimate:\n")
        for (i in 1:num_estimate) {
          if (is.na(obj_vals[i])) {
            status_str <- "Failed"
          } else {
            status_str <- if (conv_codes[i] == 0) "Converged" else "Not Converged"
          }
          best_marker <- if (i == best_idx) "  <-- BEST" else ""
          cat(sprintf("  est%d: Objective = %10.2f, Code = %s (%s)%s\n",
                      i, obj_vals[i], as.character(conv_codes[i]), status_str, best_marker))
        }
        cat("\n")
      }

      ad_obj <- best_res$ad_obj
      opt    <- best_res$opt
      ad_obj$fn(opt$par)

      # --- 以降のコード（sdreportの計算など）は元のまま変更不要です ---
      sd_rep <- tryCatch(RTMB::sdreport(ad_obj), error = function(e) NULL)

      ad_obj <- best_res$ad_obj
      opt    <- best_res$opt
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
            names_def <- self$par_names[[name]]
            f_names <- generate_flat_names(name, p_info$dim, names_def)

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

      # 1. 事前に transform と generate の結果をリストとして計算
      tran_list <- NULL
      if (!is.null(self$transform)) {
        tran_list <- tryCatch(self$transform(self$data, con_est_list), error = function(e) NULL)
      }

      gq_list <- NULL
      if (!is.null(self$generate)) {
        tmp_con_list <- con_est_list
        if (!is.null(tran_list)) tmp_con_list <- c(tmp_con_list, tran_list)
        gq_list <- tryCatch(self$generate(self$data, tmp_con_list), error = function(e) NULL)
      }

      # 2. base_out を受け取るように変更した関数
      build_derived_df <- function(func, base_out, is_generate = FALSE) {
        if (is.null(func) || is.null(base_out) || length(base_out) == 0) return(NULL)

        u_base <- unc_est_vec
        L_u <- length(u_base)

        calc_derived <- function(u) {
          tmp_unc_list <- unconstrained_vector_to_list(u, self$par_list)
          tmp_con_list <- to_constrained(tmp_unc_list, self$par_list)

          if (is_generate && !is.null(self$transform)) {
            user_tran <- tryCatch(self$transform(self$data, tmp_con_list), error = function(e) NULL)
            if (!is.null(user_tran)) tmp_con_list <- c(tmp_con_list, user_tran)
          }

          res <- tryCatch(func(self$data, tmp_con_list), error = function(e) NULL)
          return(res)
        }

        flat_base <- unlist(base_out, use.names = FALSE)
        L_out <- length(flat_base)

        eps_diff <- 1e-5
        J <- matrix(0, nrow = L_out, ncol = L_u)
        for (i in 1:L_u) {
          u_tmp <- u_base
          u_tmp[i] <- u_tmp[i] + eps_diff
          tmp_out <- calc_derived(u_tmp)
          flat_tmp <- unlist(tmp_out, use.names = FALSE)
          J[, i] <- (flat_tmp - flat_base) / eps_diff
        }

        Cov_u <- diag(unc_se_vec^2, nrow = L_u, ncol = L_u)
        Cov_u[is.na(Cov_u)] <- 0

        if (!is.null(sd_rep) && !is.null(sd_rep$cov.fixed)) {
          idx_ran <- ad_obj$env$random
          idx_fix <- if (length(idx_ran) > 0) (1:L_u)[-idx_ran] else 1:L_u
          Cov_u[idx_fix, idx_fix] <- sd_rep$cov.fixed
        }

        se_out <- numeric(L_out)
        for (j in 1:L_out) {
          se_out[j] <- sqrt(sum((J[j, ] %*% Cov_u) * J[j, ]))
        }

        z_95 <- qnorm(0.975)
        low_out <- flat_base - z_95 * se_out
        up_out <- flat_base + z_95 * se_out

        names_vec <- c()
        for (name in names(base_out)) {
          val <- base_out[[name]]
          dim_val <- dim(val)
          if (is.null(dim_val)) dim_val <- length(val)
          names_def <- self$par_names[[name]]
          names_vec <- c(names_vec, generate_flat_names(name, dim_val, names_def))
        }

        df <- data.frame(
          Estimate     = flat_base,
          `Std. Error` = se_out,
          `Lower 95%`  = low_out,
          `Upper 95%`  = up_out,
          row.names    = names_vec,
          check.names  = FALSE
        )
        return(df)
      }

      # 3. 呼び出し時に第2引数として tran_list を渡す
      df_transform <- build_derived_df(self$transform, tran_list, is_generate = FALSE)

      if (!is.null(gq_list) && length(gq_list) > 0) {
        flat_base <- unlist(gq_list, use.names = FALSE)
        names_vec <- c()
        for (name in names(gq_list)) {
          val <- gq_list[[name]]
          dim_val <- dim(val)
          if (is.null(dim_val)) dim_val <- length(val)
          names_vec <- c(names_vec, generate_flat_names(name, dim_val, self$par_names[[name]]))
        }
        df_generate <- data.frame(
          Estimate     = flat_base,
          `Std. Error` = NA,
          `Lower 95%`  = NA,
          `Upper 95%`  = NA,
          row.names    = names_vec,
          check.names  = FALSE
        )
      } else {
        df_generate <- NULL
      }

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

      opt_history <- data.frame(estimate = 1:num_estimate, objective = obj_vals, code = conv_codes)

      # 4. 最後にMAP_Fitへ渡す
      res_obj <- MAP_Fit$new(
        model          = self,
        par_vec        = con_est_vec,
        par            = con_est_list,
        objective      = opt$objective,
        log_ml = log_ml,
        convergence = opt$convergence,
        sd_rep = sd_rep,
        df_fixed = df_fixed,
        random_effects = df_random,
        df_transform = df_transform,
        df_generate = df_generate,
        opt_history = opt_history,
        transform = tran_list,
        generate = gq_list
      )

      return(res_obj)
    },


    # 4. MCMC NUTS サンプリング メソッド
    #' @description Draw posterior samples from the model.
    #' @param sampling Number of sampling iterations. Default is 1000.
    #' @param warmup Number of warmup iterations. Default is 1000.
    #' @param chains Number of MCMC chains. Default is 4.
    #' @param thin Thinning interval. Default is 1.
    #' @param seed Random seed.
    #' @param delta Target acceptance rate for HMC/NUTS. Default is 0.8.
    #' @param max_treedepth Maximum tree depth for HMC/NUTS. Default is 10.
    #' @param parallel Logical; whether to run chains in parallel. Default is TRUE.
    #' @param laplace Logical; whether to use Laplace approximation. Default is FALSE.
    #' @param init Optional initial values for parameters.
    #' @param init_jitter sd of randomize initial values for parameters.
    #' @param save_csv Optional list for saving MCMC results. e.g., list(name = "model", dir = "BayesRTMB_mcmc").
    #' @return A fitted `MCMC_Fit` object.
    sample = function(sampling=1000, warmup=1000, chains=4,
                      thin=1, seed=sample.int(1e6,1),
                      delta=0.8, max_treedepth = 10,
                      parallel = FALSE, laplace = FALSE,
                      init = NULL, init_jitter = 0.1, save_csv = NULL) {

      if (!is.null(init)) init <- as.numeric(init)
      set.seed(seed)
      orig_pl <- self$par_list

      # --- CSV保存用の情報整理とディレクトリ作成 ---
      if (!is.null(save_csv)) {
        if (!is.list(save_csv)) stop("save_csv は list(name='...', dir='...') の形式で指定してください。")
        save_name <- if (!is.null(save_csv$name)) save_csv$name else "model"
        save_dir <- if (!is.null(save_csv$dir)) save_csv$dir else "BayesRTMB_mcmc"

        # 追加: 保存頻度（0ならチェイン終了時にまとめて1回で保存）
        save_freq <- if (!is.null(save_csv$freq)) save_csv$freq else 0

        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
        save_info <- list(name = save_name, dir = save_dir, freq = save_freq)
      } else {
        save_info <- NULL
      }

      random_flags <- sapply(orig_pl, function(x) isTRUE(x$random))

      if (laplace && any(random_flags)) {
        pl_fixed  <- parse_parameters(orig_pl[!random_flags], self$par_names)
        pl_random <- parse_parameters(orig_pl[random_flags], self$par_names)

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

        base_init <- self$prepare_init(init)

        unc_init_list <- to_unconstrained(constrained_vector_to_list(base_init, self$par_list), self$par_list)
        unc_init_vec <- unlist(unc_init_list, use.names = FALSE)

        if (init_jitter > 0) {
          unc_init_vec <- unc_init_vec + rnorm(length(unc_init_vec), mean = 0, sd = init_jitter)
        }
        unc_init_list_new <- unconstrained_vector_to_list(unc_init_vec, self$par_list)
        init_full_list <- to_constrained(unc_init_list_new, self$par_list)
        init_full <- unlist(init_full_list, use.names = FALSE)

        ad_setup <- self$build_ad_obj(init = init_full, laplace = laplace, jacobian_target = "all")
        ad_obj <- ad_setup$ad_obj

        # NUTS を実行 (save_info を渡す)
        res <- NUTS_method(
          model = ad_obj,
          sampling = sampling,
          warmup = warmup,
          delta = delta,
          max_treedepth = max_treedepth,
          chain = c,
          update_progress = p_callback,
          laplace = laplace,
          save_info = save_info
        )

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
            run_chain(c, p_callback = function(msg = "") p(message = msg))
          }, future.seed = TRUE,
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
      names(eps_chains) <- names(accept_chains) <- names(treedepth_chains) <- paste0("chain", 1:chains)

      posterior_mean <- numeric(length(self$pl_full$names))
      names(posterior_mean) <- self$pl_full$names
      fixed_mean <- apply(fit[, , -1, drop = FALSE], 3, mean)
      posterior_mean[names(fixed_mean)] <- fixed_mean
      if (!is.null(random_fit)) {
        random_mean <- apply(random_fit, 3, mean)
        posterior_mean[names(random_mean)] <- random_mean
      }

      if (!is.null(save_info)) {
        for (c in 1:chains) {
          backup_file <- file.path(save_info$dir, paste0(save_info$name, "-", c, ".csv"))

          # サンプラーの指標をまとめたデータフレームを作成
          df_metrics <- data.frame(
            iteration = mcmc_index,
            accept    = accept_mat[, c],
            treedepth = td_mat[, c],
            eps       = eps_vec[c]
          )

          if (!is.null(random_fit)) {
            df_out <- cbind(df_metrics,
                            as.data.frame(fit[, c, ]),
                            as.data.frame(random_fit[, c, ]))
          } else {
            df_out <- cbind(df_metrics,
                            as.data.frame(fit[, c, ]))
          }

          write.csv(df_out, file = backup_file, row.names = FALSE)
        }
      }

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
      #has_cf_corr <- any(sapply(self$par_list, function(x) x$type == "CF_corr"))

      if (has_tran) res_obj$transformed_draws(self$transform)
      if (has_generate) res_obj$generated_quantities(self$generate)

      return(res_obj)
    },

    #' @description Run Automatic Differentiation Variational Inference (ADVI).
    #' @param iter Integer; fixed number of iterations for the optimization. Default is 10000.
    #' @param tol_rel_obj Numeric; relative tolerance for the ELBO change to determine convergence. Default is 0.001.
    #' @param window_size Integer; window size for median smoothing in the convergence check. Default is 100.
    #' @param num_samples Integer; number of posterior samples to generate from the fitted variational distribution. Default is 1000.
    #' @param num_estimate Integer; number of times to run the VB estimation (treated as chains). Default is 4.
    #' @param alpha Numeric; learning rate for the Adam optimizer. Default is 0.01.
    #' @param laplace Logical; whether to use Laplace approximation to marginalize random effects. Default is TRUE.
    #' @param print_freq Integer; iterations interval for progress output. Set to 0 to disable. Default is 100.
    #' @param method Vector; method of Variational Inference Default is meanfield.
    #' @param parallel Logical; whether to run estimations in parallel. Default is FALSE.
    #' @param seed Integer; random seed for reproducibility.
    #' @param init Optional numeric vector or list for initial parameter values. Default is NULL.
    #' @param save_csv Optional list for saving VB results. e.g., list(name = "model", dir = "BayesRTMB_vb").
    #' @return A fitted `VB_Fit` object containing posterior samples and diagnostic information.
    variational = function(iter = 3000,
                           tol_rel_obj = 0.005, window_size = 100,
                           num_samples = 1000, num_estimate = 4, alpha = 0.01,
                           laplace = FALSE, print_freq = 1000,
                           method = c("meanfield", "fullrank", "hybrid"), parallel = FALSE,
                           seed = sample.int(1e6, 1), init = NULL, save_csv = NULL) {

      set.seed(seed)
      method <- match.arg(method)

      # --- CSV保存用の情報整理とディレクトリ作成 ---
      if (!is.null(save_csv)) {
        if (!is.list(save_csv)) stop("save_csv は list(name='...', dir='...') の形式で指定してください。")
        save_name <- if (!is.null(save_csv$name)) save_csv$name else "model_vb"
        save_dir <- if (!is.null(save_csv$dir)) save_csv$dir else "BayesRTMB_vb"
        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
        save_info <- list(name = save_name, dir = save_dir)
      } else {
        save_info <- NULL
      }

      run_advi <- function(c) {
        if (print_freq > 0) cat(sprintf("\n--- VB推定開始: est%d ---\n", c))

        ad_setup <- self$build_ad_obj(init = init, laplace = laplace, jacobian_target = "all")

        res <- ADVI_method(
          model = ad_setup$ad_obj, par_list = self$par_list, pl_full = self$pl_full,
          iter = iter, tol_rel_obj = tol_rel_obj,
          window_size = window_size, num_samples = num_samples, alpha = alpha,
          laplace = laplace, print_freq = print_freq, method = method
        )
        return(res)
      }

      results_list <- list()
      if (parallel && num_estimate > 1) {
        future::plan(future::multisession, workers = num_estimate)
        cat(paste0("並列VB推定を開始します (num_estimate = ", num_estimate, ")...\n"))

        if (requireNamespace("progressr", quietly = TRUE)) {
          progressr::handlers(global = TRUE)

          update_interval <- 100
          steps_per_chain <- ceiling(iter / update_interval)
          total_steps <- steps_per_chain * num_estimate

          results_list <- progressr::with_progress({
            p <- progressr::progressor(steps = total_steps)

            run_advi_prog <- function(c) {

              ad_setup <- self$build_ad_obj(init = init, laplace = laplace, jacobian_target = "all")

              update_prog_fn <- function(amount = 1) {
                p(amount = amount)
              }

              res <- ADVI_method(
                model = ad_setup$ad_obj, par_list = self$par_list, pl_full = self$pl_full,
                iter = iter, tol_rel_obj = tol_rel_obj,
                window_size = window_size, num_samples = num_samples, alpha = alpha,
                laplace = laplace, print_freq = 0, method = method,
                update_progress = update_prog_fn, update_interval = update_interval
              )
              return(res)
            }

            future.apply::future_lapply(1:num_estimate, function(c) {
              run_advi_prog(c)
            }, future.seed = TRUE, future.packages = "RTMB")
          })

        } else {
          cat("※プログレスバーを表示するには 'progressr' パッケージをインストールしてください。\n")
          results_list <- future.apply::future_lapply(1:num_estimate, function(c) {
            run_advi(c)
          }, future.seed = TRUE, future.packages = c("RTMB","BayesRTMB"))
        }

        future::plan(future::sequential)
      } else {
        cat(paste0("直列VB推定を開始します (num_estimate = ", num_estimate, ")...\n"))
        results_list <- lapply(1:num_estimate, function(c) {
          run_advi(c)
        })
      }

      # 推定結果をまとめる
      P_fixed <- dim(results_list[[1]]$fit)[3] - 1
      P_random <- if (!is.null(results_list[[1]]$random_fit)) dim(results_list[[1]]$random_fit)[3] else 0

      fit <- array(NA, dim = c(num_samples, num_estimate, P_fixed + 1))
      dimnames(fit) <- list(
        iteration = NULL,
        chain = paste0("est", 1:num_estimate),
        variable = dimnames(results_list[[1]]$fit)[[3]]
      )

      if (P_random > 0) {
        random_fit <- array(NA, dim = c(num_samples, num_estimate, P_random))
        dimnames(random_fit) <- list(
          iteration = NULL,
          chain = paste0("est", 1:num_estimate),
          variable = dimnames(results_list[[1]]$random_fit)[[3]]
        )
      } else {
        random_fit <- NULL
      }

      elbo_history_list <- list()
      elbo_final_vec <- numeric(num_estimate)
      rel_obj_vec <- numeric(num_estimate)

      for (c in 1:num_estimate) {
        res <- results_list[[c]]
        fit[, c, ] <- res$fit[, 1, ]
        if (P_random > 0) random_fit[, c, ] <- res$random_fit[, 1, ]
        elbo_history_list[[c]] <- res$elbo_history
        elbo_final_vec[c] <- res$elbo_final
        rel_obj_vec[c] <- res$rel_obj_final
      }

      # --- 推定完了後に各estimateの事後サンプルをCSVへ一括保存 ---
      if (!is.null(save_info)) {
        for (c in 1:num_estimate) {
          backup_file <- file.path(save_info$dir, paste0(save_info$name, "-", c, ".csv"))

          if (!is.null(random_fit)) {
            df_out <- cbind(iteration = 1:num_samples, as.data.frame(fit[, c, ]), as.data.frame(random_fit[, c, ]))
          } else {
            df_out <- cbind(iteration = 1:num_samples, as.data.frame(fit[, c, ]))
          }

          write.csv(df_out, file = backup_file, row.names = FALSE)
        }
      }

      best_chain <- which.max(elbo_final_vec)

      posterior_mean <- numeric(length(self$pl_full$names))
      names(posterior_mean) <- self$pl_full$names
      fixed_mean <- apply(fit[, best_chain, -1, drop = FALSE], 3, mean)
      posterior_mean[names(fixed_mean)] <- fixed_mean

      if (!is.null(random_fit)) {
        random_mean <- apply(random_fit[, best_chain, , drop=FALSE], 3, mean)
        posterior_mean[names(random_mean)] <- random_mean
      }

      cat("\nConvergence Diagnostics per estimate:\n")
      for (c in 1:num_estimate) {
        status <- if (!is.na(rel_obj_vec[c]) && rel_obj_vec[c] < tol_rel_obj) "Converged" else "Not Converged"
        best_marker <- if (c == best_chain) "  <-- BEST" else ""

        cat(sprintf("  est%d: ELBO = %10.2f, Final rel_obj = %.5f (%s)%s\n",
                    c, elbo_final_vec[c], rel_obj_vec[c], status, best_marker))
      }
      cat("\n")

      best_mu_history <- results_list[[best_chain]]$mu_history

      res_obj <- VB_Fit$new(
        model          = self,
        fit            = fit,
        random_fit     = random_fit,
        elbo_history   = elbo_history_list,
        laplace        = laplace,
        posterior_mean = posterior_mean,
        ELBO           = elbo_final_vec,
        rel_obj_vals   = rel_obj_vec,
        best_chain     = best_chain,
        mu_history     = best_mu_history
      )

      has_tran <- !is.null(self$transform)
      has_generate <- !is.null(self$generate)
      #has_cf_corr <- any(sapply(self$par_list, function(x) x$type == "CF_corr"))

      if (has_tran){
        cat("Calculating transformed parameters...\n")
        res_obj$transformed_draws(self$transform)
      }
      if (has_generate) {
        cat("Calculating generated quantities...\n")
        res_obj$generated_quantities(self$generate)
      }


      return(res_obj)
    },


    # 5. モデルコードの表示メソッド
    #' @description Print model code or model structure.
    #' @return The object itself, invisibly.
    print_code = function() {
      if (is.null(self$code)) {
        cat("Model code is not available (Not built with rtmb_code).\n")
        return(invisible(self))
      }

      cat("=== RTMB Model Code ===\n\n")
      cat("rtmb_code(\n")

      blocks <- names(self$code)
      n_blocks <- length(blocks)

      for (i in seq_along(blocks)) {
        block <- blocks[i]
        expr <- self$code[[block]]

        # 修正箇所: width.cutoff = 500L を追加して勝手な改行を防ぐ
        lines <- deparse(expr, width.cutoff = 500L, control = "useSource")

        # インデントの調整
        lines <- paste0("  ", lines)

        # 最初の "{" を "ブロック名 = {" に書き換える
        if (trimws(lines[1]) == "{") {
          lines[1] <- paste0("  ", block, " = {")
        } else {
          # 中身が1行だけで {} で囲まれていなかった場合の安全対策
          lines <- c(paste0("  ", block, " = {"), paste0("  ", lines), "  }")
        }

        # 最後のブロック以外は、閉じカッコ "}" の後にカンマを追加する
        if (i < n_blocks) {
          lines[length(lines)] <- paste0(lines[length(lines)], ",")
        }

        cat(paste(lines, collapse = "\n"), "\n")
      }

      cat(")\n\n")
      invisible(self)
    }
  )
)
