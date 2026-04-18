#' Restore MCMC Fit from CSV
#'
#' @param model An RTMB_Model object.
#' @param name Base name of the saved CSVs.
#' @param dir Directory where CSVs are saved. Default is "BayesRTMB_mcmc".
#' @param chains Number of chains. Default is 4.
#' @param laplace Logical; whether Laplace approximation was used. Default is FALSE.
#' @return An MCMC_Fit object.
#' @export
read_mcmc_csv <- function(model, name, dir = "BayesRTMB_mcmc", chains = 4, laplace = FALSE) {

  test_file <- file.path(dir, paste0(name, "-1.csv"))
  if (!file.exists(test_file)) stop(paste("File not found:", test_file))

  # --- 修正: check.names = FALSE を追加 ---
  test_dat <- read.csv(test_file, header = TRUE, check.names = FALSE)
  n_samples <- nrow(test_dat)

  orig_pl <- model$par_list
  random_flags <- sapply(orig_pl, function(x) isTRUE(x$random))

  if (laplace && any(random_flags)) {
    pl_fixed  <- BayesRTMB:::parse_parameters(orig_pl[!random_flags], model$par_names)
    pl_random <- BayesRTMB:::parse_parameters(orig_pl[random_flags], model$par_names)
  } else {
    pl_fixed   <- model$pl_full
    pl_random  <- NULL
  }

  P_fixed <- length(pl_fixed$names)
  P_random <- if(!is.null(pl_random)) length(pl_random$names) else 0

  fit <- array(NA, dim=c(n_samples, chains, P_fixed + 1))
  dimnames(fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = c("lp", pl_fixed$names))

  if (P_random > 0) {
    random_fit <- array(NA, dim=c(n_samples, chains, P_random))
    dimnames(random_fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = pl_random$names)
  } else {
    random_fit <- NULL
  }

  accept_mat <- array(NA, dim=c(n_samples, chains))
  td_mat <- array(NA, dim=c(n_samples, chains))
  eps_vec <- numeric(chains)

  for (c in 1:chains) {
    file_path <- file.path(dir, paste0(name, "-", c, ".csv"))
    if (!file.exists(file_path)) stop(paste("File not found:", file_path))

    # --- 修正: check.names = FALSE を追加 ---
    dat <- read.csv(file_path, header = TRUE, check.names = FALSE)

    if (nrow(dat) != n_samples) {
      warning(sprintf("Chain %d has %d iterations, expected %d", c, nrow(dat), n_samples))
    }

    accept_mat[, c] <- dat$accept
    td_mat[, c] <- dat$treedepth
    eps_vec[c] <- dat$eps[1]

    fit[, c, "lp"] <- dat$lp
    for (p_name in pl_fixed$names) {
      fit[, c, p_name] <- dat[[p_name]]
    }

    if (P_random > 0) {
      for (p_name in pl_random$names) {
        random_fit[, c, p_name] <- dat[[p_name]]
      }
    }
  }

  eps_chains <- eps_vec
  accept_chains <- apply(accept_mat, 2, mean)
  treedepth_chains <- apply(td_mat, 2, max)
  names(eps_chains) <- names(accept_chains) <- names(treedepth_chains) <- paste0("chain", 1:chains)

  posterior_mean <- numeric(length(model$pl_full$names))
  names(posterior_mean) <- model$pl_full$names
  fixed_mean <- apply(fit[, , -1, drop = FALSE], 3, mean)
  posterior_mean[names(fixed_mean)] <- fixed_mean

  if (!is.null(random_fit)) {
    random_mean <- apply(random_fit, 3, mean)
    posterior_mean[names(random_mean)] <- random_mean
  }

  res_obj <- MCMC_Fit$new(
    model          = model,
    fit            = fit,
    random_fit     = random_fit,
    eps            = eps_chains,
    accept         = accept_chains,
    treedepth      = treedepth_chains,
    laplace        = laplace,
    posterior_mean = posterior_mean
  )

  has_tran <- !is.null(model$transform)
  has_generate <- !is.null(model$generate)
  has_cf_corr <- any(sapply(model$par_list, function(x) x$type == "CF_corr"))

  if (has_tran || has_cf_corr) res_obj$transformed_draws(model$transform)
  if (has_generate) res_obj$generated_quantities(model$generate)

  return(res_obj)
}
# summaryの結果にクラスを付与するヘルパー
as_summary_df <- function(df) {
  class(df) <- c("summary_BayesRTMB", "data.frame")
  return(df)
}

#' print for summary_BayesRTMB class
#' @param digits integer
#' @export
print.summary_BayesRTMB <- function(x, digits = NULL,...) {
  df <- x

  if (is.null(digits)) {
    digits <- attr(x, "digits")
    if (is.null(digits)) digits <- 2
  }

  # 1. 各列を文字列に変換 (表示形式の指定)
  out_char <- as.data.frame(lapply(names(df), function(cn) {
    val <- df[[cn]]
    if (is.numeric(val)) {
      if (grepl("ess|iter|count", cn, ignore.case = TRUE)) {
        # ESSなどは整数 (NA対策込み)
        ifelse(is.na(val), "NA", sprintf("%.0f", val))
      } else {
        fmt <- paste0("%.", digits, "f")
        ifelse(is.na(val), "NA", sprintf(fmt, val))
      }
    } else {
      ifelse(is.na(val), "NA", as.character(val))
    }
  }), stringsAsFactors = FALSE)
  colnames(out_char) <- names(df)

  col_widths <- sapply(seq_along(names(out_char)), function(i) {
    max(nchar(names(out_char)[i]), nchar(out_char[, i]), na.rm = TRUE)
  })

  # 2. 列幅の計算
  col_widths <- sapply(seq_along(names(out_char)), function(i) {
    max(nchar(names(out_char)[i]), nchar(out_char[, i]))
  })

  # 3. ヘッダーの表示 (タイトルはすべて右揃え)
  header_parts <- sapply(seq_along(col_widths), function(i) {
    sprintf(paste0("%", col_widths[i], "s"), names(out_char)[i])
  })
  cat(paste(header_parts, collapse = "  "), "\n")

  # 4. データの表示 (1列目のデータのみ左揃え)
  for (r in seq_len(nrow(out_char))) {
    row_parts <- sapply(seq_along(col_widths), function(c) {
      if (c == 1) {
        # 変数名データのみ左揃え
        sprintf(paste0("%-", col_widths[c], "s"), out_char[r, c])
      } else {
        # 数値データは右揃え
        sprintf(paste0("%", col_widths[c], "s"), out_char[r, c])
      }
    })
    cat(paste(row_parts, collapse = "  "), "\n")
  }

  return(invisible(x))
}

#' 周辺効果 (Conditional Effects) を計算する
#' @export
conditional_effects <- function(fit, effect, ...) {
  UseMethod("conditional_effects")
}

#' ce_rtmb クラス専用のプロットメソッド (Base R)
#' @method conditional_effects mcmc_fit
#' @export
conditional_effects.mcmc_fit <- function(fit, effect, resolution = 100, prob = 0.95, ...) {
  model_obj <- fit$model
  if (is.null(model_obj$formula) || is.null(model_obj$raw_data)) {
    stop("このモデルオブジェクトには formula または元のデータが含まれていません。")
  }

  form <- model_obj$formula
  raw_data <- model_obj$raw_data
  fam <- model_obj$family

  # 交互作用かどうかの判定と分割
  eff_vars <- strsplit(effect, ":")[[1]]
  if (length(eff_vars) > 2) {
    stop("3つ以上の変数の交互作用プロットには現在対応していません。")
  }
  eff1 <- eff_vars[1]
  eff2 <- if (length(eff_vars) == 2) eff_vars[2] else NULL

  # 1. ベースとなるデータフレームの作成 (他の変数を平均/最頻値に固定)
  base_data <- lapply(raw_data, function(x) {
    if (is.numeric(x)) {
      mean(x, na.rm = TRUE)
    } else if (is.factor(x) || is.character(x)) {
      tbl <- table(x)
      factor(names(tbl)[which.max(tbl)], levels = levels(as.factor(x)))
    } else {
      x[1]
    }
  })
  base_data <- as.data.frame(base_data)

  # 2. effect変数のみを動かした newdata を作成
  val1 <- raw_data[[eff1]]
  if (is.null(val1)) stop(sprintf("データの中に変数 '%s' が見つかりません。", eff1))

  is_numeric1 <- is.numeric(val1)
  if (is_numeric1) {
    seq1 <- seq(min(val1, na.rm = TRUE), max(val1, na.rm = TRUE), length.out = resolution)
  } else {
    seq1 <- sort(unique(val1))
  }

  if (!is.null(eff2)) {
    val2 <- raw_data[[eff2]]
    if (is.null(val2)) stop(sprintf("データの中に変数 '%s' が見つかりません。", eff2))

    # 第2変数が連続値で種類が多い場合は、代表的な3点(Mean-1SD, Mean, Mean+1SD)に絞る
    if (is.numeric(val2) && length(unique(val2)) > 5) {
      seq2 <- c(mean(val2, na.rm=TRUE) - sd(val2, na.rm=TRUE),
                mean(val2, na.rm=TRUE),
                mean(val2, na.rm=TRUE) + sd(val2, na.rm=TRUE))
      seq2 <- round(seq2, 2) # 表示が見やすいように丸める
    } else {
      seq2 <- sort(unique(val2))
    }

    grid_data <- expand.grid(eff1 = seq1, eff2 = seq2)
    names(grid_data) <- c(eff1, eff2)

    newdata <- base_data[rep(1, nrow(grid_data)), , drop = FALSE]
    newdata[[eff1]] <- grid_data[[eff1]]
    newdata[[eff2]] <- grid_data[[eff2]]
  } else {
    newdata <- base_data[rep(1, length(seq1)), , drop = FALSE]
    newdata[[eff1]] <- seq1
  }

  # 3. デザイン行列の作成
  rhs <- delete.response(terms(form))
  X_new <- model.matrix(rhs, data = newdata)

  # 4. 事後サンプルの取得 (固定効果のみ)
  beta_samples <- fit$draws(pars = "beta", inc_random = FALSE, inc_transform = FALSE, inc_generate = FALSE)
  I <- dim(beta_samples)[1]
  C <- dim(beta_samples)[2]
  P <- dim(beta_samples)[3]
  beta_flat <- matrix(beta_samples, nrow = I * C, ncol = P)

  # 5. 線形予測子の計算
  eta <- X_new %*% t(beta_flat)

  # 6. 逆リンク関数で期待値に変換
  if (is.null(fam)) fam <- "gaussian"
  inv_link <- switch(fam,
                     "gaussian" = , "lognormal" = , "student_t" = function(x) x,
                     "poisson" = , "neg_binomial" = , "gamma" = exp,
                     "bernoulli" = , "binomial" = plogis,
                     function(x) x
  )
  mu <- inv_link(eta)

  # 7. 事後分布の要約
  alpha <- 1 - prob
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2

  res_df <- data.frame(
    estimate   = apply(mu, 1, mean),
    lower      = apply(mu, 1, quantile, probs = lower_q),
    upper      = apply(mu, 1, quantile, probs = upper_q)
  )
  res_df <- cbind(newdata[, eff_vars, drop = FALSE], res_df)

  # 8. 結果をリストにまとめてクラスを付与
  res <- list(data = res_df, effect_vars = eff_vars, is_numeric = is_numeric1)
  class(res) <- "ce_rtmb"
  return(res)
}

#' ce_rtmb クラス専用のプロットメソッド (Base R)
#' @method plot ce_rtmb
#' @export
plot.ce_rtmb <- function(x, ...) {
  df <- x$data
  eff_vars <- x$effect_vars
  eff1 <- eff_vars[1]
  has_interaction <- length(eff_vars) > 1

  x_val <- df[[eff1]]
  y_est <- df$estimate
  y_low <- df$lower
  y_up  <- df$upper

  if (!has_interaction) {
    # --- 交互作用なしの場合 ---
    col_line <- rgb(0, 0.45, 0.7)
    col_ribbon <- rgb(0, 0.45, 0.7, 0.2)

    if (x$is_numeric) {
      plot(x_val, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlab = eff1, ylab = "Predicted value", main = paste("Conditional effect of", eff1), ...)
      polygon(c(x_val, rev(x_val)), c(y_low, rev(y_up)), col = col_ribbon, border = NA)
      lines(x_val, y_est, col = col_line, lwd = 2)
    } else {
      x_num <- as.numeric(as.factor(x_val))
      plot(x_num, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlim = c(min(x_num) - 0.5, max(x_num) + 0.5), xaxt = "n",
           xlab = eff1, ylab = "Predicted value", main = paste("Conditional effect of", eff1), ...)
      axis(1, at = x_num, labels = as.character(x_val))
      segments(x0 = x_num, y0 = y_low, x1 = x_num, y1 = y_up, col = col_line, lwd = 2)
      points(x_num, y_est, col = col_line, pch = 16, cex = 1.5)
    }

  } else {
    # --- 交互作用ありの場合 ---
    eff2 <- eff_vars[2]
    groups <- unique(df[[eff2]])
    n_groups <- length(groups)

    # グループ数に応じたカラーパレットの作成
    cols_line <- hcl.colors(n_groups, palette = "Dark 2")
    cols_ribbon <- sapply(cols_line, function(col) {
      rgb_val <- col2rgb(col) / 255
      rgb(rgb_val[1], rgb_val[2], rgb_val[3], 0.2) # 透過度0.2
    })

    if (x$is_numeric) {
      plot(x_val, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlab = eff1, ylab = "Predicted value",
           main = paste("Conditional effect of", eff1, "by", eff2), ...)

      for (i in seq_along(groups)) {
        idx <- df[[eff2]] == groups[i]
        xv <- df[[eff1]][idx]
        polygon(c(xv, rev(xv)), c(df$lower[idx], rev(df$upper[idx])), col = cols_ribbon[i], border = NA)
        lines(xv, df$estimate[idx], col = cols_line[i], lwd = 2)
      }
    } else {
      x_fct <- as.factor(x_val)
      x_num <- as.numeric(x_fct)

      plot(x_num, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlim = c(min(x_num) - 0.5, max(x_num) + 0.5), xaxt = "n",
           xlab = eff1, ylab = "Predicted value",
           main = paste("Conditional effect of", eff1, "by", eff2), ...)
      axis(1, at = unique(x_num), labels = levels(x_fct))

      # エラーバーが重ならないようにX軸を少しずらす
      offset_step <- 0.1
      offsets <- seq(-offset_step * (n_groups-1)/2, offset_step * (n_groups-1)/2, length.out = n_groups)

      for (i in seq_along(groups)) {
        idx <- df[[eff2]] == groups[i]
        xv <- x_num[idx] + offsets[i]
        segments(x0 = xv, y0 = df$lower[idx], x1 = xv, y1 = df$upper[idx], col = cols_line[i], lwd = 2)
        points(xv, df$estimate[idx], col = cols_line[i], pch = 16, cex = 1.5)
      }
    }

    # 凡例の追加
    legend("topright", title = eff2, legend = format(groups, digits = 3),
           col = cols_line, lty = 1, pch = ifelse(x$is_numeric, NA, 16), lwd = 2, bty = "n")
  }
  invisible(x)
}

#' ce_rtmb クラスのprintメソッド (自動的にplotを呼ぶ)
#' @method print ce_rtmb
#' @export
print.ce_rtmb <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}
