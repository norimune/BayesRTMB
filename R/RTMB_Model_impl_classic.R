.validate_classic_model <- function(self, private) {
  if (is.null(self$extra) || !isTRUE(self$extra$source == "wrapper")) {
    stop(
      "classic() is available only for models created by wrapper functions.",
      call. = FALSE
    )
  }

  if (!identical(self$extra$prior_type, "flat")) {
    stop(
      "classic() is available only for models created with prior_flat(). ",
      "Use optimize() or sample() for models with informative priors.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.resolve_classic_settings <- function(self, private, df_method = "auto") {
  type <- self$type %||% "generic"
  has_random_parameters <- any(vapply(self$par_list, function(p) isTRUE(p$random), logical(1)))

  if (is.numeric(df_method)) {
    use_reml <- type %in% c("lm", "lmer", "ttest", "corr", "mediation")
    return(list(
      use_reml = use_reml,
      use_laplace = use_reml || has_random_parameters,
      df_method = "numeric",
      df_value = as.numeric(df_method),
      target_vars = if (use_reml && !is.null(self$extra$marginal)) self$extra$marginal else character(0),
      show_df = TRUE
    ))
  }

  if (is.null(df_method)) df_method <- "auto"
  df_method <- tolower(df_method)

  if (df_method %in% c("infinite", "none")) {
    df_method <- "inf"
  }

  valid <- c("auto", "residual", "satterthwaite", "inf")
  if (!df_method %in% valid) {
    stop(
      "Invalid df_method in classic(). ",
      "Use 'auto', 'residual', 'satterthwaite', or 'inf'.",
      call. = FALSE
    )
  }

  use_reml <- type %in% c("lm", "lmer", "ttest", "corr", "mediation")

  if (identical(df_method, "auto")) {
    df_method <- if (type %in% c("lmer")) {
      "satterthwaite"
    } else if (type == "ttest") {
      # Satterthwaite for independent unequal-variance t-tests
      if (isFALSE(self$extra$paired) && isFALSE(self$extra$var_equal)) {
        "satterthwaite"
      } else {
        "residual"
      }
    } else if (type %in% c("lm", "corr", "mediation")) {
      "residual"
    } else {
      "inf"
    }
  }

  target_vars <- if (use_reml && !is.null(self$extra$marginal)) {
    self$extra$marginal
  } else {
    character(0)
  }

  list(
    use_reml = use_reml,
    use_laplace = use_reml || has_random_parameters,
    df_method = df_method,
    df_value = NULL,
    target_vars = target_vars,
    show_df = !identical(df_method, "inf")
  )
}

.classic_pipeline <- function(self, private,
                             df_method = "auto",
                             se_method = c("wald", "sampling"),
                             num_samples = 1000,
                             seed = 123,
                             view = NULL,
                             map = NULL,
                             fixed = NULL,
                             ...) {
  se_method <- match.arg(se_method)

  if (!is.null(fixed)) {
    return(private$.dispatch_fixed(
      .method_to_call = "classic",
      df_method = df_method,
      se_method = se_method,
      num_samples = num_samples,
      seed = seed,
      view = view,
      map = map,
      fixed = fixed,
      ...
    ))
  }

  private$.validate_classic_model()

  if (isTRUE(self$type == "corr") &&
      (self$extra$corr_method %||% "reml") %in% c("pearson", "spearman")) {
    return(.classic_corr_test_fit(self, method = self$extra$corr_method, view = view))
  }

  settings <- private$.resolve_classic_settings(df_method)
  has_random_parameters <- any(vapply(self$par_list, function(p) isTRUE(p$random), logical(1)))

  # For classic mode, we use .fit_rtmb with apply_prior_correction = FALSE.
  # Constraint targets (fixed parameters for REML) are handled via target_vars.
  model_to_fit <- self

  raw <- model_to_fit$.__enclos_env__$private$.fit_rtmb(
    laplace = settings$use_laplace,
    init = NULL,
    num_estimate = 1,
    control = list(),
    optimizer = "nlminb",
    method = "BFGS",
    map = map,
    se = TRUE,
    reml = settings$use_reml,
    target_vars = settings$target_vars,
    jacobian_target = NULL,
    apply_prior_correction = FALSE
  )

  info_log_lik <- NULL
  restricted_log_lik <- NULL
  if (isTRUE(settings$use_reml)) {
    restricted_log_lik <- -raw$opt$objective
    raw_ic <- model_to_fit$.__enclos_env__$private$.fit_rtmb(
      laplace = has_random_parameters,
      init = raw$par,
      num_estimate = 1,
      control = list(),
      optimizer = "nlminb",
      method = "BFGS",
      map = map,
      se = FALSE,
      reml = FALSE,
      target_vars = character(0),
      jacobian_target = NULL,
      apply_prior_correction = FALSE,
      verbose = FALSE
    )
    info_log_lik <- -raw_ic$opt$objective
  }

  comp <- private$.build_classic_components(
    raw = raw,
    settings = settings,
    se_method = se_method,
    num_samples = num_samples,
    seed = seed,
    info_log_lik = info_log_lik
  )

  fit <- Classic_Fit$new(
    model = self,
    par_vec = unlist(raw$par, use.names = FALSE),
    par = raw$par,
    objective = raw$opt$objective,
    log_lik = comp$log_lik,
    restricted_log_lik = restricted_log_lik,
    convergence = raw$opt$convergence,
    sd_rep = raw$sd_rep,
    df_fixed = comp$df_fixed,
    random_effects = comp$random_effects,
    df_transform = comp$df_transform,
    df_generate = comp$df_generate,
    opt_history = raw$opt_history,
    transform = comp$transform,
    generate = comp$generate,
    se_samples = comp$se_samples,
    par_unc = raw$par_unc,
    vcov_unc = raw$vcov_unc,
    ci_method = se_method,
    laplace = settings$use_laplace,
    map = raw$map,
    test_results = comp$test_results,
    view = view,
    vcov = comp$vcov,
    df_method = settings$df_method,
    idx_fix_active = raw$idx_fix_active,
    show_df = settings$show_df,
    rss = comp$rss,
    df_residual = comp$df_residual
  )

  class(fit) <- unique(c(
    paste0("rtmb_", self$type),
    "Classic_Fit",
    class(fit)
  ))

  fit
}

.classic_corr_test_fit <- function(self, method = c("pearson", "spearman"), view = NULL) {
  method <- match.arg(method)
  Y <- as.matrix(self$data$Y)
  P_y <- self$data$P_y %||% ncol(Y)
  P_x <- self$data$P_x %||% 0
  Y_target <- Y[, seq_len(P_y), drop = FALSE]
  var_names <- colnames(Y_target)
  if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(Y_target)))

  if (ncol(Y_target) < 2) {
    stop("rtmb_corr() requires at least two numeric variables for correlation tests.", call. = FALSE)
  }

  pairs <- utils::combn(seq_len(ncol(Y_target)), 2)

  if (P_x > 0) {
    missing_arg <- self$extra$missing %||% "listwise"
    if (missing_arg == "pairwise") {
      n_complete <- nrow(Y)
      Y_corr <- if (method == "spearman") {
        apply(Y, 2, rank, na.last = "keep", ties.method = "average")
      } else {
        Y
      }
      R <- stats::cor(Y_corr, use = "pairwise.complete.obs")
    } else {
      ok <- stats::complete.cases(Y)
      Y_complete <- Y[ok, , drop = FALSE]
      n_complete <- nrow(Y_complete)
      if (n_complete <= P_x + 2) {
        stop("Not enough complete observations to compute partial correlations.", call. = FALSE)
      }
      Y_corr <- if (method == "spearman") {
        apply(Y_complete, 2, rank, ties.method = "average")
      } else {
        Y_complete
      }
      R <- stats::cor(Y_corr, use = "complete.obs")
    }
    R_yy <- R[seq_len(P_y), seq_len(P_y), drop = FALSE]
    x_idx <- P_y + seq_len(P_x)
    R_yx <- R[seq_len(P_y), x_idx, drop = FALSE]
    R_xx <- R[x_idx, x_idx, drop = FALSE]
    P_cov <- R_yy - R_yx %*% MASS::ginv(R_xx) %*% t(R_yx)
    D <- diag(1 / sqrt(diag(P_cov)), nrow = P_y, ncol = P_y)
    P_corr <- D %*% P_cov %*% D

    rows <- vector("list", ncol(pairs))
    pcorr_values <- numeric(ncol(pairs))
    row_names <- character(ncol(pairs))
    df_val <- n_complete - P_x - 2

    for (k in seq_len(ncol(pairs))) {
      i <- pairs[1, k]
      j <- pairs[2, k]
      est <- P_corr[i, j]
      est <- pmax(pmin(est, 0.999999), -0.999999)
      pcorr_values[k] <- est

      t_val <- est * sqrt(df_val / pmax(1 - est^2, 1e-12))
      p_val <- 2 * stats::pt(-abs(t_val), df = df_val)
      z_se <- if (n_complete - P_x - 3 > 0) 1 / sqrt(n_complete - P_x - 3) else NA_real_
      if (is.na(z_se)) {
        lower <- upper <- NA_real_
      } else {
        z_val <- atanh(est)
        lower <- tanh(z_val - stats::qnorm(0.975) * z_se)
        upper <- tanh(z_val + stats::qnorm(0.975) * z_se)
      }

      rows[[k]] <- data.frame(
        Estimate = est,
        `Std. Error` = NA_real_,
        `Lower 95%` = lower,
        `Upper 95%` = upper,
        df = df_val,
        `t value` = t_val,
        Pr = p_val,
        check.names = FALSE
      )
      row_names[k] <- if (ncol(Y_target) == 2) {
        "pcorr[rho]"
      } else {
        paste0("pcorr[", var_names[i], ",", var_names[j], "]")
      }
    }

    df_fixed <- do.call(rbind, rows)
    rownames(df_fixed) <- row_names

    pcorr_out <- if (ncol(Y_target) == 2) {
      stats::setNames(pcorr_values, "rho")
    } else {
      out <- diag(1, ncol(Y_target))
      colnames(out) <- rownames(out) <- var_names
      for (k in seq_len(ncol(pairs))) {
        i <- pairs[1, k]
        j <- pairs[2, k]
        out[i, j] <- out[j, i] <- pcorr_values[k]
      }
      out
    }

    fit <- Classic_Fit$new(
      model = self,
      par = list(pcorr = pcorr_out),
      log_lik = NA_real_,
      df_fixed = df_fixed,
      transform = list(pcorr = pcorr_out),
      generate = NULL,
      ci_method = "wald",
      laplace = FALSE,
      test_results = list(cor_test_method = method, partial = TRUE),
      view = view,
      df_method = method,
      show_df = TRUE
    )

    class(fit) <- unique(c(
      "rtmb_corr",
      "Classic_Fit",
      class(fit)
    ))

    return(fit)
  }

  rows <- vector("list", ncol(pairs))
  corr_values <- numeric(ncol(pairs))
  row_names <- character(ncol(pairs))

  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]
    j <- pairs[2, k]
    ok <- stats::complete.cases(Y_target[, i], Y_target[, j])
    x <- Y_target[ok, i]
    y <- Y_target[ok, j]
    ct <- suppressWarnings(stats::cor.test(x, y, method = method))
    est <- unname(ct$estimate)
    corr_values[k] <- est

    lower <- upper <- NA_real_
    if (!is.null(ct$conf.int)) {
      lower <- unname(ct$conf.int[1])
      upper <- unname(ct$conf.int[2])
    }

    row <- data.frame(
      Estimate = est,
      `Std. Error` = NA_real_,
      `Lower 95%` = lower,
      `Upper 95%` = upper,
      check.names = FALSE
    )

    if (method == "pearson") {
      row$df <- unname(ct$parameter %||% (length(x) - 2))
      row$`t value` <- unname(ct$statistic)
    } else {
      row$`S value` <- unname(ct$statistic)
    }
    row$Pr <- ct$p.value

    row_names[k] <- if (ncol(Y_target) == 2) {
      "corr[rho]"
    } else {
      paste0("corr[", var_names[i], ",", var_names[j], "]")
    }
    rows[[k]] <- row
  }

  df_fixed <- do.call(rbind, rows)
  rownames(df_fixed) <- row_names

  corr_out <- if (ncol(Y_target) == 2) {
    stats::setNames(corr_values, "rho")
  } else {
    out <- diag(1, ncol(Y_target))
    colnames(out) <- rownames(out) <- var_names
    for (k in seq_len(ncol(pairs))) {
      i <- pairs[1, k]
      j <- pairs[2, k]
      out[i, j] <- out[j, i] <- corr_values[k]
    }
    out
  }

  fit <- Classic_Fit$new(
    model = self,
    par = list(corr = corr_out),
    log_lik = NA_real_,
    df_fixed = df_fixed,
    transform = NULL,
    generate = NULL,
    ci_method = "wald",
    laplace = FALSE,
    test_results = list(cor_test_method = method),
    view = view,
    df_method = method,
    show_df = method == "pearson"
  )

  class(fit) <- unique(c(
    "rtmb_corr",
    "Classic_Fit",
    class(fit)
  ))

  fit
}

.build_classic_components <- function(self, private, raw, settings, se_method, num_samples, seed, info_log_lik = NULL) {
  # Extract components from raw
  ad_obj         <- raw$ad_obj
  opt            <- raw$opt
  sd_rep         <- raw$sd_rep
  unc_est_vec    <- raw$par_unc
  unc_se_vec     <- raw$se_unc
  Cov_u          <- raw$vcov_unc
  con_est_list   <- raw$par
  target_map     <- raw$map
  idx_fix_active <- raw$idx_fix_active
  idx_ran_full   <- raw$idx_ran_full
  idx_ran        <- raw$idx_ran
  factor_levels_seen <- raw$factor_levels_seen
  unc_est_list   <- raw$par_unc_list
  unc_se_list    <- raw$se_unc_list
  L_u_total      <- raw$L_u_total
  
  se_sampling <- (se_method == "sampling")
  
  # Initialize sensitivity variables for Satterthwaite/delta-method
  dH_list <- NULL
  V_opt <- NULL
  active_idx <- NULL
  dH_theta <- NULL
  V_theta <- NULL
  dH_beta <- NULL
  V_beta <- NULL
  
  # 1. Degrees of Freedom
  est_dfs_all <- rep(Inf, L_u_total)
  
  if (identical(settings$df_method, "numeric")) {
    est_dfs_all <- rep(settings$df_value, L_u_total)
  } else if (identical(settings$df_method, "inf")) {
    est_dfs_all <- rep(Inf, L_u_total)
  } else if (identical(settings$df_method, "residual")) {
    n_obs <- self$get_n_obs()
    
    if (is.null(n_obs) || length(n_obs) == 0 || is.na(n_obs) || n_obs <= 0) {
      if (!is.null(self$extra$nobs)) {
        n_obs <- self$extra$nobs
      } else if (!is.null(self$raw_data)) {
        n_obs <- nrow(as.data.frame(self$raw_data))
      } else if (!is.null(self$data$Y)) {
        n_obs <- length(self$data$Y)
      } else {
        n_obs <- NA_integer_
      }
    }
    
    if (is.na(n_obs) || n_obs <= 0) {
      warning("Could not determine residual degrees of freedom; using infinite df.", call. = FALSE)
      df_res <- Inf
    } else {
      # Calculate effective parameters (fixed effects)
      p_eff <- NA_integer_
      if (!is.null(self$type) && self$type %in% c("lm", "lmer") &&
          !is.null(self$formula) && !is.null(self$raw_data)) {
        p_eff <- tryCatch({
          mf <- stats::model.frame(nobars(self$formula), data = as.data.frame(self$raw_data))
          qr(stats::model.matrix(nobars(self$formula), data = mf))$rank
        }, error = function(e) NA_integer_)
      }
      if (is.na(p_eff)) {
        p_eff <- length(settings$target_vars)
      }
      if (p_eff == 0 && !is.null(self$extra$marginal)) {
        p_eff <- sum(vapply(self$extra$marginal, function(nm) {
          if (nm %in% names(self$par_list)) self$par_list[[nm]]$unc_length else 0
        }, numeric(1)))
      }
      
      if (p_eff == 0) {
        # Fallback to sum of lengths of fixed parameters
        p_eff <- sum(vapply(names(self$par_list), function(nm) {
          p <- self$par_list[[nm]]
          if (!isTRUE(p$random)) p$unc_length else 0
        }, numeric(1)))
      }
      df_res <- pmax(n_obs - p_eff, 1)
    }
    
    # Apply to all non-random parameters
    for (name in names(self$par_list)) {
      if (!isTRUE(self$par_list[[name]]$random)) {
        u_indices <- private$.get_unc_indices(name)
        if (length(u_indices) > 0) est_dfs_all[u_indices] <- df_res
      }
    }
    
    # Apply df_map if present (e.g. for mediation)
    if (!is.null(self$extra$df_map)) {
      df_map <- self$extra$df_map
      df_list_unc <- unconstrained_vector_to_list(est_dfs_all, self$par_list)
      for (v in names(df_list_unc)) {
        v_display <- sub("^b_c", "b", v)
        if (v %in% names(df_map)) {
          df_list_unc[[v]] <- rep(as.numeric(df_map[[v]]), length(df_list_unc[[v]]))
        } else if (v_display %in% names(df_map)) {
          df_list_unc[[v]] <- rep(as.numeric(df_map[[v_display]]), length(df_list_unc[[v]]))
        }
      }
      est_dfs_all <- unlist(df_list_unc, use.names = FALSE)
    }

    if (isTRUE(self$type == "corr") && identical(self$extra$corr_method %||% "reml", "reml")) {
      n_corr <- if (!is.null(self$data$Y) && is.matrix(self$data$Y)) {
        nrow(self$data$Y)
      } else {
        n_obs
      }
      if (!is.na(n_corr) && n_corr > 2) {
        p_y <- self$data$P_y %||% 0
        p_x <- if (!is.null(self$data$Y) && is.matrix(self$data$Y) && p_y > 0) {
          max(ncol(self$data$Y) - p_y, 0)
        } else {
          self$data$P_x %||% 0
        }
        df_reml <- list(
          mean = n_corr - 1,
          sd = n_corr - 2,
          corr = n_corr - 2,
          CF_corr = n_corr - 2,
          pcorr = pmax(n_corr - p_x - 2, 1)
        )
        df_list_unc <- unconstrained_vector_to_list(est_dfs_all, self$par_list)
        for (v in intersect(names(df_list_unc), names(df_reml))) {
          df_list_unc[[v]] <- rep(df_reml[[v]], length(df_list_unc[[v]]))
        }
        est_dfs_all <- unlist(df_list_unc, use.names = FALSE)

        if (is.null(self$extra$df_map)) self$extra$df_map <- list()
        self$extra$df_map$mean <- df_reml$mean
        self$extra$df_map$sd <- df_reml$sd
        self$extra$df_map$corr <- df_reml$corr
        self$extra$df_map$CF_corr <- df_reml$CF_corr
        self$extra$df_map$pcorr <- df_reml$pcorr
      }
    }
  } else if (identical(settings$df_method, "satterthwaite")) {
    satt_res <- self$calculate_satterthwaite_df(
      ad_obj, idx_fix_active, L_u_total, opt$par, 
      silent = FALSE, return_sensitivities = TRUE
    )
    if (is.list(satt_res)) {
      est_dfs_all <- satt_res$df
      dH_list <- satt_res$sensitivities
      V_opt <- satt_res$V
    } else {
      est_dfs_all <- satt_res
    }
    dH_theta <- dH_list
    V_theta <- V_opt
    
    if (settings$use_reml && length(settings$target_vars) > 0 && length(idx_ran) > 0) {
      full_par_names_mapped <- names(ad_obj$env$last.par)[idx_ran]
      target_ran_idx <- which(full_par_names_mapped %in% settings$target_vars)
      if (length(target_ran_idx) > 0) {
        reml_res <- self$calculate_reml_satterthwaite_df(
          ad_obj, opt$par, target_ran_idx, 
          silent = FALSE, return_sensitivities = TRUE
        )
        if (is.list(reml_res)) {
          active_idx <- idx_ran_full[target_ran_idx]
          est_dfs_all[active_idx] <- reml_res$df
          dH_beta <- reml_res$sensitivities
          V_beta <- reml_res$V_beta
        }
      }
    }
    # Important: Reset ad_obj after sensitivity runs
    ad_obj$fn(opt$par)
  }

  # 2. Standard Errors and Confidence Intervals
  con_se_list <- list()
  con_lower_list <- list()
  con_upper_list <- list()
  samps_con <- list()
  samps_tran <- list()
  samps_gq <- list()

  if (se_sampling) {
    # Ported from .inference_pipeline (simplified for readability)
    set.seed(seed)
    # Ensure PD covariance for sampling
    Cov_u_valid <- Cov_u[idx_fix_active, idx_fix_active, drop = FALSE]
    mu_valid <- unc_est_vec[idx_fix_active]
    eig <- eigen(Cov_u_valid, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-8)
    Cov_u_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    
    raw_samples <- MASS::mvrnorm(num_samples, mu = mu_valid, Sigma = Cov_u_pd)
    active_dfs <- est_dfs_all[idx_fix_active]
    
    for (j in seq_len(ncol(raw_samples))) {
      df_j <- active_dfs[j]
      if (is.finite(df_j)) {
        sigma_j <- sqrt(Cov_u_pd[j, j])
        z <- (raw_samples[, j] - mu_valid[j]) / (sigma_j + 1e-12)
        w <- rchisq(num_samples, df = df_j)
        raw_samples[, j] <- mu_valid[j] + z * sigma_j * sqrt(df_j / pmax(w, 1e-6))
      }
    }
    
    u_samples_mat <- matrix(rep(unc_est_vec, each = num_samples), nrow = num_samples)
    u_samples_mat[, idx_fix_active] <- raw_samples
    
    # Handle mapped parameters in samples
    idx_curr <- 1
    for (name in names(self$par_list)) {
      L_u <- self$par_list[[name]]$unc_length
      if (L_u > 0) {
        map_f <- if (!is.null(target_map[[name]])) target_map[[name]] else NULL
        for (i in seq_len(L_u)) {
          pos_full <- idx_curr + i - 1
          if (pos_full %in% idx_ran_full) next
          if (!is.null(map_f) && !is.na(map_f[i])) {
            lvl <- as.character(map_f[i])
            mapped_opt_idx <- factor_levels_seen[[lvl]]
            if (!is.null(mapped_opt_idx)) u_samples_mat[, pos_full] <- raw_samples[, mapped_opt_idx]
          }
        }
        idx_curr <- idx_curr + L_u
      }
    }
    
    # Handle random effects in samples (Laplace)
    if (settings$use_laplace && length(idx_ran_full) > 0) {
      for (ridx in idx_ran_full) {
        u_samples_mat[, ridx] <- rnorm(num_samples, mean = unc_est_vec[ridx], sd = unc_se_vec[ridx])
      }
    }
    
    # Transform samples to constrained space and derived variables
    local_data <- self$data
    local_par_list <- self$par_list
    local_transform <- self$transform
    local_generate <- self$generate
    
    for (s in seq_len(num_samples)) {
      u_vec <- u_samples_mat[s, ]
      tmp_con_list <- self$to_constrained(unconstrained_vector_to_list(u_vec, local_par_list))
      
      if (s == 1) {
        for (name in names(tmp_con_list)) samps_con[[name]] <- matrix(NA_real_, num_samples, length(as.numeric(tmp_con_list[[name]])))
      }
      for (name in names(tmp_con_list)) samps_con[[name]][s, ] <- as.numeric(tmp_con_list[[name]])
      
      if (!is.null(local_transform)) {
        tmp_tran <- tryCatch(local_transform(local_data, tmp_con_list), error = function(e) NULL)
        if (!is.null(tmp_tran)) {
          if (s == 1) {
            for (name in names(tmp_tran)) samps_tran[[name]] <- matrix(NA_real_, num_samples, length(as.numeric(tmp_tran[[name]])))
          }
          for (name in names(tmp_tran)) samps_tran[[name]][s, ] <- as.numeric(tmp_tran[[name]])
          tmp_con_list <- c(tmp_con_list, tmp_tran)
        }
      }
      
      if (!is.null(local_generate)) {
        tmp_gq <- tryCatch(local_generate(local_data, tmp_con_list), error = function(e) NULL)
        if (!is.null(tmp_gq)) {
          if (s == 1) {
            for (name in names(tmp_gq)) samps_gq[[name]] <- matrix(NA_real_, num_samples, length(as.numeric(tmp_gq[[name]])))
          }
          for (name in names(tmp_gq)) samps_gq[[name]][s, ] <- as.numeric(tmp_gq[[name]])
        }
      }
    }
    
    # Aggregate sample statistics
    for (name in names(self$par_list)) {
      mat <- samps_con[[name]]
      p_info <- self$par_list[[name]]
      con_se_list[[name]] <- apply(mat, 2, sd, na.rm = TRUE)
      con_lower_list[[name]] <- apply(mat, 2, quantile, 0.025, na.rm = TRUE)
      con_upper_list[[name]] <- apply(mat, 2, quantile, 0.975, na.rm = TRUE)
      if (length(p_info$dim) > 1) {
        dim(con_se_list[[name]]) <- dim(con_lower_list[[name]]) <- dim(con_upper_list[[name]]) <- p_info$dim
      }
    }
  } else {
    # Wald CI
    u_idx_current <- 1
    for (name in names(self$par_list)) {
      p_info <- self$par_list[[name]]
      u_val <- unc_est_list[[name]]
      u_se <- unc_se_list[[name]]
      c_val <- con_est_list[[name]]
      L_u <- p_info$unc_length
      L_c <- p_info$length
      u_indices <- if (L_u > 0) u_idx_current:(u_idx_current + L_u - 1) else integer(0)
      u_idx_current <- u_idx_current + L_u
      
      # Delta method for SE
      transform_single <- function(u) {
        tmp_list <- unc_est_list
        tmp_list[[name]] <- u
        as.numeric(self$to_constrained(tmp_list)[[name]])
      }
      
      J <- matrix(0, L_c, L_u)
      for (i in seq_len(L_u)) {
        u_tmp <- u_val
        u_tmp[i] <- u_tmp[i] + 1e-5
        J[, i] <- (transform_single(u_tmp) - as.numeric(c_val)) / 1e-5
      }
      
      c_se <- numeric(L_c)
      if (L_u > 0) {
        Cov_sub <- Cov_u[u_indices, u_indices, drop = FALSE]
        for (j in seq_len(L_c)) {
          c_se[j] <- sqrt(sum((J[j, ] %*% Cov_sub) * J[j, ], na.rm = TRUE))
        }
      }
      con_se_list[[name]] <- if (length(p_info$dim) > 1) structure(c_se, dim = p_info$dim) else c_se
      
      # CI using t-distribution
      u_dfs <- est_dfs_all[u_indices]
      u_q_95 <- ifelse(is.na(u_dfs) | is.infinite(u_dfs), qnorm(0.975), qt(0.975, df = pmax(u_dfs, 2.1)))
      
      if (p_info$bounds == "corr_matrix" || p_info$type == "corr_matrix") {
        # Fisher's z transform for correlations
        rho_val <- as.numeric(c_val)
        rho_clipped <- pmax(pmin(rho_val, 0.9999), -0.9999)
        z_val <- atanh(rho_clipped)
        z_se <- c_se / (1 - rho_clipped^2 + 1e-10)
        u_q_95_m <- if (any(is.na(u_dfs) | is.infinite(u_dfs))) qnorm(0.975) else qt(0.975, df = pmax(mean(u_dfs, na.rm = TRUE), 2.1))
        c_low <- tanh(z_val - u_q_95_m * z_se)
        c_up <- tanh(z_val + u_q_95_m * z_se)
        is_fixed <- c_se < 1e-10 & abs(rho_val - 1) < 1e-8
        if (any(is_fixed)) { c_low[is_fixed] <- 1.0; c_up[is_fixed] <- 1.0 }
      } else {
        u_low <- u_val - u_q_95 * u_se
        u_up <- u_val + u_q_95 * u_se
        c_low_raw <- transform_single(u_low)
        c_up_raw <- transform_single(u_up)
        c_low <- pmin(c_low_raw, c_up_raw)
        c_up <- pmax(c_low_raw, c_up_raw)
      }
      con_lower_list[[name]] <- if (length(p_info$dim) > 1) structure(c_low, dim = p_info$dim) else c_low
      con_upper_list[[name]] <- if (length(p_info$dim) > 1) structure(c_up, dim = p_info$dim) else c_up
    }
  }

  # 3. Build Summary Tables
  df_list_con <- list()
  df_list_unc <- unconstrained_vector_to_list(est_dfs_all, self$par_list)
  for (name in names(self$par_list)) {
    p <- self$par_list[[name]]
    if (p$unc_length == p$length) {
      df_list_con[[name]] <- df_list_unc[[name]]
    } else {
      df_val <- rep(Inf, p$length)
      if (p$unc_length > 0) {
        m_df <- mean(df_list_unc[[name]], na.rm = TRUE)
        if (is.finite(m_df)) df_val[] <- m_df
      }
      df_list_con[[name]] <- df_val
    }
  }

  build_summary_df <- function(target_random = FALSE) {
    names_vec <- c(); est_vec <- c(); se_vec <- c(); low_vec <- c(); up_vec <- c(); df_vec <- c()
    for (name in names(self$par_list)) {
      p_info <- self$par_list[[name]]
      if (isTRUE(p_info$random) == target_random) {
        f_names <- generate_flat_names(name, p_info$dim, self$par_names[[name]])
        names_vec <- c(names_vec, f_names)
        est_vec <- c(est_vec, as.numeric(con_est_list[[name]]))
        se_vec <- c(se_vec, as.numeric(con_se_list[[name]]))
        low_vec <- c(low_vec, as.numeric(con_lower_list[[name]]))
        up_vec <- c(up_vec, as.numeric(con_upper_list[[name]]))
        df_vec <- c(df_vec, if (!is.null(df_list_con[[name]])) as.numeric(df_list_con[[name]]) else rep(Inf, p_info$length))
      }
    }
    if (length(names_vec) == 0) return(NULL)
    res_df <- data.frame(
      Estimate = est_vec, 
      `Std. Error` = se_vec, 
      `Lower 95%` = low_vec, 
      `Upper 95%` = up_vec, 
      check.names = FALSE
    )
    if (settings$show_df) res_df$df <- df_vec
    rownames(res_df) <- names_vec
    return(res_df)
  }
  
  df_fixed <- build_summary_df(FALSE)
  df_random <- if (settings$use_laplace) build_summary_df(TRUE) else NULL

  # 4. Transform and Generate Tables
  tran_list <- if (!is.null(self$transform)) tryCatch(self$transform(self$data, con_est_list), error = function(e) NULL) else NULL
  gq_list <- if (!is.null(self$generate)) {
    full_con <- if (!is.null(tran_list)) c(con_est_list, tran_list) else con_est_list
    tryCatch(self$generate(self$data, full_con), error = function(e) NULL)
  } else NULL

  # Helper for derived summaries (aligned with .inference_pipeline)
  build_derived_summary <- function(
    func, base_out, is_generate = FALSE, 
    dH_list_in = NULL, V_opt_in = NULL, is_reml_in = FALSE, 
    active_idx_in = NULL, dH_theta_in = NULL, V_theta_in = NULL, 
    theta_idx_in = NULL
  ) {
    if (is.null(func) || is.null(base_out) || length(base_out) == 0) return(NULL)
    flat_base <- unlist(base_out, use.names = FALSE)
    L_out <- length(flat_base)
    names_vec <- c()
    for (name in names(base_out)) {
      names_vec <- c(names_vec, generate_flat_names(name, if (is.null(dim(base_out[[name]]))) length(base_out[[name]]) else dim(base_out[[name]]), self$par_names[[name]]))
    }
    
    # J for Delta method
    J <- matrix(0, L_out, L_u_total)
    for (i in seq_len(L_u_total)) {
      u_tmp <- unc_est_vec
      u_tmp[i] <- u_tmp[i] + 1e-5
      tmp_con <- self$to_constrained(unconstrained_vector_to_list(u_tmp, self$par_list))
      if (is_generate && !is.null(self$transform)) {
        ut <- tryCatch(self$transform(self$data, tmp_con), error = function(e) NULL)
        if (!is.null(ut)) tmp_con <- c(tmp_con, ut)
      }
      tmp_out <- tryCatch(func(self$data, tmp_con), error = function(e) NULL)
      if (!is.null(tmp_out)) J[, i] <- (unlist(tmp_out, use.names = FALSE) - flat_base) / 1e-5
    }
    
    # 1. Degrees of Freedom for derived quantities
    derived_dfs <- rep(Inf, L_out)
    if (settings$show_df) {
      if (identical(settings$df_method, "satterthwaite") && !is.null(dH_theta_in) && !is.null(V_theta_in)) {
        # Use the Satterthwaite delta method helper
        derived_dfs <- .satterthwaite_df_delta(
          J_full = J,
          V_theta = V_theta_in,
          dH_theta = dH_theta_in,
          idx_theta_full = theta_idx_in %||% idx_fix_active,
          V_beta = if (is_reml_in) V_opt_in else NULL,
          dH_beta = if (is_reml_in) dH_list_in else NULL,
          idx_beta_full = if (is_reml_in) active_idx_in else NULL,
          max_df = NULL
        )
      } else if (identical(settings$df_method, "residual")) {
        for (j in seq_len(L_out)) {
          c_dfs <- est_dfs_all[abs(J[j, ]) > 1e-8]
          if (length(c_dfs) > 0) derived_dfs[j] <- min(c_dfs)
        }
      } else if (identical(settings$df_method, "numeric")) {
        derived_dfs[] <- settings$df_value
      }
    }

    if (se_sampling) {
      samps_list <- if (is_generate) samps_gq else samps_tran
      target_names <- names(base_out)
      samps_sub <- samps_list[target_names]
      samps_sub <- samps_sub[!vapply(samps_sub, is.null, logical(1))]

      if (length(samps_sub) == 0) {
        se_out <- low_out <- up_out <- rep(NA_real_, L_out)
      } else {
        mat_all <- do.call(cbind, unname(samps_sub))
        if (ncol(mat_all) != L_out) {
          se_out <- low_out <- up_out <- rep(NA_real_, L_out)
        } else {
          se_out <- apply(mat_all, 2, sd, na.rm = TRUE)
          low_out <- apply(mat_all, 2, quantile, 0.025, na.rm = TRUE)
          up_out <- apply(mat_all, 2, quantile, 0.975, na.rm = TRUE)
        }
      }
    } else {
      se_out <- sqrt(abs(rowSums((J %*% Cov_u) * J)))
      
      # Use t-critical if finite df is available
      crit <- ifelse(is.finite(derived_dfs), 
                     qt(0.975, df = pmax(derived_dfs, 2.1)), 
                     qnorm(0.975))
      low_out <- flat_base - crit * se_out
      up_out <- flat_base + crit * se_out
    }
    
    res_df <- data.frame(
      Estimate = flat_base, 
      `Std. Error` = se_out, 
      `Lower 95%` = low_out, 
      `Upper 95%` = up_out, 
      check.names = FALSE
    )
    if (settings$show_df) {
      res_df$df <- derived_dfs
    }
    rownames(res_df) <- names_vec
    return(res_df)
  }

  df_transform <- build_derived_summary(
    self$transform, tran_list, FALSE,
    if (settings$use_reml) dH_beta else dH_list,
    if (settings$use_reml) V_beta else V_opt,
    settings$use_reml,
    if (settings$use_reml) active_idx else NULL,
    if (settings$use_reml) dH_theta else NULL,
    if (settings$use_reml) V_theta else NULL,
    if (settings$use_reml) idx_fix_active else NULL
  )

  df_generate <- build_derived_summary(
    self$generate, gq_list, TRUE,
    if (settings$use_reml) dH_beta else dH_list,
    if (settings$use_reml) V_beta else V_opt,
    settings$use_reml,
    if (settings$use_reml) active_idx else NULL,
    if (settings$use_reml) dH_theta else NULL,
    if (settings$use_reml) V_theta else NULL,
    if (settings$use_reml) idx_fix_active else NULL
  )

  apply_df_map_to_summary <- function(df) {
    if (is.null(df) || is.null(self$extra$df_map) || !("df" %in% names(df))) {
      return(df)
    }
    rn <- rownames(df)
    if (is.null(rn)) return(df)
    base_names <- gsub("\\[.*\\]$", "", rn)
    for (i in seq_along(base_names)) {
      if (base_names[i] %in% names(self$extra$df_map)) {
        df$df[i] <- as.numeric(self$extra$df_map[[base_names[i]]])
      }
    }
    df
  }

  df_fixed <- apply_df_map_to_summary(df_fixed)
  df_random <- apply_df_map_to_summary(df_random)
  df_transform <- apply_df_map_to_summary(df_transform)
  df_generate <- apply_df_map_to_summary(df_generate)

  # 5. Log-likelihood and Vcov
  log_lik <- if (!is.null(info_log_lik)) info_log_lik else -opt$objective
  rss <- NULL
  df_residual <- NULL
  if (identical(self$type, "lm") && !is.null(con_est_list$sigma) && !is.null(self$data$Y)) {
    n_y <- length(as.numeric(self$data$Y))
    p_fixed <- if (is.null(df_fixed)) 0L else sum(gsub("\\[.*\\]$", "", rownames(df_fixed)) != "sigma")
    df_residual <- n_y - p_fixed
    if (is.finite(df_residual) && df_residual > 0) {
      rss <- as.numeric(con_est_list$sigma)^2 * df_residual
    }
  }


  # Build V_fixed_full (vcov)
  if (is.null(df_fixed) || nrow(df_fixed) == 0) {
    V_fixed_full <- matrix(0, 0, 0)
  } else {
    L_fixed <- nrow(df_fixed)
    J_fixed <- matrix(0, L_fixed, L_u_total)
    base_f <- df_fixed$Estimate
    for (i in seq_len(L_u_total)) {
      u_tmp <- unc_est_vec
      u_tmp[i] <- u_tmp[i] + 1e-5
      tmp_con <- self$to_constrained(unconstrained_vector_to_list(u_tmp, self$par_list))
      tmp_vec <- c()
      for (name in names(self$par_list)) {
        if (!isTRUE(self$par_list[[name]]$random)) {
          tmp_vec <- c(tmp_vec, as.numeric(tmp_con[[name]]))
        }
      }
      if (length(tmp_vec) == L_fixed) J_fixed[, i] <- (tmp_vec - base_f) / 1e-5
    }
    V_fixed_full <- J_fixed %*% Cov_u %*% t(J_fixed)
  }
  rownames(V_fixed_full) <- colnames(V_fixed_full) <- rownames(df_fixed)

  # 6. Additional test results (e.g. for table)
  test_results <- list()
  if (isTRUE(self$type == "table") && !is.null(self$extra$tab)) {
    tab <- self$extra$tab
    test_results$chisq <- stats::chisq.test(tab)
    test_results$fisher <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
  }

  list(
    df_fixed = df_fixed,
    random_effects = df_random,
    df_transform = df_transform,
    df_generate = df_generate,
    transform = tran_list,
    generate = gq_list,
    se_samples = if (se_sampling) list(con = samps_con, tran = samps_tran, gq = samps_gq) else NULL,
    vcov = V_fixed_full,
    log_lik = log_lik,
    rss = rss,
    df_residual = df_residual,
    test_results = test_results
  )
}
