.inference_pipeline <- function(self, private, laplace, init, num_estimate, control,
                               optimizer, method, map, fixed, se, se_method,
                               num_samples, seed, marginal, df_method, view, ...) {

  se_method <- match.arg(se_method, c("wald", "sampling", "none"))
  se_enabled <- !identical(se_method, "none")
  se <- isTRUE(se) && se_enabled

  if (is.null(marginal)) marginal <- "none"
  target_vars <- private$.resolve_marginal_vars(marginal)
  marginal_flag <- length(target_vars) > 0L

  dot_args <- list(...)
  requested_contrasts <- dot_args$contrasts %||% self$contrasts
  if (!is.null(requested_contrasts)) {
    old_contrasts <- options()$contrasts
    on.exit(options(contrasts = old_contrasts), add = TRUE)
    if (requested_contrasts == "sum") options(contrasts = c("contr.sum", "contr.poly"))
    else if (requested_contrasts == "treatment") options(contrasts = c("contr.treatment", "contr.poly"))
    self$contrasts <- requested_contrasts
  }

  # Delegate to the corrected apply_constraints
  if (!is.null(fixed)) {
    new_model <- self$fixed_model(fixed)

    return(.inference_pipeline(
      self = new_model, private = new_model$.__enclos_env__$private,
      laplace = laplace, init = new_model$init, num_estimate = num_estimate, control = control,
      optimizer = optimizer, method = method, map = map, fixed = NULL,
      se = se, se_method = se_method, num_samples = num_samples, seed = seed,
      marginal = marginal, df_method = df_method, view = view, ...
    ))
  }

  if (se_method == "sampling") se_sampling <- TRUE else se_sampling <- FALSE

  # Resolve DF settings using helper
  df_settings <- .resolve_optimize_df_settings(df_method, marginal_flag, se_method, marginal)
  auto_df <- df_settings$auto_df
  show_df <- df_settings$show_df
  df_t <- df_settings$df_t
  df_method_l <- df_settings$df_method_norm

  # --- Warnings for ignored df_method ---
  if (!se_enabled) {
    explicit_df_requested <- is.numeric(df_method) || (is.character(df_method) && !(tolower(df_method) %in% c("auto", "inf", "infinite", "none")))
    if (explicit_df_requested) {
      warning("df_method is ignored when se_method = 'none'.", call. = FALSE)
    }
  } else if (!marginal_flag) {
    explicit_df_requested <- is.numeric(df_method) || (is.character(df_method) && !(tolower(df_method) %in% c("auto", "inf", "infinite", "none")))
    if (explicit_df_requested) {
      if (identical(marginal, "auto")) {
        warning("df_method is ignored because marginal = 'auto' resolved to no parameters. Set model$extra$marginal or specify marginal = c(...) explicitly.", call. = FALSE)
      } else if (identical(marginal, "none")) {
        warning("df_method is ignored when marginal = 'none' in optimize(). Use marginal to enable finite-df adjustment, or set df_method = 'none'.", call. = FALSE)
      }
    }
  }



  # --- Logic for using the common RTMB fitting engine ---
  raw <- private$.fit_rtmb(
    laplace = laplace, init = init, num_estimate = num_estimate, control = control,
    optimizer = optimizer, method = method, map = map, se = se,
    reml = marginal_flag, target_vars = target_vars,
    jacobian_target = NULL,  # Let .fit_rtmb decide after marginalization handling
    apply_prior_correction = TRUE
  )

  # Extract results and essential dimensions
  ad_setup       <- raw$ad_setup
  ad_obj         <- raw$ad_obj
  opt            <- raw$opt
  sd_rep         <- if (se_enabled) raw$sd_rep else NULL
  unc_est_vec    <- raw$par_unc
  L_u_total      <- raw$L_u_total
  unc_se_vec     <- if (se_enabled) raw$se_unc else rep(NA_real_, L_u_total)
  Cov_u          <- if (se_enabled) raw$vcov_unc else matrix(NA_real_, L_u_total, L_u_total)
  con_est_list   <- raw$par
  target_map     <- raw$map
  idx_fix_active <- raw$idx_fix_active
  idx_ran_full   <- raw$idx_ran_full
  idx_ran        <- raw$idx_ran
  factor_levels_seen <- raw$factor_levels_seen
  unc_est_list   <- raw$par_unc_list
  unc_se_list    <- raw$se_unc_list
  opt_history    <- raw$opt_history
  fallback_needed <- raw$fallback_needed

  con_se_list <- list(); con_lower_list <- list(); con_upper_list <- list(); samps_con <- list(); samps_tran <- list(); samps_gq <- list()

  dH_list <- NULL; V_opt <- NULL; active_idx <- NULL; dH_theta <- dH_beta <- V_theta <- V_beta <- NULL
  if (any(!is.infinite(df_t))) {
    est_dfs_all <- rep(df_t[1], L_u_total)
  } else if (auto_df) {
    # Perform Satterthwaite approximation
    satt_res <- self$calculate_satterthwaite_df(ad_obj, idx_fix_active, L_u_total, opt$par, max_df = 1e6, silent = FALSE, return_sensitivities = TRUE)
    if (is.list(satt_res)) {
      est_dfs_all <- satt_res$df
      dH_list <- satt_res$sensitivities
      V_opt <- satt_res$V
    } else {
      est_dfs_all <- satt_res
    }
    dH_theta <- dH_list
    V_theta <- V_opt

    if (marginal_flag && length(ad_obj$env$random) > 0) {
      full_par_names_mapped <- names(ad_obj$env$last.par)[idx_ran]
      target_ran_idx <- which(full_par_names_mapped %in% target_vars)
      if (length(target_ran_idx) > 0) {
        reml_res <- self$calculate_reml_satterthwaite_df(ad_obj, opt$par, target_ran_idx, silent = FALSE, return_sensitivities = TRUE)
        if (is.list(reml_res)) {
          active_idx <- idx_ran_full[target_ran_idx]
          est_dfs_all[active_idx] <- reml_res$df
          dH_beta <- reml_res$sensitivities
          V_beta <- reml_res$V_beta
        }
      }
    }
    ad_obj$fn(opt$par)
  } else {
    est_dfs_all <- rep(Inf, L_u_total)
  }

  if (se_sampling && se_enabled) {
    set.seed(seed)
    if (!getOption("BayesRTMB.silent", FALSE)) cat(sprintf("Using simulation-based error propagation (%d samples)...\n", num_samples))
    Cov_u_valid <- Cov_u[idx_fix_active, idx_fix_active, drop = FALSE]; mu_valid <- unc_est_vec[idx_fix_active]
    eig <- eigen(Cov_u_valid, symmetric = TRUE); eig$values <- pmax(eig$values, 1e-8); Cov_u_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
    raw_samples <- MASS::mvrnorm(num_samples, mu = mu_valid, Sigma = Cov_u_pd); active_dfs <- est_dfs_all[idx_fix_active]
    for (j in 1:ncol(raw_samples)) {
      df_j <- active_dfs[j]
      if (is.finite(df_j)) { sigma_j <- sqrt(Cov_u_pd[j, j]); z <- (raw_samples[, j] - mu_valid[j]) / (sigma_j + 1e-12); w <- rchisq(num_samples, df = df_j); raw_samples[, j] <- mu_valid[j] + z * sigma_j * sqrt(df_j / pmax(w, 1e-6)) }
    }
    u_samples_mat <- matrix(rep(unc_est_vec, each = num_samples), nrow = num_samples); u_samples_mat[, idx_fix_active] <- raw_samples
    idx_curr <- 1
    for (name in names(self$par_list)) {
      L_u <- self$par_list[[name]]$unc_length
      if (L_u > 0) {
        map_f <- if (!is.null(target_map[[name]])) target_map[[name]] else NULL
        for (i in 1:L_u) { pos_last <- idx_curr + i - 1; if (pos_last %in% idx_ran_full) next; if (!is.null(map_f) && !is.na(map_f[i])) { lvl <- as.character(map_f[i]); mapped_opt_idx <- factor_levels_seen[[lvl]]; if (!is.null(mapped_opt_idx)) u_samples_mat[, pos_last] <- raw_samples[, mapped_opt_idx] } }
        idx_curr <- idx_curr + L_u
      }
    }
    if (laplace && length(idx_ran_full) > 0) for (ridx in idx_ran_full) u_samples_mat[, ridx] <- rnorm(num_samples, mean = unc_est_vec[ridx], sd = unc_se_vec[ridx])

    local_data <- self$data; local_par_list <- self$par_list; local_transform <- self$transform; local_generate <- self$generate
    test_con_list <- self$to_constrained(unconstrained_vector_to_list(unc_est_vec, local_par_list, map = target_map))
    for (name in names(test_con_list)) samps_con[[name]] <- matrix(NA_real_, num_samples, length(as.numeric(test_con_list[[name]])))
    if (!is.null(local_transform)) {
      test_tran <- local_transform(local_data, test_con_list)
      if (!is.null(test_tran)) { for (name in names(test_tran)) samps_tran[[name]] <- matrix(NA_real_, num_samples, length(as.numeric(test_tran[[name]]))); test_con_list <- c(test_con_list, test_tran) }
    }
    if (!is.null(local_generate)) {
      test_gq <- local_generate(local_data, test_con_list)
      if (!is.null(test_gq)) for (name in names(test_gq)) samps_gq[[name]] <- matrix(NA_real_, num_samples, length(as.numeric(test_gq[[name]])))
    }

    process_sample <- function(s, u_vec) {
      tmp_con_list <- self$to_constrained(unconstrained_vector_to_list(u_vec, local_par_list))
      for (name in names(tmp_con_list)) samps_con[[name]][s, ] <<- as.numeric(tmp_con_list[[name]])
      if (!is.null(local_transform)) {
        tmp_tran <- local_transform(local_data, tmp_con_list); if (!is.null(tmp_tran)) { for (name in names(tmp_tran)) samps_tran[[name]][s, ] <<- as.numeric(tmp_tran[[name]]); tmp_con_list <- c(tmp_con_list, tmp_tran) }
      }
      if (!is.null(local_generate)) {
        tmp_gq <- local_generate(local_data, tmp_con_list); if (!is.null(tmp_gq)) for (name in names(tmp_gq)) samps_gq[[name]][s, ] <<- as.numeric(tmp_gq[[name]])
      }
    }
    for (s in 1:num_samples) process_sample(s, u_samples_mat[s, ])

    for (name in names(self$par_list)) {
      mat <- samps_con[[name]]; p_info <- self$par_list[[name]]
      con_se_list[[name]] <- if (is.null(mat)) rep(NA, p_info$length) else apply(mat, 2, sd, na.rm = TRUE)
      con_lower_list[[name]] <- if (is.null(mat)) rep(NA, p_info$length) else apply(mat, 2, quantile, 0.025, na.rm = TRUE)
      con_upper_list[[name]] <- if (is.null(mat)) rep(NA, p_info$length) else apply(mat, 2, quantile, 0.975, na.rm = TRUE)
      if (length(p_info$dim) > 1) { dim(con_se_list[[name]]) <- dim(con_lower_list[[name]]) <- dim(con_upper_list[[name]]) <- p_info$dim }
    }
  } else if (se_enabled) {
    u_idx_current <- 1
    for (name in names(self$par_list)) {
      p_info <- self$par_list[[name]]; u_val <- unc_est_list[[name]]; u_se <- unc_se_list[[name]]; c_val <- con_est_list[[name]]
      L_u <- p_info$unc_length; L_c <- p_info$length; u_indices <- if (L_u > 0) u_idx_current:(u_idx_current + L_u - 1) else integer(0); u_idx_current <- u_idx_current + L_u
      transform_single <- function(u) { tmp_list <- unc_est_list; tmp_list[[name]] <- u; as.numeric(self$to_constrained(tmp_list)[[name]]) }
      J <- matrix(0, L_c, L_u); for (i in seq_len(L_u)) { u_tmp <- u_val; u_tmp[i] <- u_tmp[i] + 1e-5; J[, i] <- (transform_single(u_tmp) - as.numeric(c_val)) / 1e-5 }
      c_se <- numeric(L_c); if (L_u > 0) { Cov_sub <- Cov_u[u_indices, u_indices, drop = FALSE]; for (j in 1:L_c) c_se[j] <- sqrt(sum((J[j, ] %*% Cov_sub) * J[j, ], na.rm = TRUE)) }
      con_se_list[[name]] <- if (length(p_info$dim) > 1) structure(c_se, dim = p_info$dim) else c_se
      u_dfs <- est_dfs_all[u_indices]; u_q_95 <- ifelse(is.na(u_dfs) | is.infinite(u_dfs), qnorm(0.975), qt(0.975, df = pmax(u_dfs, 2.1)))
      if (p_info$bounds == "corr_matrix" || p_info$type == "corr_matrix") {
        rho_val <- as.numeric(c_val); rho_clipped <- pmax(pmin(rho_val, 0.9999), -0.9999); z_val <- atanh(rho_clipped); z_se <- c_se / (1 - rho_clipped^2 + 1e-10); u_q_95_m <- if (any(is.na(u_dfs) | is.infinite(u_dfs))) qnorm(0.975) else qt(0.975, df = pmax(mean(u_dfs, na.rm = TRUE), 2.1))
        c_low <- tanh(z_val - u_q_95_m * z_se); c_up <- tanh(z_val + u_q_95_m * z_se); is_fixed <- c_se < 1e-10 & abs(rho_val - 1) < 1e-8; if (any(is_fixed)) { c_low[is_fixed] <- 1.0; c_up[is_fixed] <- 1.0 }
      } else {
        u_low <- u_val - u_q_95 * u_se; u_up <- u_val + u_q_95 * u_se; c_low_raw <- transform_single(u_low); c_up_raw <- transform_single(u_up); c_low <- pmin(c_low_raw, c_up_raw); c_up <- pmax(c_low_raw, c_up_raw)
      }
      con_lower_list[[name]] <- if (length(p_info$dim) > 1) structure(c_low, dim = p_info$dim) else c_low
      con_upper_list[[name]] <- if (length(p_info$dim) > 1) structure(c_up, dim = p_info$dim) else c_up
    }
  } else {
    # SE disabled: use point estimates and NA for uncertainty
    for (name in names(self$par_list)) {
      p_info <- self$par_list[[name]]
      L_c <- p_info$length
      con_se_list[[name]] <- rep(NA_real_, L_c)
      con_lower_list[[name]] <- rep(NA_real_, L_c)
      con_upper_list[[name]] <- rep(NA_real_, L_c)
      if (length(p_info$dim) > 1) {
        dim(con_se_list[[name]]) <- dim(con_lower_list[[name]]) <- dim(con_upper_list[[name]]) <- p_info$dim
      }
    }
  }


  df_list_con <- list(); df_list_unc <- unconstrained_vector_to_list(est_dfs_all, self$par_list)
  for (name in names(self$par_list)) {
    p <- self$par_list[[name]]; if (p$unc_length == p$length) df_list_con[[name]] <- df_list_unc[[name]]
    else { df_val <- rep(Inf, p$length); if (p$unc_length > 0) { m_df <- mean(df_list_unc[[name]], na.rm = TRUE); if (is.finite(m_df)) df_val[] <- m_df }; df_list_con[[name]] <- df_val }
  }

  .apply_df_map <- function(df, df_map) {
    if (is.null(df) || is.null(df_map) || length(df_map) == 0 || !("df" %in% names(df))) return(df)

    rn <- rownames(df)
    if (is.null(rn)) return(df)

    for (i in seq_along(rn)) {
      nm <- rn[i]
      base <- sub("\\[.*$", "", nm)

      if (nm %in% names(df_map)) {
        df$df[i] <- as.numeric(df_map[[nm]])
      } else if (base %in% names(df_map)) {
        df$df[i] <- as.numeric(df_map[[base]])
      }
    }
    df
  }

  build_summary_df <- function(target_random = FALSE) {
    names_vec <- c(); est_vec <- c(); se_vec <- c(); low_vec <- c(); up_vec <- c(); df_vec <- c()
    for (name in names(self$par_list)) {
      p_info <- self$par_list[[name]]; if (isTRUE(p_info$random) == target_random) {
        f_names <- generate_flat_names(name, p_info$dim, self$par_names[[name]]); names_vec <- c(names_vec, f_names); est_vec <- c(est_vec, as.numeric(con_est_list[[name]])); se_vec <- c(se_vec, as.numeric(con_se_list[[name]])); low_vec <- c(low_vec, as.numeric(con_lower_list[[name]])); up_vec <- c(up_vec, as.numeric(con_upper_list[[name]])); df_vec <- c(df_vec, if (!is.null(df_list_con[[name]])) as.numeric(df_list_con[[name]]) else rep(Inf, p_info$length))
      }
    }
    if (length(names_vec) == 0) return(NULL)
    res_df <- data.frame(Estimate = est_vec, `Std. Error` = se_vec, `Lower 95%` = low_vec, `Upper 95%` = up_vec, check.names = FALSE); if (show_df) res_df$df <- df_vec; rownames(res_df) <- names_vec; return(res_df)
  }
  df_fixed <- build_summary_df(FALSE); df_random <- if (laplace) build_summary_df(TRUE) else NULL

  df_map <- self$extra$df_map
  if (!is.null(df_map)) {
    df_fixed    <- .apply_df_map(df_fixed, df_map)
    df_random   <- .apply_df_map(df_random, df_map)
  }
  con_est_vec <- unlist(con_est_list, use.names = FALSE); tran_list <- if (!is.null(self$transform)) tryCatch(self$transform(self$data, con_est_list), error = function(e) NULL) else NULL
  gq_list <- if (!is.null(self$generate)) tryCatch(self$generate(self$data, if (!is.null(tran_list)) c(con_est_list, tran_list) else con_est_list), error = function(e) NULL) else NULL

  build_derived_summary <- function(func, base_out, is_generate = FALSE, dH_list_in = NULL, V_opt_in = NULL, is_reml_in = FALSE, active_idx_in = NULL, dH_theta_in = NULL, V_theta_in = NULL, theta_idx_in = NULL) {
    if (is.null(func) || is.null(base_out) || length(base_out) == 0) return(NULL)
    flat_base <- unlist(base_out, use.names = FALSE); L_out <- length(flat_base); names_vec <- c()
    for (name in names(base_out)) names_vec <- c(names_vec, generate_flat_names(name, if (is.null(dim(base_out[[name]]))) length(base_out[[name]]) else dim(base_out[[name]]), self$par_names[[name]]))
    
    # 1. Delta method Jacobian calculation (needed for both SE and DF if not sampling)
    J <- NULL
    if (se_enabled && (show_df || !se_sampling)) {
      J <- matrix(0, L_out, L_u_total)
      for (i in 1:L_u_total) {
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
    }

    # 2. Degrees of Freedom
    derived_dfs <- rep(Inf, L_out)
    if (show_df && se_enabled && auto_df && !is.null(dH_theta_in) && !is.null(V_theta_in) && !is.null(J)) {
       # Use the new Satterthwaite delta method helper
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
    } else if (show_df && !is.infinite(df_t[1])) {
       derived_dfs <- rep(df_t[1], L_out)
    }

    # 3. Uncertainty (SE/CI)
    if (!se_enabled) {
      se_out <- low_out <- up_out <- rep(NA_real_, L_out)
    } else if (se_sampling) {
      # FIXED Sampling Branch: match samples to specific output names
      samps_list <- if (is_generate) samps_gq else samps_tran
      target_names <- names(base_out)
      samps_sub <- samps_list[target_names]
      samps_sub <- samps_sub[!vapply(samps_sub, is.null, logical(1))]
      
      if (length(samps_sub) == 0L) {
        se_out <- low_out <- up_out <- rep(NA_real_, L_out)
      } else {
        mat_all <- do.call(cbind, unname(samps_sub))
        if (ncol(mat_all) != L_out) {
          se_out <- low_out <- up_out <- rep(NA_real_, L_out)
        } else {
          se_out <- apply(mat_all, 2, stats::sd, na.rm = TRUE)
          low_out <- apply(mat_all, 2, stats::quantile, 0.025, na.rm = TRUE)
          up_out <- apply(mat_all, 2, stats::quantile, 0.975, na.rm = TRUE)
        }
      }
    } else {
      # delta method SE
      se_out <- rep(NA_real_, L_out)
      smry_rep <- if (!is.null(sd_rep)) {
        tryCatch(
          suppressWarnings(as.data.frame(summary(sd_rep, select = "report"))),
          error = function(e) NULL
        )
      } else {
        NULL
      }
      if (!is.null(smry_rep)) {
        m_idx <- match(names_vec, rownames(smry_rep))
        valid <- !is.na(m_idx)
        if (any(valid)) se_out[valid] <- smry_rep$`Std. Error`[m_idx[valid]]
      }
      for (j in which(is.na(se_out))) {
        se_out[j] <- sqrt(abs(sum((J[j, ] %*% Cov_u) * J[j, ])))
      }
      
      crit <- ifelse(is.finite(derived_dfs), 
                     qt(0.975, df = pmax(derived_dfs, 2.1)), 
                     qnorm(0.975))
      low_out <- flat_base - crit * se_out
      up_out <- flat_base + crit * se_out
      
      corr_pat <- "^(corr|pcorr|B_corr|W_corr)(\\[|$)"
      is_corr <- grepl(corr_pat, names_vec)
      for (j in which(is_corr)) {
        r <- pmax(pmin(flat_base[j], 0.9999), -0.9999)
        z <- atanh(r)
        se_z <- se_out[j] / (1 - r^2 + 1e-10)
        
        cj <- crit[j]
        low_out[j] <- tanh(z - cj * se_z)
        up_out[j] <- tanh(z + cj * se_z)
      }
    }
    
    res_df <- data.frame(Estimate = flat_base, `Std. Error` = se_out, `Lower 95%` = low_out, `Upper 95%` = up_out, check.names = FALSE)
    if (show_df) {
      res_df$df <- derived_dfs
    }
    rownames(res_df) <- names_vec
    return(res_df)
  }

  df_transform <- build_derived_summary(self$transform, tran_list, FALSE, if (marginal_flag) dH_beta else dH_list, if (marginal_flag) V_beta else V_opt, marginal_flag, if (marginal_flag) active_idx else NULL, if (marginal_flag) dH_theta else NULL, if (marginal_flag) V_theta else NULL, if (marginal_flag) idx_fix_active else NULL)
  if (!is.null(df_map)) df_transform <- .apply_df_map(df_transform, df_map)


  # Execute calculation if a user-defined generate block exists
  df_generate <- build_derived_summary(self$generate, gq_list, TRUE, if (marginal_flag) dH_beta else dH_list, if (marginal_flag) V_beta else V_opt, marginal_flag, if (marginal_flag) active_idx else NULL, if (marginal_flag) dH_theta else NULL, if (marginal_flag) V_theta else NULL, if (marginal_flag) idx_fix_active else NULL)
  if (!is.null(df_map)) df_generate <- .apply_df_map(df_generate, df_map)

  log_ml <- NA
  if (!is.null(sd_rep) && !is.null(sd_rep$cov.fixed) && !fallback_needed) {
    eig <- tryCatch(eigen(sd_rep$cov.fixed, symmetric = TRUE), error = function(e) NULL)
    if (!is.null(eig) && all(eig$values > 1e-8)) {
      lj_missing <- if (laplace) calc_log_jacobian(unc_est_list, self$par_list, FALSE) - calc_log_jacobian(unc_est_list, self$par_list, TRUE) else calc_log_jacobian(unc_est_list, self$par_list, FALSE)
      log_ml <- -opt$objective + lj_missing + (length(opt$par) / 2) * log(2 * pi) + 0.5 * sum(log(eig$values)) - self$prior_correction
    }
  }



  is_se_sampling <- identical(se_method, "sampling")

  res_obj <- MAP_Fit$new(
    model = self, par_vec = con_est_vec, par = con_est_list,
    objective = opt$objective, log_ml = log_ml, convergence = opt$convergence,
    sd_rep = sd_rep, df_fixed = df_fixed, random_effects = df_random,
    df_transform = df_transform, df_generate = df_generate,
    opt_history = opt_history, transform = tran_list, generate = gq_list,
    se_samples = if (is_se_sampling && se_enabled) list(con = samps_con, tran = samps_tran, gq = samps_gq) else NULL,
    par_unc = unc_est_vec, ci_method = se_method, laplace = laplace,
    vcov_unc = Cov_u,
    map = target_map,
    marginal_vars = target_vars,
    laplace_random_vars = ad_setup$use_random,
    idx_fix_active = idx_fix_active, show_df = show_df,
    view = view,
    fallback_needed = fallback_needed
  )
  return(res_obj)
}

.resolve_optimize_df_settings <- function(df_method, marginal_flag, se_method, marginal) {
  if (is.null(df_method)) df_method <- "auto"
  df_method_l <- if (is.character(df_method)) tolower(df_method) else "numeric"

  if (identical(df_method_l, "bw")) {
    stop(
      "`df_method = 'bw'` is not supported in optimize(). ",
      "Use `df_method = 'satterthwaite'`, 'auto', 'inf', or 'none'.",
      call. = FALSE
    )
  }

  auto_df <- FALSE
  show_df <- FALSE
  df_t <- Inf

  if (identical(se_method, "none")) {
    # se_method = "none" always disables finite-df
  } else if (!marginal_flag) {
    # No marginalization means no finite-df by default in MAP
  } else if (is.numeric(df_method)) {
    df_t <- as.numeric(df_method)
    show_df <- TRUE
  } else if (df_method_l %in% c("auto", "satterthwaite")) {
    auto_df <- TRUE
    show_df <- TRUE
  } else if (df_method_l %in% c("inf", "infinite", "none")) {
    # Keep defaults (Inf, FALSE)
  } else {
    stop("Invalid df_method.", call. = FALSE)
  }

  return(list(
    auto_df = auto_df,
    show_df = show_df,
    df_t = df_t,
    df_method_norm = if (is.numeric(df_method)) "numeric" else df_method_l
  ))
}

.satterthwaite_df_delta <- function(J_full, V_theta, dH_theta, idx_theta_full, V_beta = NULL, dH_beta = NULL, idx_beta_full = NULL, max_df = NULL) {
  L_out <- nrow(J_full)
  df_out <- rep(Inf, L_out)
  rel_tol <- 1e-10

  # Optimized block (theta)
  j_theta <- J_full[, idx_theta_full, drop = FALSE]
  P_theta <- length(idx_theta_full)
  
  # Optional marginalized block (beta)
  has_beta <- !is.null(V_beta) && is.matrix(V_beta) && 
              !is.null(dH_beta) && is.list(dH_beta) &&
              !is.null(idx_beta_full) && length(idx_beta_full) > 0L

  if (has_beta) {
    j_beta <- J_full[, idx_beta_full, drop = FALSE]
    P_beta <- length(idx_beta_full)
  }

  for (i in 1:L_out) {
    ji_t <- j_theta[i, , drop = FALSE]
    v_i <- as.numeric(ji_t %*% V_theta %*% t(ji_t))
    
    grad_v <- rep(0, P_theta)
    for (k in 1:P_theta) {
      dh <- dH_theta[[k]]
      
      # Guard against NULL or dimension mismatch in numerical sensitivities
      if (is.null(dh)) next
      if (!is.matrix(dh)) next
      if (nrow(dh) != nrow(V_theta) || ncol(dh) != ncol(V_theta)) next

      # dV_theta/dtheta_k = -V_theta %*% dh %*% V_theta
      dV_k <- - V_theta %*% dh %*% V_theta
      grad_v[k] <- as.numeric(ji_t %*% dV_k %*% t(ji_t))
    }
    
    if (has_beta) {
      ji_b <- j_beta[i, , drop = FALSE]
      v_i <- v_i + as.numeric(ji_b %*% V_beta %*% t(ji_b))
      
      for (k in 1:P_theta) {
         # dV_beta/dtheta_k is directly provided in dH_beta
         dV_b_k <- dH_beta[[k]]
         
         if (is.null(dV_b_k)) next
         if (!is.matrix(dV_b_k)) next
         if (nrow(dV_b_k) != nrow(V_beta) || ncol(dV_b_k) != ncol(V_beta)) next

         grad_v[k] <- grad_v[k] + as.numeric(ji_b %*% dV_b_k %*% t(ji_b))
      }
    }

    if (all(abs(grad_v) < 1e-30)) next
    
    var_vi <- as.numeric(t(grad_v) %*% V_theta %*% grad_v)
    
    # Numerical stability: If variance of v_i is extremely small relative to v_i^2, 
    # it means sensitivity is negligible, so DF should be Infinite.
    if (!is.na(var_vi) && is.finite(var_vi) && var_vi > rel_tol * (v_i^2) && v_i > 0) {
      df_out[i] <- 2 * (v_i^2) / var_vi
    }
  }

  if (!is.null(max_df)) df_out[df_out > max_df] <- Inf
  df_out[!is.finite(df_out)] <- Inf
  df_out <- pmax(df_out, 2.1)
  
  return(df_out)
}
