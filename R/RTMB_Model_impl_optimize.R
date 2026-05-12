#' @noRd
.optimize_impl <- function(self, private, laplace, init, num_estimate, control,
                           optimizer, method, map, fixed, se, se_method, 
                           num_samples, seed, df, df_method, 
                           empirical, REML, marginal, view, 
                           is_classic, target_vars, .return_object,
                           ci_method, se_sampling, df_pars, views, ...) {
  
  if (!is.null(ci_method)) se_method <- ci_method
  if (!is.null(se_sampling)) {
    if (isTRUE(se_sampling)) se_method <- "sampling"
  }
  se_method <- match.arg(se_method, c("wald", "sampling"))
  
  if (!is.null(df_pars)) {
    if (is.null(marginal) || identical(marginal, "auto")) marginal <- df_pars
  }
  if (is.null(view) && !is.null(views)) view <- views

  reml_flag <- isTRUE(REML) || isTRUE(empirical)
  if (is.character(empirical)) {
    if (is.null(marginal) || identical(marginal, "none")) marginal <- empirical
    reml_flag <- TRUE
  }
  if (is.null(marginal)) marginal <- "none"

  dot_args <- list(...)
  requested_contrasts <- dot_args$contrasts
  if (!is.null(requested_contrasts)) {
    old_contrasts <- options()$contrasts
    on.exit(options(contrasts = old_contrasts), add = TRUE)
    if (requested_contrasts == "sum") options(contrasts = c("contr.sum", "contr.poly"))
    else if (requested_contrasts == "treatment") options(contrasts = c("contr.treatment", "contr.poly"))
    self$contrasts <- requested_contrasts
  } else if (!is.null(self$contrasts)) {
    old_contrasts <- options()$contrasts
    on.exit(options(contrasts = old_contrasts), add = TRUE)
    if (self$contrasts == "sum") options(contrasts = c("contr.sum", "contr.poly"))
    else if (self$contrasts == "treatment") options(contrasts = c("contr.treatment", "contr.poly"))
  }

  # Delegate to the corrected apply_constraints
  if (!is.null(fixed)) {
    new_model <- self$fixed_model(fixed)
    reml_to_pass <- if (is.character(empirical)) TRUE else if (!is.null(REML)) REML else isTRUE(empirical)
    
    return(new_model$optimize(
      laplace = laplace, init = new_model$init, num_estimate = num_estimate, control = control,
      optimizer = optimizer, method = method, map = map, se = se, se_method = se_method,
      num_samples = num_samples, seed = seed, df = df, df_method = df_method,
      empirical = empirical, REML = reml_to_pass, marginal = marginal, 
      is_classic = is_classic, target_vars = target_vars, .return_object = .return_object, 
      view = view, fixed = NULL, ...
    ))
  }

  if (is.null(target_vars)) {
    if (identical(marginal, "auto")) {
      if (!is.null(self$extra$marginal)) target_vars <- self$extra$marginal
      else if (!is.null(self$extra$df_pars)) target_vars <- self$extra$df_pars
      else {
        target_vars <- names(self$par_list)[sapply(names(self$par_list), function(n) {
          p_info <- self$par_list[[n]]
          is_fixed_name <- any(sapply(private$fix_patterns, function(v) grepl(paste0("^", v, "($|\\[|\\.|_)"), n)))
          is_unbounded <- identical(p_info$bounds, "none")
          (is_fixed_name || is_unbounded) && !isTRUE(p_info$random)
        })]
      }
    } else if (identical(marginal, "none")) target_vars <- character(0)
    else if (identical(marginal, "all")) target_vars <- names(self$par_list)
    else target_vars <- marginal
  }

  if (is.null(target_vars)) {
    target_vars <- character(0)
  }

  if (isTRUE(is_classic)) {
    model_obj <- self$null_model(target_vars = target_vars)
    rec_args <- dot_args; rec_args$df_method <- NULL; rec_args$view <- NULL; rec_args$views <- NULL
    return(do.call(model_obj$optimize, c(list(
      laplace = laplace, init = init, num_estimate = num_estimate, control = control,
      optimizer = optimizer, method = method, map = map, se = se, se_method = se_method,
      num_samples = num_samples, seed = seed, df = df, REML = REML, marginal = marginal,
      target_vars = target_vars, is_classic = FALSE, .return_object = "Classic", view = view, df_method = df_method
    ), rec_args)))
  }

  ci_method <- se_method
  if (ci_method == "sampling") se_sampling <- TRUE else se_sampling <- FALSE
  if (!getOption("BayesRTMB.silent", FALSE)) cat("Starting optimization...\n")

  auto_df <- FALSE; df_t <- Inf
  if (!is.null(df)) {
    if (identical(df, "auto")) auto_df <- TRUE else if (is.numeric(df)) df_t <- df
  } else {
    # --- Effective degrees of freedom calculation for Satterthwaite approximation ---
    # Calculated when the user specifies df = TRUE and REML / ML estimation is performed
    if (identical(df_method, "satterthwaite")) auto_df <- TRUE
    else if (identical(df_method, "Inf")) df_t <- Inf
    else if (is.numeric(df_method)) df_t <- df_method
    else if (identical(df_method, "bw")) {
      if (isTRUE(self$type %in% c("lmer", "glmer"))) auto_df <- TRUE else df_t <- Inf
    }
  }

  opt_results <- list(); obj_vals <- numeric(num_estimate); conv_codes <- numeric(num_estimate)

  original_par_list <- self$par_list
  if (reml_flag) {
    modified_par_list <- self$par_list
    target_map <- if (!is.null(map)) map else self$map
    for (name in names(modified_par_list)) {
      is_mapped_out <- !is.null(target_map[[name]]) && all(is.na(target_map[[name]]))
      is_fix_pattern <- any(sapply(private$fix_patterns, function(v) grepl(paste0("^", v, "($|\\[|\\.|_)"), name)))
      is_target <- name %in% target_vars
      is_var <- any(sapply(private$var_patterns, function(v) grepl(paste0("^", v, "(_[0-9]+)?$"), name)))
      is_ran_original <- isTRUE(modified_par_list[[name]]$random)
      
      if (!is_mapped_out && (is_target || (is_fix_pattern && !is_var) || is_ran_original)) {
        modified_par_list[[name]]$random <- TRUE
      } else {
        modified_par_list[[name]]$random <- FALSE
      }
    }
    self$par_list <- modified_par_list
    on.exit({ self$par_list <- original_par_list }, add = TRUE)
    laplace <- TRUE
  }

  # --- Logic for switching Jacobian adjustment ---
  # target = "all" : Necessary for finding the posterior mode (MAP)
  # target = "random" : When calculating Laplace approximation (marginal likelihood), only the random effect parameters need adjustment
  # target = "none" : Necessary for pure likelihood maximization (ML/REML)
  jac_target <- if (isTRUE(laplace)) "random" else "none"

  if (isTRUE(self$silent)) {
     ad_setup <- suppressMessages(suppressWarnings(self$build_ad_obj(init = init, laplace = laplace, jacobian_target = jac_target, map = map)))
  } else {
     ad_setup <- self$build_ad_obj(init = init, laplace = laplace, jacobian_target = jac_target, map = map)
  }
  base_ad_obj <- ad_setup$ad_obj

  for (i in 1:num_estimate) {
    if (num_estimate > 1 && !getOption("BayesRTMB.silent", FALSE)) cat(sprintf("Optimization run %d/%d...\r", i, num_estimate))
    res <- tryCatch({
      if (i > 1 && is.null(init)) base_ad_obj$par <- base_ad_obj$par + rnorm(length(base_ad_obj$par), mean = 0, sd = 0.5)
      if (length(base_ad_obj$par) == 0) {
        opt <- list(par = numeric(0), objective = base_ad_obj$fn(numeric(0)), convergence = 0, message = "all parameters fixed", iterations = 0, evaluations = c(obj = 1, grad = 0))
      } else if (optimizer == "nlminb") {
        if (is.null(control$iter.max)) control$iter.max <- 5000; if (is.null(control$eval.max)) control$eval.max <- 5000; if (is.null(control$rel.tol)) control$rel.tol <- 1e-8
        opt <- suppressWarnings(nlminb(start = base_ad_obj$par, objective = base_ad_obj$fn, gradient = base_ad_obj$gr, control = control))
      } else if (optimizer == "optim") {
        if (is.null(control$maxit)) control$maxit <- 5000
        opt <- optim(par = base_ad_obj$par, fn = base_ad_obj$fn, gr = base_ad_obj$gr, method = method, control = control); opt$objective <- opt$value
      } else stop("optimizer must be either 'optim' or 'nlminb'.")
      list(opt = opt, ad_obj = base_ad_obj)
    }, error = function(e) { warning("Optimization error: ", e$message, call. = FALSE); NULL })
    if (!is.null(res)) { opt_results[[i]] <- res; obj_vals[i] <- res$opt$objective; conv_codes[i] <- res$opt$convergence }
    else { opt_results[[i]] <- NULL; obj_vals[i] <- NA; conv_codes[i] <- NA }
  }
  if (!getOption("BayesRTMB.silent", FALSE)) cat("\n")

  valid_idx <- which(!is.na(obj_vals))
  if (length(valid_idx) == 0) stop("All optimization attempts failed.")
  best_idx <- valid_idx[which.min(obj_vals[valid_idx])]; best_res <- opt_results[[best_idx]]
  
  if (!is.na(self$prior_correction) && self$prior_correction != 0) {
    for (i in seq_along(obj_vals)) if (!is.na(obj_vals[i])) obj_vals[i] <- obj_vals[i] + self$prior_correction
    best_res$opt$objective <- best_res$opt$objective + self$prior_correction
  }

  if (!getOption("BayesRTMB.silent", FALSE)) {
    if (num_estimate == 1) {
      cat(sprintf("Optimization %s. Final objective: %.2f\n", if (conv_codes[1] == 0) "converged" else "Not Converged", obj_vals[1]))
    } else {
      cat("\nOptimization Diagnostics per estimate:\n")
      for (i in 1:num_estimate) cat(sprintf("  est%d: Objective = %10.2f, Code = %s (%s)%s\n", i, obj_vals[i], as.character(conv_codes[i]), if (conv_codes[i] == 0) "Converged" else "Not Converged", if (i == best_idx) "  <-- BEST" else ""))
      cat("\n")
    }
  }

  ad_obj <- best_res$ad_obj; opt <- best_res$opt; ad_obj$fn(opt$par)
  sd_rep <- if (se) tryCatch(RTMB::sdreport(ad_obj, getJointPrecision = reml_flag), warning = function(w) NULL, error = function(e) NULL) else NULL
  unc_est_vec <- unlist(ad_obj$env$parList(), use.names = FALSE); L_u_total <- length(unc_est_vec)
  idx_ran <- ad_obj$env$random; if (is.null(idx_ran)) idx_ran <- integer(0)
  target_map <- if (!is.null(map)) map else self$map
  
  idx_fix_active <- integer(0); idx_ran_full <- integer(0); factor_levels_seen <- list()
  opt_par_curr <- 1; idx_mapped_curr <- 1; idx_full_curr <- 1
  for (name in names(self$par_list)) {
    L_u <- self$par_list[[name]]$unc_length
    if (L_u > 0) {
      map_f <- if (!is.null(target_map[[name]])) target_map[[name]] else NULL
      for (i in 1:L_u) {
        pos_full <- idx_full_curr + i - 1
        if (!is.null(map_f) && is.na(map_f[i])) next
        if (idx_mapped_curr %in% idx_ran) { idx_ran_full <- c(idx_ran_full, pos_full); idx_mapped_curr <- idx_mapped_curr + 1; next }
        if (!is.null(map_f)) {
          lvl <- as.character(map_f[i])
          if (!(lvl %in% names(factor_levels_seen))) { factor_levels_seen[[lvl]] <- opt_par_curr; idx_fix_active <- c(idx_fix_active, pos_full); opt_par_curr <- opt_par_curr + 1 }
        } else { idx_fix_active <- c(idx_fix_active, pos_full); opt_par_curr <- opt_par_curr + 1 }
        idx_mapped_curr <- idx_mapped_curr + 1
      }
      idx_full_curr <- idx_full_curr + L_u
    }
  }

  unc_est_vec[idx_fix_active] <- opt$par
  if (laplace && !is.null(sd_rep) && length(idx_ran_full) > 0) {
    smry_ran <- tryCatch(summary(sd_rep, select = "random"), error = function(e) NULL)
    if (!is.null(smry_ran) && nrow(smry_ran) == length(idx_ran_full)) unc_est_vec[idx_ran_full] <- smry_ran[, "Estimate"]
  }

  unc_se_vec <- rep(0, L_u_total); Cov_u <- diag(0, nrow = L_u_total, ncol = L_u_total); fallback_needed <- FALSE
  
  if (reml_flag && !is.null(sd_rep) && !is.null(sd_rep$jointPrecision)) {
    V_joint <- tryCatch(
      solve(as.matrix(sd_rep$jointPrecision)), 
      error = function(e) if (requireNamespace("MASS", quietly=TRUE)) MASS::ginv(as.matrix(sd_rep$jointPrecision)) else NULL
    )
    if (!is.null(V_joint) && length(idx_ran_full) == nrow(V_joint)) {
      Cov_u[idx_ran_full, idx_ran_full] <- as.matrix(V_joint)
      unc_se_vec[idx_ran_full] <- sqrt(pmax(diag(Cov_u[idx_ran_full, idx_ran_full, drop = FALSE]), 0))
    }
  }

  if (!is.null(sd_rep) && !is.null(sd_rep$cov.fixed)) {
    se_fix <- sqrt(pmax(diag(sd_rep$cov.fixed), 0))
    if (!isTRUE(sd_rep$pdHess) || any(is.na(se_fix)) || any(is.nan(se_fix)) || length(idx_fix_active) != length(se_fix)) fallback_needed <- TRUE
    else {
      unc_se_vec[idx_fix_active] <- se_fix; Cov_u[idx_fix_active, idx_fix_active] <- sd_rep$cov.fixed
      idx_curr <- 1
      for (name in names(self$par_list)) {
        L_u <- self$par_list[[name]]$unc_length
        if (L_u > 0) {
          map_f <- if (!is.null(target_map[[name]])) target_map[[name]] else NULL
          for (i in 1:L_u) { 
            pos_last <- idx_curr + i - 1; 
            if (pos_last %in% idx_ran) next; 
            if (!is.null(map_f) && !is.na(map_f[i])) { 
              lvl <- as.character(map_f[i]); 
              mapped_opt_idx <- factor_levels_seen[[lvl]]; 
              if (!is.null(mapped_opt_idx)) unc_se_vec[pos_last] <- se_fix[mapped_opt_idx] 
            } 
          }
          idx_curr <- idx_curr + L_u
        }
      }
      # 1. Extract values from the sdreport object
      # 2. Calculate degrees of freedom, t-values, and p-values (if df = TRUE)
      # 3. Construct confidence intervals
      # 4. Restore original dimensions (dim)
      if (laplace && length(idx_ran_full) > 0) { smry_ran <- tryCatch(summary(sd_rep, select = "random"), error = function(e) NULL); if (!is.null(smry_ran)) unc_se_vec[idx_ran_full] <- smry_ran[, "Std. Error"] }
    }
  } else fallback_needed <- TRUE

  if (se && fallback_needed) {
    H <- tryCatch(ad_obj$he(opt$par + 1e-6), error = function(e) NULL)
    if (!is.null(H) && requireNamespace("MASS", quietly = TRUE)) {
      Cov_pseudo <- tryCatch(MASS::ginv(H), error = function(e) NULL)
      if (!is.null(Cov_pseudo)) {
        se_pseudo <- sqrt(pmax(diag(Cov_pseudo), 0))
        if (length(idx_fix_active) == length(se_pseudo)) {
          unc_se_vec[idx_fix_active] <- se_pseudo; Cov_u[idx_fix_active, idx_fix_active] <- Cov_pseudo
          idx_curr <- 1
          for (name in names(self$par_list)) {
            L_u <- self$par_list[[name]]$unc_length
            if (L_u > 0) {
              map_f <- if (!is.null(target_map[[name]])) target_map[[name]] else NULL
              for (i in 1:L_u) { pos_last <- idx_curr + i - 1; if (pos_last %in% idx_ran) next; if (!is.null(map_f) && !is.na(map_f[i])) { lvl <- as.character(map_f[i]); mapped_opt_idx <- factor_levels_seen[[lvl]]; if (!is.null(mapped_opt_idx)) unc_se_vec[pos_last] <- se_pseudo[mapped_opt_idx] } }
              idx_curr <- idx_curr + L_u
            }
          }
        }
      }
    }
  }

  unc_est_list <- unconstrained_vector_to_list(unc_est_vec, self$par_list); unc_se_list <- unconstrained_vector_to_list(unc_se_vec, self$par_list)
  con_est_list <- self$to_constrained(unc_est_list); diag(Cov_u) <- pmax(diag(Cov_u), unc_se_vec^2)
  con_se_list <- list(); con_lower_list <- list(); con_upper_list <- list(); samps_con <- list(); samps_tran <- list(); samps_gq <- list()

  dH_list <- NULL; V_opt <- NULL; active_idx <- NULL; dH_theta <- dH_beta <- V_theta <- V_beta <- NULL
  if (any(!is.infinite(df_t))) est_dfs_all <- rep(df_t[1], L_u_total)
  else if (auto_df || (reml_flag && (!identical(df_method, "bw") || isTRUE(self$type %in% c("lmer", "glmer"))))) {
    is_silent <- identical(.return_object, "Classic") && identical(df_method, "bw")
    satt_res <- self$calculate_satterthwaite_df(ad_obj, idx_fix_active, L_u_total, opt$par, silent = is_silent, return_sensitivities = TRUE)
    if (is.list(satt_res)) { est_dfs_all <- satt_res$df; dH_list <- satt_res$sensitivities; V_opt <- satt_res$V } else est_dfs_all <- satt_res
    dH_theta <- dH_list; V_theta <- V_opt
    if (reml_flag && length(ad_obj$env$random) > 0) {
      full_par_names_mapped <- names(ad_obj$env$last.par)[idx_ran]; target_ran_idx <- which(full_par_names_mapped %in% target_vars)
      if (length(target_ran_idx) > 0) {
        reml_res <- self$calculate_reml_satterthwaite_df(ad_obj, opt$par, target_ran_idx, silent = is_silent, return_sensitivities = TRUE)
        if (is.list(reml_res)) { active_idx <- idx_ran_full[target_ran_idx]; est_dfs_all[active_idx] <- reml_res$df; dH_beta <- reml_res$sensitivities; V_beta <- reml_res$V_beta }
      }
    }
    ad_obj$fn(opt$par)
  } else if (identical(df_method, "bw")) {
    df_list_unc <- unconstrained_vector_to_list(rep(Inf, L_u_total), self$par_list); type <- self$type
    
    raw_df <- if (!is.null(self$raw_data)) as.data.frame(self$raw_data) else NULL
    N_obs <- 0
    if (!is.null(raw_df)) {
      N_obs <- nrow(raw_df)
    } else if (!is.null(self$data)) {
      if (!is.null(self$data$Y)) N_obs <- length(self$data$Y)
      else if (!is.null(self$data$X) && is.matrix(self$data$X)) N_obs <- nrow(self$data$X)
      else if (!is.null(self$data$Y1) && !is.null(self$data$Y2)) N_obs <- length(self$data$Y1) + length(self$data$Y2)
      else if (!is.null(self$data$diffs)) N_obs <- length(self$data$diffs)
      else {
        for (item in self$data) {
          if (is.numeric(item) && is.vector(item) && length(item) > 1) { N_obs <- length(item); break }
        }
      }
    }
    if (is.null(N_obs) || N_obs == 0) N_obs <- 0

    if (type == "mediation" && !is.null(self$extra$df_map)) {
      df_map <- self$extra$df_map

      for (v in names(df_list_unc)) {
        v_display <- sub("^b_c", "b", v)

        if (v %in% names(df_map)) {
          df_list_unc[[v]] <- rep(as.numeric(df_map[[v]]), length(df_list_unc[[v]]))
        } else if (v_display %in% names(df_map)) {
          df_list_unc[[v]] <- rep(as.numeric(df_map[[v_display]]), length(df_list_unc[[v]]))
        }
      }

    } else if (type %in% c("lm", "ttest", "mediation")) {
      K_fixed <- sum(sapply(target_vars, function(v) if (v %in% names(self$par_list)) self$par_list[[v]]$unc_length else 0))
      for (v in target_vars) if (v %in% names(df_list_unc)) df_list_unc[[v]] <- rep(if (v == "delta") Inf else pmax(N_obs - K_fixed, 1), length(df_list_unc[[v]]))
    } else if (type == "corr") {
      if ("corr" %in% names(df_list_unc)) df_list_unc[["corr"]] <- rep(N_obs - 2, length(df_list_unc[["corr"]]))
      for (v in c("mean", "sd")) if (v %in% names(df_list_unc)) df_list_unc[[v]] <- rep(N_obs - 1, length(df_list_unc[[v]]))
    }
    est_dfs_all <- unlist(df_list_unc, use.names = FALSE)
  } else est_dfs_all <- rep(Inf, L_u_total)

  if (se_sampling) {
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
        for (i in 1:L_u) { pos_last <- idx_curr + i - 1; if (pos_last %in% idx_ran) next; if (!is.null(map_f) && !is.na(map_f[i])) { lvl <- as.character(map_f[i]); mapped_opt_idx <- factor_levels_seen[[lvl]]; if (!is.null(mapped_opt_idx)) u_samples_mat[, pos_last] <- raw_samples[, mapped_opt_idx] } }
        idx_curr <- idx_curr + L_u
      }
    }
    if (laplace && length(idx_ran) > 0) for (ridx in idx_ran) u_samples_mat[, ridx] <- rnorm(num_samples, mean = unc_est_vec[ridx], sd = unc_se_vec[ridx])
    
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
  } else {
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
  }

  if (reml_flag) self$par_list <- original_par_list

  show_df <- auto_df || reml_flag || !is.infinite(df_t)
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
    if (show_df || !se_sampling) {
      J <- matrix(0, L_out, L_u_total); for (i in 1:L_u_total) { u_tmp <- unc_est_vec; u_tmp[i] <- u_tmp[i] + 1e-5; tmp_con <- self$to_constrained(unconstrained_vector_to_list(u_tmp, self$par_list)); if (is_generate && !is.null(self$transform)) { ut <- tryCatch(self$transform(self$data, tmp_con), error = function(e) NULL); if (!is.null(ut)) tmp_con <- c(tmp_con, ut) }; tmp_out <- tryCatch(func(self$data, tmp_con), error = function(e) NULL); if (!is.null(tmp_out)) J[, i] <- (unlist(tmp_out, use.names = FALSE) - flat_base) / 1e-5 }
    }
    if (se_sampling) {
      samps_list <- if (is_generate) samps_gq else samps_tran
      if (length(samps_list) == 0) { se_out <- low_out <- up_out <- rep(NA, L_out) }
      else { mat_all <- do.call(cbind, unname(samps_list)); se_out <- apply(mat_all, 2, sd, na.rm = TRUE); low_out <- apply(mat_all, 2, quantile, 0.025, na.rm = TRUE); up_out <- apply(mat_all, 2, quantile, 0.975, na.rm = TRUE) }
    } else {
      # Fallback when the inverse Hessian is not positive definite
      # Either estimate SE using only diagonal components or return NA
      se_out <- rep(NA_real_, L_out); smry_rep <- if (!is.null(sd_rep)) tryCatch(as.data.frame(summary(sd_rep, select = "report")), error = function(e) NULL) else NULL
      if (!is.null(smry_rep)) { m_idx <- match(names_vec, rownames(smry_rep)); valid <- !is.na(m_idx); if (any(valid)) se_out[valid] <- smry_rep$`Std. Error`[m_idx[valid]] }
      for (j in which(is.na(se_out))) se_out[j] <- sqrt(abs(sum((J[j, ] %*% Cov_u) * J[j, ])))
      low_out <- flat_base - 1.96 * se_out; up_out <- flat_base + 1.96 * se_out
      corr_pat <- "^(corr|pcorr|B_corr|W_corr)(\\[|$)"; is_corr <- grepl(corr_pat, names_vec)
      for (j in which(is_corr)) { r <- pmax(pmin(flat_base[j], 0.9999), -0.9999); z <- atanh(r); se_z <- se_out[j] / (1 - r^2 + 1e-10); low_out[j] <- tanh(z - 1.96 * se_z); up_out[j] <- tanh(z + 1.96 * se_z) }
    }
    res_df <- data.frame(Estimate = flat_base, `Std. Error` = se_out, `Lower 95%` = low_out, `Upper 95%` = up_out, check.names = FALSE)
    if (show_df) {
      derived_dfs <- rep(Inf, L_out)
      if (!is.null(dH_list_in) && !is.null(V_opt_in)) {
        if (!is_reml_in && length(dH_list_in) == ncol(J)) {
          A <- J %*% V_opt_in; for (j in 1:L_out) { Aj <- A[j, , drop = FALSE]; grad_v <- sapply(dH_list_in, function(dh) -as.numeric(Aj %*% dh %*% t(Aj))); var_vj <- as.numeric(t(grad_v) %*% V_opt_in %*% grad_v); vj <- as.numeric(Aj %*% J[j, ]); if (var_vj > 1e-30 && vj > 0) derived_dfs[j] <- 2 * (vj^2) / var_vj }
        } else if (is_reml_in) {
          for (j in 1:L_out) {
            var_vj <- 0
            vj <- se_out[j]^2

            if (!is.null(active_idx_in) &&
                !is.null(dH_list_in) &&
                !is.null(V_theta_in) &&
                length(dH_list_in) > 0) {

              Ja <- J[j, active_idx_in, drop = FALSE]
              first_dh <- dH_list_in[[which(!vapply(dH_list_in, is.null, logical(1)))[1]]]

              if (!is.null(first_dh) && ncol(Ja) == nrow(first_dh)) {
                grad_beta <- numeric(length(dH_list_in))

                for (k in seq_along(dH_list_in)) {
                  dh <- dH_list_in[[k]]
                  if (!is.null(dh) && ncol(Ja) == nrow(dh)) {
                    grad_beta[k] <- as.numeric(Ja %*% dh %*% t(Ja))
                  }
                }

                if (length(grad_beta) == nrow(V_theta_in)) {
                  var_vj <- var_vj + as.numeric(t(grad_beta) %*% V_theta_in %*% grad_beta)
                }
              }
            }

            if (!is.null(dH_theta_in) &&
                !is.null(V_theta_in) &&
                length(dH_theta_in) > 0) {

              theta_idx <- theta_idx_in
              if (is.null(theta_idx)) {
                theta_idx <- seq_len(nrow(V_theta_in))
              }

              Jt <- J[j, theta_idx, drop = FALSE]
              first_dh <- dH_theta_in[[which(!vapply(dH_theta_in, is.null, logical(1)))[1]]]

              if (!is.null(first_dh) && ncol(Jt) == nrow(first_dh)) {
                At <- Jt %*% V_theta_in
                grad_theta <- numeric(length(dH_theta_in))

                for (k in seq_along(dH_theta_in)) {
                  dh <- dH_theta_in[[k]]
                  if (!is.null(dh) && ncol(At) == nrow(dh)) {
                    grad_theta[k] <- -as.numeric(At %*% dh %*% t(At))
                  }
                }

                if (length(grad_theta) == nrow(V_theta_in)) {
                  var_vj <- var_vj + as.numeric(t(grad_theta) %*% V_theta_in %*% grad_theta)
                }
              }
            }

            if (!is.na(var_vj) && var_vj > 1e-30 && vj > 0) {
              derived_dfs[j] <- 2 * (vj^2) / var_vj
            }
          }
        }
      } else if (identical(df_method, "bw")) {
        for (j in 1:L_out) {
          if (self$type == "corr") { if (grepl(corr_pat, names_vec[j])) derived_dfs[j] <- N_obs - 2 else derived_dfs[j] <- N_obs - 1 }
          else { c_dfs <- est_dfs_all[abs(J[j, ]) > 1e-8]; if (length(c_dfs) > 0) derived_dfs[j] <- min(c_dfs) }
        }
      }
      
      if (!is.null(df_map)) {
        res_df$df <- derived_dfs
        res_df <- .apply_df_map(res_df, df_map)
      } else {
        res_df$df <- derived_dfs
      }
    }
    rownames(res_df) <- names_vec; return(res_df)
  }

  df_transform <- build_derived_summary(self$transform, tran_list, FALSE, if (reml_flag) dH_beta else dH_list, if (reml_flag) V_beta else V_opt, reml_flag, if (reml_flag) active_idx else NULL, if (reml_flag) dH_theta else NULL, if (reml_flag) V_theta else NULL, if (reml_flag) idx_fix_active else NULL)
  if (!is.null(df_map)) df_transform <- .apply_df_map(df_transform, df_map)
  
  if (identical(df_method, "bw") && self$type == "ttest") {
    if (!is.null(df_transform) && "df" %in% names(df_transform)) {
      K_val <- sum(sapply(target_vars, function(v) if (v %in% names(self$par_list)) self$par_list[[v]]$unc_length else 0))
      for (i in 1:nrow(df_transform)) {
        nm <- rownames(df_transform)[i]
        if (grepl("^delta", nm)) df_transform$df[i] <- Inf
        else df_transform$df[i] <- pmax(N_obs - K_val, 1)
      }
    }
  }

  # Execute calculation if a user-defined generate block exists
  df_generate <- build_derived_summary(self$generate, gq_list, TRUE, if (reml_flag) dH_beta else dH_list, if (reml_flag) V_beta else V_opt, reml_flag, if (reml_flag) active_idx else NULL, if (reml_flag) dH_theta else NULL, if (reml_flag) V_theta else NULL, if (reml_flag) idx_fix_active else NULL)
  if (!is.null(df_map)) df_generate <- .apply_df_map(df_generate, df_map)

  log_ml <- NA
  if (!is.null(sd_rep) && !is.null(sd_rep$cov.fixed) && !fallback_needed) {
    eig <- tryCatch(eigen(sd_rep$cov.fixed, symmetric = TRUE), error = function(e) NULL)
    if (!is.null(eig) && all(eig$values > 1e-8)) {
      lj_missing <- if (laplace) calc_log_jacobian(unc_est_list, self$par_list, FALSE) - calc_log_jacobian(unc_est_list, self$par_list, TRUE) else calc_log_jacobian(unc_est_list, self$par_list, FALSE)
      log_ml <- -opt$objective + lj_missing + (length(opt$par) / 2) * log(2 * pi) + 0.5 * sum(log(eig$values)) - self$prior_correction
    }
  }

  if (.return_object == "Classic") {
    V_fixed_full <- NULL; test_results <- list()
    if (self$type == "table") { tab <- self$extra$tab; test_results$chisq <- chisq.test(tab); test_results$fisher <- try(fisher.test(tab), silent = TRUE) }
    
    if (reml_flag && !is.null(sd_rep$jointPrecision)) {
      V_joint <- tryCatch(solve(as.matrix(sd_rep$jointPrecision)), error = function(e) MASS::ginv(as.matrix(sd_rep$jointPrecision)))
      V_fixed_full <- matrix(0, nrow(df_fixed), nrow(df_fixed), dimnames = list(rownames(df_fixed), rownames(df_fixed)))
      idx_joint <- unlist(lapply(target_vars, function(nm) which(names(sd_rep$par.random) == nm)))
      idx_df_fe <- unlist(lapply(target_vars, function(nm) grep(paste0("^", nm, "(\\[|$)"), rownames(df_fixed))))
      if (length(idx_joint) > 0 && length(idx_joint) == length(idx_df_fe)) V_fixed_full[idx_df_fe, idx_df_fe] <- V_joint[idx_joint, idx_joint]
      idx_df_var <- setdiff(1:nrow(df_fixed), idx_df_fe); if (length(idx_df_var) > 0 && !is.null(sd_rep$cov.fixed)) V_fixed_full[idx_df_var, idx_df_var] <- sd_rep$cov.fixed
    }
    
    if (is.null(V_fixed_full)) {
      L_fixed <- nrow(df_fixed)
      J_fixed <- matrix(0, L_fixed, L_u_total)
      base_f <- df_fixed$Estimate
      for (i in 1:L_u_total) {
        u_tmp <- unc_est_vec; u_tmp[i] <- u_tmp[i] + 1e-5
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
      rownames(V_fixed_full) <- colnames(V_fixed_full) <- rownames(df_fixed)
    }
    
    cf_args <- dot_args; cf_args$df_method <- NULL; cf_args$contrasts <- NULL
    res_obj <- do.call(Classic_Fit$new, c(list(
      model = self, par_vec = con_est_vec, par = con_est_list, objective = opt$objective, 
      log_ml = log_ml, convergence = opt$convergence, sd_rep = sd_rep, 
      df_fixed = df_fixed, random_effects = df_random, df_transform = df_transform, 
      df_generate = df_generate, opt_history = data.frame(estimate = 1:num_estimate, objective = obj_vals, code = conv_codes), 
      transform = tran_list, generate = gq_list, se_samples = if (se_sampling) list(con = samps_con, tran = samps_tran, gq = samps_gq) else NULL, 
      par_unc = unc_est_vec, vcov_unc = Cov_u, ci_method = ci_method, laplace = laplace, map = target_map, 
      test_results = test_results, view = view, vcov = V_fixed_full, df_method = df_method), cf_args))
  } else {
    res_obj <- MAP_Fit$new(model = self, par_vec = con_est_vec, par = con_est_list, objective = opt$objective, log_ml = log_ml, convergence = opt$convergence, sd_rep = sd_rep, df_fixed = df_fixed, random_effects = df_random, df_transform = df_transform, df_generate = df_generate, opt_history = data.frame(estimate = 1:num_estimate, objective = obj_vals, code = conv_codes), transform = tran_list, generate = gq_list, se_samples = if (se_sampling) list(con = samps_con, tran = samps_tran, gq = samps_gq) else NULL, par_unc = unc_est_vec, ci_method = ci_method, laplace = laplace, map = target_map, vcov_unc = Cov_u, marginal_vars = ad_setup$use_random)
  }
  return(res_obj)
}

#' @noRd
.classic_impl <- function(self, private, df, df_method, REML, marginal, 
                          view, map, fixed, df_pars, views, ...) {
  if (is.null(df_method)) {
     df_method <- if (isTRUE(self$type == "ttest") || isTRUE(self$type %in% c("lmer", "glmer"))) "satterthwaite" else "bw"
  }
  if (is.null(marginal) && !is.null(df_pars)) marginal <- df_pars
  if (is.null(view) && !is.null(views)) view <- views

  return(self$optimize(
    is_classic = TRUE, 
    df = df, 
    df_method = df_method,
    # --- Hessian correction for REML (when REML=TRUE) ---
    REML = REML, 
    marginal = marginal, 
    view = view, 
    map = map, 
    fixed = fixed, 
    .return_object = "Classic", 
    ...
  ))
}
