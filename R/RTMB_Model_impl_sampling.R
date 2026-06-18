#' @noRd
.sample_impl <- function(self, private, sampling, warmup, chains, thin, seed, delta, 
                         max_treedepth, nuts_variant, metric,
                         metric_init,
                         metric_adaptation,
                         metric_regularization,
                         metric_shrinkage, metric_min, metric_max,
                         parallel, laplace, init, init_jitter,
                         save_csv, map, fixed, globals, progress) {
  nuts_variant <- match.arg(nuts_variant, c("multinomial", "slice"))
  metric <- match.arg(metric, c("auto", "diag", "dense", "hybrid"))
  metric_init <- match.arg(metric_init, c("identity", "hessian"))
  metric_adaptation <- match.arg(metric_adaptation, c("stan_window", "cumulative", "window"))
  progress_mode <- .rtmb_resolve_progress(progress)
  if (!is.null(fixed)) {
    return(private$.dispatch_fixed(.method_to_call = "sample", sampling = sampling, warmup = warmup, chains = chains, thin = thin,
      seed = seed, delta = delta, max_treedepth = max_treedepth,
      nuts_variant = nuts_variant,
      metric = metric,
      metric_init = metric_init,
      metric_adaptation = metric_adaptation,
      metric_regularization = metric_regularization,
      metric_shrinkage = metric_shrinkage,
      metric_min = metric_min,
      metric_max = metric_max,
      parallel = parallel, laplace = laplace, init = init,
      init_jitter = init_jitter, save_csv = save_csv, map = map, fixed = fixed,
      globals = globals, progress = progress))
  }

  set.seed(seed)
  orig_pl <- self$par_list
  data_na_summary <- function(dat) {
    if (is.null(dat)) return(character(0))
    out <- character(0)
    for (nm in names(dat)) {
      x <- dat[[nm]]
      if (is.atomic(x) || is.matrix(x) || is.array(x) || is.data.frame(x)) {
        n_na <- sum(is.na(x))
        if (n_na > 0L) out <- c(out, sprintf("%s (%d NA)", nm, n_na))
      }
    }
    out
  }
  stop_nonfinite_lp <- function(context) {
    na_vars <- data_na_summary(local_data)
    if (length(na_vars) > 0L) {
      stop(
        context,
        " produced non-finite log-probability values. The model data contain missing values in: ",
        paste(na_vars, collapse = ", "),
        ". If these NA values are not handled explicitly in the model code, remove/impute them or use a wrapper/missing-data option that supports NA.",
        call. = FALSE
      )
    }
    stop(context, " produced non-finite log-probability values.", call. = FALSE)
  }
  mcmc_pd_error_to_neginf <- isTRUE(getOption("BayesRTMB.mcmc_pd_error_to_neginf", TRUE))
  is_positive_definite_error <- function(msg) {
    grepl(
      paste(
        c(
          "not positive definite",
          "positive[- ]definite",
          "leading minor",
          "computationally singular",
          "reciprocal condition number",
          "\\bsingular\\b",
          "Cholesky",
          "\\bchol\\b",
          "\\bLLT\\b",
          "\\bLDLT\\b"
        ),
        collapse = "|"
      ),
      msg,
      ignore.case = TRUE
    )
  }
  wrap_mcmc_pd_errors <- function(ad_obj) {
    if (!mcmc_pd_error_to_neginf) return(ad_obj)
    pd_error_state <- new.env(parent = emptyenv())
    pd_error_state$post_warmup <- FALSE
    pd_error_state$count <- 0L
    ad_obj$rtmb_pd_error_state <- pd_error_state
    orig_fn <- ad_obj$fn
    orig_gr <- ad_obj$gr
    ad_obj$fn <- function(x, ...) {
      tryCatch(
        orig_fn(x, ...),
        error = function(e) {
          if (is_positive_definite_error(conditionMessage(e))) {
            if (isTRUE(pd_error_state$post_warmup)) pd_error_state$count <- pd_error_state$count + 1L
            return(Inf)
          }
          stop(e)
        }
      )
    }
    ad_obj$gr <- function(x, ...) {
      n <- length(x)
      tryCatch(
        orig_gr(x, ...),
        error = function(e) {
          if (is_positive_definite_error(conditionMessage(e))) {
            if (isTRUE(pd_error_state$post_warmup)) pd_error_state$count <- pd_error_state$count + 1L
            return(rep(0, n))
          }
          stop(e)
        }
      )
    }
    ad_obj
  }

  if (!is.null(save_csv)) {
    if (!is.list(save_csv)) stop("save_csv must be specified in the format list(name='...', dir='...').")
    save_name <- if (!is.null(save_csv$name)) save_csv$name else "model"
    save_dir <- if (!is.null(save_csv$dir)) save_csv$dir else "BayesRTMB_mcmc"
    save_freq <- if (!is.null(save_csv$freq)) save_csv$freq else 0
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    save_info <- list(name = save_name, dir = save_dir, freq = save_freq)
  } else {
    save_info <- NULL
  }

  random_flags <- vapply(orig_pl, function(x) isTRUE(x$random), logical(1))

  if (laplace && any(random_flags)) {
    pl_fixed <- parse_parameters(orig_pl[!random_flags], self$par_names)
    pl_random <- parse_parameters(orig_pl[random_flags], self$par_names)
    fixed_idx <- which(self$pl_full$names %in% pl_fixed$names)
    random_idx <- which(self$pl_full$names %in% pl_random$names)
  } else {
    pl_fixed <- self$pl_full
    pl_random <- NULL
    fixed_idx <- 1:length(self$pl_full$names)
    random_idx <- integer(0)
  }

  P_fixed <- length(pl_fixed$names)
  P_random <- if (!is.null(pl_random)) length(pl_random$names) else 0

  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop(
        "parallel = TRUE for sample() requires the suggested package 'future'. ",
        "Please install it or use parallel = FALSE.",
        call. = FALSE
      )
    }
    if (inherits(future::plan(), "sequential")) {
      if (.Platform$OS.type == "unix") future::plan(future::multicore, workers = chains)
      else future::plan(future::multisession, workers = chains)
    }
    .rtmb_progress_start_line(paste0("Starting parallel sampling (chains = ", chains, ")..."))
  } else {
    .rtmb_progress_start_line(paste0("Starting sequential sampling (chains = ", chains, ")..."))
  }

  # --- [IMPORTANT] Data extraction to avoid serialization ---
  local_data <- self$data
  local_par_list <- self$par_list
  local_pl_full <- self$pl_full
  local_map <- if (is.null(map)) self$map else utils::modifyList(as.list(self$map), as.list(map))
  # Normalize map factor levels for consistent MakeADFun behavior
  if (!is.null(local_map)) {
    local_map <- lapply(local_map, function(m) {
      if (is.null(m)) return(NULL)
      vals <- as.character(m)
      factor(vals, levels = unique(vals[!is.na(vals)]))
    })
  }
  local_log_prob <- self$log_prob
  local_transform <- self$transform
  local_fixed_prior_specs <- self$fixed_prior_specs
  local_code_model <- self$code$model
  
  # Pre-calculate shadow list for fixed-prior removal
  env_shadow_local <- list2env(local_data, parent = asNamespace("BayesRTMB"))
  idx_s <- 1
  for (nm in names(local_par_list)) {
    p <- local_par_list[[nm]]
    shadow_vec <- local_pl_full$names[idx_s:(idx_s + p$length - 1)]
    if (length(p$dim) > 1) dim(shadow_vec) <- p$dim
    env_shadow_local[[nm]] <- shadow_vec
    idx_s <- idx_s + p$length
  }
  if (!is.null(local_transform)) {
    par_only <- mget(names(local_par_list), envir = env_shadow_local, inherits = FALSE)
    tran_res_shadow <- tryCatch(suppressMessages(local_transform(local_data, par_only)), error = function(e) list())
    for (nm in names(tran_res_shadow)) env_shadow_local[[nm]] <- tran_res_shadow[[nm]]
  }
  local_shadow_list <- as.list(env_shadow_local)

  random_effs <- names(local_par_list)[random_flags]
  use_random <- if (laplace && length(random_effs) > 0) random_effs else NULL
  base_init <- self$prepare_init(init)

  # [Fix] Prepare f_ad beforehand (to avoid capturing private)
  f_ad_global <- private$.build_f_ad(
    data_local = local_data,
    par_list_local = local_par_list,
    log_prob_local = local_log_prob,
    transform_local = local_transform,
    jacobian_target = "all",
    adreport = FALSE,
    fixed_prior_specs_local = local_fixed_prior_specs,
    code_model_local = local_code_model,
    shadow_list_local = local_shadow_list
  )
  environment(f_ad_global) <- list2env(list(
    data_local = local_data,
    par_list_local = local_par_list,
    log_prob_local = local_log_prob,
    transform_local = local_transform,
    jacobian_target = "all",
    adreport = FALSE,
    fixed_prior_specs_local = local_fixed_prior_specs,
    code_model_local = local_code_model,
    shadow_list_local = local_shadow_list
  ), parent = asNamespace("BayesRTMB"))

  run_chain <- function(c, f_ad, p_callback = NULL) {
    unc_init_list <- to_unconstrained(constrained_vector_to_list(base_init, local_par_list), local_par_list)
    unc_init_vec <- unlist(unc_init_list, use.names = FALSE)

    if (init_jitter > 0) {
      jitter_vec <- rnorm(length(unc_init_vec), mean = 0, sd = init_jitter)
      if (!is.null(local_map)) {
        idx <- 1
        for (name in names(local_par_list)) {
          L_unc <- local_par_list[[name]]$unc_length
          if (L_unc > 0) {
            if (!is.null(local_map[[name]])) {
              is_fixed <- is.na(local_map[[name]])
              jitter_vec[idx:(idx + L_unc - 1)][is_fixed] <- 0
            }
            idx <- idx + L_unc
          }
        }
      }
      unc_init_vec <- unc_init_vec + jitter_vec
    }
    unc_init_list_new <- unconstrained_vector_to_list(unc_init_vec, local_par_list)

    # Use the f_ad provided from outside
    ad_obj <- tryCatch({
      RTMB::MakeADFun(func = f_ad, parameters = unc_init_list_new, random = use_random, map = local_map, silent = TRUE)
    }, error = function(e) stop(.rtmb_format_makeadfun_error(e$message, context = "MakeADFun in parallel worker"), call. = FALSE))
    ad_obj <- wrap_mcmc_pd_errors(ad_obj)

    metric_random_idx <- integer(0)
    if (metric %in% c("auto", "hybrid") && !isTRUE(laplace)) {
      active_is_random <- rep(FALSE, length(ad_obj$par))
      active_pos <- 1L
      for (name in names(local_par_list)) {
        p <- local_par_list[[name]]
        L_unc <- p$unc_length
        if (L_unc <= 0L) next
        map_i <- local_map[[name]]
        n_active <- L_unc
        if (!is.null(map_i)) {
          non_na <- !is.na(map_i)
          n_active <- length(unique(as.character(map_i[non_na])))
        }
        if (n_active > 0L) {
          idx <- active_pos:(active_pos + n_active - 1L)
          idx <- idx[idx <= length(active_is_random)]
          if (isTRUE(p$random)) active_is_random[idx] <- TRUE
          active_pos <- active_pos + n_active
        }
      }
      if ((active_pos - 1L) == length(ad_obj$par)) {
        metric_random_idx <- which(active_is_random)
      }
    }

    init_lp <- tryCatch(-ad_obj$fn(ad_obj$par), error = function(e) NA_real_)
    if (!is.finite(init_lp)) {
      stop_nonfinite_lp(sprintf("MCMC initialization for chain %d", c))
    }

    res <- NUTS_method(model = ad_obj, sampling = sampling, warmup = warmup, delta = delta, max_treedepth = max_treedepth, chain = c, update_progress = p_callback, laplace = laplace, save_info = save_info, nuts_variant = nuts_variant, metric = metric, metric_random_idx = metric_random_idx, metric_init = metric_init, metric_adaptation = metric_adaptation, metric_regularization = metric_regularization, metric_shrinkage = metric_shrinkage, metric_min = metric_min, metric_max = metric_max)
    if (is.null(res$lp) || all(!is.finite(res$lp))) {
      stop_nonfinite_lp(sprintf("MCMC sampling for chain %d", c))
    }

    P_all_true <- length(local_pl_full$names)
    retained_index <- seq(from = (warmup + 1), to = (warmup + sampling), by = thin)
    para_final <- array(NA, dim = c(length(retained_index), P_all_true))

    for (ii in seq_along(retained_index)) {
      i <- retained_index[ii]
      x_in <- as.numeric(res$para_fixed[i, ])
      if (laplace && length(ad_obj$env$random) > 0) {
        ad_obj$fn(x_in)
        para_list <- ad_obj$env$parList()
      } else {
        para_list <- ad_obj$env$parList(x = x_in)
      }
      con_list <- to_constrained(para_list, local_par_list)
      para_final[ii, ] <- unlist(con_list, use.names = FALSE)
    }
    res$para <- para_final
    res$retained_index <- retained_index
    res$para_full <- NULL
    if (!is.null(p_callback)) {
      p_callback(msg = paste0("chain ", c, " done (100%)"), amt = 1)
    }
    return(res)
  }

  worker_env <- list2env(list(
    base_init = base_init,
    local_par_list = local_par_list,
    init_jitter = init_jitter,
    local_map = local_map,
    use_random = use_random,
    metric = metric,
    laplace = laplace,
    sampling = sampling,
    warmup = warmup,
    delta = delta,
    max_treedepth = max_treedepth,
    save_info = save_info,
    nuts_variant = nuts_variant,
    metric_init = metric_init,
    metric_adaptation = metric_adaptation,
    metric_regularization = metric_regularization,
    metric_shrinkage = metric_shrinkage,
    metric_min = metric_min,
    metric_max = metric_max,
    thin = thin,
    local_pl_full = local_pl_full,
    local_data = local_data,
    mcmc_pd_error_to_neginf = mcmc_pd_error_to_neginf
  ), parent = asNamespace("BayesRTMB"))
  environment(data_na_summary) <- worker_env
  worker_env$data_na_summary <- data_na_summary
  environment(stop_nonfinite_lp) <- worker_env
  worker_env$stop_nonfinite_lp <- stop_nonfinite_lp
  environment(is_positive_definite_error) <- worker_env
  worker_env$is_positive_definite_error <- is_positive_definite_error
  environment(wrap_mcmc_pd_errors) <- worker_env
  worker_env$wrap_mcmc_pd_errors <- wrap_mcmc_pd_errors
  environment(run_chain) <- worker_env

  future_globals <- function(extra = list()) {
    if (isTRUE(globals)) return(TRUE)
    c(list(run_chain = run_chain, f_ad_global = f_ad_global), extra)
  }

  results_list <- list()
  if (parallel) {
    iter <- sampling + warmup
    total_updates <- chains * (2 + floor(iter / 200))
    if (identical(progress_mode, "none")) {
      results_list <- withCallingHandlers({
        futures <- lapply(seq_len(chains), function(c) {
          chain_id <- c
          future::future({
            run_chain(chain_id, f_ad = f_ad_global, p_callback = function(...) invisible(NULL))
          }, seed = TRUE, packages = "RTMB", globals = future_globals(list(chain_id = chain_id)))
        })
        lapply(futures, future::value)
      }, warning = function(w) { if (grepl("package:BayesRTMB", conditionMessage(w))) invokeRestart("muffleWarning") })
    } else {
      progress_dir <- .rtmb_progress_file_dir()
      on.exit(unlink(progress_dir, recursive = TRUE, force = TRUE), add = TRUE)
      progress_files <- file.path(progress_dir, paste0("chain_", seq_len(chains), ".txt"))
      .rtmb_progress_line("Preparing parallel workers...", progress_mode)
      futures <- vector("list", chains)
      line_counts <- integer(length(progress_files))
      for (c in seq_len(chains)) {
        futures[[c]] <- local({
          chain_id <- c
          progress_file <- progress_files[chain_id]
          future::future({
          write_progress_file <- function(path, msg) {
            if (is.numeric(msg)) msg <- ""
            msg <- as.character(msg[1])
            if (!nzchar(msg)) return(invisible(FALSE))

            ok <- tryCatch({
              cat(msg, "\n", file = path, append = TRUE, sep = "")
              TRUE
            }, error = function(e) FALSE, warning = function(w) FALSE)
            invisible(isTRUE(ok))
          }
          write_progress <- function(msg = "", amt = 1, ...) {
            if (is.numeric(msg)) { amt <- msg; msg <- "" }
            write_progress_file(progress_file, msg)
          }
          withCallingHandlers({
            run_chain(chain_id, f_ad = f_ad_global, p_callback = write_progress)
          }, warning = function(w) {
            if (grepl("package:BayesRTMB", conditionMessage(w))) invokeRestart("muffleWarning")
          })
          }, seed = TRUE, packages = "RTMB", globals = future_globals(list(chain_id = chain_id, progress_file = progress_file)))
        })
        line_counts <- .rtmb_report_progress_files(progress_files, line_counts)
      }
      results_list <- .rtmb_collect_progress_futures(futures, progress_files, line_counts = line_counts)
    }
  } else {
    iter <- sampling + warmup
    total_updates <- chains * (2 + floor(iter / 200))
    meter <- .rtmb_progress_meter(total_updates, progress_mode, label = "sampling")
    on.exit(meter$finish(), add = TRUE)
    results_list <- lapply(1:chains, function(c) {
      run_chain(c, f_ad = f_ad_global, p_callback = function(msg = "", amt = 1, ...) {
        if (is.numeric(msg)) { amt <- msg; msg <- "" }
        meter$advance(amt, msg = as.character(msg))
      })
    })
    meter$finish()
  }

  mcmc_index <- seq(from = (warmup + 1), to = (warmup + sampling), by = thin)
  accept_mat <- array(NA, dim = c(length(mcmc_index), chains))
  td_mat <- array(NA, dim = c(length(mcmc_index), chains))
  leapfrog_mat <- array(NA, dim = c(length(mcmc_index), chains))
  divergent_mat <- array(FALSE, dim = c(length(mcmc_index), chains))
  energy_mat <- array(NA, dim = c(length(mcmc_index), chains))
  eps_vec <- numeric(chains)
  metric_list <- vector("list", chains)
  metric_effective <- character(chains)
  metric_auto <- vector("list", chains)
  warmup_diagnostics <- vector("list", chains)
  pd_error_counts <- integer(chains)
  fit <- array(NA, dim = c(length(mcmc_index), chains, P_fixed + 1))
  dimnames(fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = c("lp", pl_fixed$names))

  if (P_random > 0) {
    random_fit <- array(NA, dim = c(length(mcmc_index), chains, P_random))
    dimnames(random_fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = pl_random$names)
  } else { random_fit <- NULL }

  for (c in 1:chains) {
    res <- results_list[[c]]
    pd_error_counts[c] <- if (!is.null(res$pd_error_count)) res$pd_error_count else 0L
    fit[, c, 1] <- res$lp[mcmc_index]
    for (j in 1:P_fixed) fit[, c, j + 1] <- res$para[, fixed_idx[j]]
    if (P_random > 0) { for (j in 1:P_random) random_fit[, c, j] <- res$para[, random_idx[j]] }
    accept_mat[, c] <- res$accept[mcmc_index]; td_mat[, c] <- res$treedepth[mcmc_index]; eps_vec[c] <- res$eps
    leapfrog_mat[, c] <- res$n_leapfrog[mcmc_index]
    divergent_mat[, c] <- res$divergent[mcmc_index]
    energy_mat[, c] <- res$energy[mcmc_index]
    metric_list[[c]] <- res$metric
    metric_effective[c] <- if (!is.null(res$metric_type)) res$metric_type else metric
    metric_auto[c] <- list(res$metric_auto)
    warmup_diagnostics[[c]] <- res$warmup_diagnostics
  }
  if (all(!is.finite(fit[, , "lp"]))) {
    stop_nonfinite_lp("MCMC retained samples")
  }
  if (any(!is.finite(fit[, , "lp"]))) {
    warning(
      "Some retained MCMC draws have non-finite log-probability values. ",
      "They will be ignored by summaries that use finite values, but diagnostics should be interpreted carefully.",
      call. = FALSE
    )
  }
  n_divergent <- sum(divergent_mat, na.rm = TRUE)
  n_transition <- sum(!is.na(divergent_mat))
  if (n_divergent > 0L && n_transition > 0L) {
    warning(
      sprintf(
        "%d of %d (%.2f%%) post-warmup MCMC transitions ended with a divergence. ",
        n_divergent, n_transition, 100 * n_divergent / n_transition
      ),
      "Try increasing delta, checking parameterization, or inspecting divergent draws with diagnose().",
      call. = FALSE
    )
  }
  if (sum(pd_error_counts) > 0L) {
    warning(
      "During post-warmup MCMC sampling, ",
      sum(pd_error_counts),
      " matrix positive-definite/singularity errors were treated as lp = -Inf.",
      call. = FALSE
    )
  }

  eps_chains <- eps_vec; accept_chains <- apply(accept_mat, 2, mean); treedepth_chains <- apply(td_mat, 2, max)
  names(eps_chains) <- names(accept_chains) <- names(treedepth_chains) <- paste0("chain", 1:chains)
  names(metric_effective) <- paste0("chain", 1:chains)
  metric_type <- if (length(unique(metric_effective)) == 1L) unname(unique(metric_effective)) else "mixed"

  posterior_mean <- numeric(length(self$pl_full$names)); names(posterior_mean) <- self$pl_full$names
  fixed_mean <- apply(fit[, , -1, drop = FALSE], 3, mean); posterior_mean[names(fixed_mean)] <- fixed_mean
  if (!is.null(random_fit)) { random_mean <- apply(random_fit, 3, mean); posterior_mean[names(random_mean)] <- random_mean }

  if (!is.null(save_info)) {
    for (c in 1:chains) {
      backup_file <- file.path(save_info$dir, paste0(save_info$name, "-", c, ".csv"))
      df_metrics <- data.frame(iteration = mcmc_index, accept = accept_mat[, c], treedepth = td_mat[, c], n_leapfrog = leapfrog_mat[, c], divergent = divergent_mat[, c], energy = energy_mat[, c], eps = eps_vec[c])
      df_out <- if (!is.null(random_fit)) cbind(df_metrics, as.data.frame(fit[, c, ]), as.data.frame(random_fit[, c, ])) else cbind(df_metrics, as.data.frame(fit[, c, ]))
      write.csv(df_out, file = backup_file, row.names = FALSE)
      warmup_file <- file.path(save_info$dir, paste0(save_info$name, "-", c, "-warmup.csv"))
      if (is.data.frame(warmup_diagnostics[[c]]) && nrow(warmup_diagnostics[[c]]) > 0L) {
        write.csv(warmup_diagnostics[[c]], file = warmup_file, row.names = FALSE)
      }
    }
  }

  res_obj <- MCMC_Fit$new(
    model = self, fit = fit, random_fit = random_fit,
    eps = eps_chains, accept = accept_chains, treedepth = treedepth_chains,
    laplace = laplace, posterior_mean = posterior_mean,
    max_treedepth = max_treedepth,
    pd_error_count = pd_error_counts,
    n_leapfrog = leapfrog_mat,
    divergent = divergent_mat,
    energy = energy_mat,
    metric = metric_list,
    metric_type = metric_type,
    metric_requested = metric,
    metric_effective = metric_effective,
    metric_auto = metric_auto,
    metric_init = metric_init,
    metric_adaptation = metric_adaptation,
    nuts_variant = nuts_variant,
    warmup_diagnostics = warmup_diagnostics
  )
  if (!is.null(self$transform)) res_obj$transformed_draws(self$transform, progress = progress_mode)
  if (!is.null(self$generate)) res_obj$generated_quantities(self$code$generate, progress = progress_mode)
  return(res_obj)
}

#' @noRd
.rtmb_vb_chain_score <- function(elbo_final_vec, rel_obj_vec, fit) {
  score <- elbo_final_vec
  nonfinite_score <- !is.finite(score)
  if (any(nonfinite_score)) {
    lp_idx <- match("lp", dimnames(fit)[[3]])
    if (!is.na(lp_idx)) {
      for (c in which(nonfinite_score)) {
        lp_vals <- fit[, c, lp_idx]
        lp_vals <- lp_vals[is.finite(lp_vals)]
        if (length(lp_vals) > 0L) {
          score[c] <- mean(lp_vals)
        }
      }
    }
  }
  if (!any(is.finite(score))) {
    finite_rel <- is.finite(rel_obj_vec)
    if (any(finite_rel)) {
      score[finite_rel] <- -rel_obj_vec[finite_rel]
    }
  }
  score
}

#' @noRd
.variational_impl <- function(self, private, iter, tol_rel_obj, window_size, num_samples, 
                              num_estimate, alpha, laplace, print_freq, 
                              method = c("meanfield", "fullrank", "hybrid"), 
                              parallel, seed, init, save_csv, map, fixed, globals, progress) {
  if (!is.null(fixed)) {
    return(private$.dispatch_fixed(.method_to_call = "variational", iter = iter, tol_rel_obj = tol_rel_obj, window_size = window_size,
      num_samples = num_samples, num_estimate = num_estimate, alpha = alpha,
      laplace = laplace, print_freq = print_freq, method = method,
      parallel = parallel, seed = seed, init = init, save_csv = save_csv, map = map, fixed = fixed, globals = globals, progress = progress))
  }

  set.seed(seed); method <- match.arg(method)
  progress_mode <- .rtmb_resolve_progress(progress)
  if (!is.null(save_csv)) {
    save_name <- if (!is.null(save_csv$name)) save_csv$name else "model_vb"
    save_dir <- if (!is.null(save_csv$dir)) save_csv$dir else "BayesRTMB_vb"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    save_info <- list(name = save_name, dir = save_dir)
  } else { save_info <- NULL }

  if (parallel && num_estimate > 1) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop(
        "parallel = TRUE for variational() requires the suggested package 'future'. ",
        "Please install it or use parallel = FALSE.",
        call. = FALSE
      )
    }
    if (inherits(future::plan(), "sequential")) {
      if (.Platform$OS.type == "unix") future::plan(future::multicore, workers = num_estimate)
      else future::plan(future::multisession, workers = num_estimate)
    }
    .rtmb_progress_start_line(paste0("Starting parallel VB estimation (num_estimate = ", num_estimate, ")..."))
  } else {
    .rtmb_progress_start_line(paste0("Starting sequential VB estimation (num_estimate = ", num_estimate, ")..."))
  }

  # --- Data extraction to avoid serialization ---
  local_data <- self$data; local_par_list <- self$par_list; local_pl_full <- self$pl_full
  local_map <- if (is.null(map)) self$map else utils::modifyList(as.list(self$map), as.list(map))
  # Normalize map factor levels for consistent MakeADFun behavior
  if (!is.null(local_map)) {
    local_map <- lapply(local_map, function(m) {
      if (is.null(m)) return(NULL)
      vals <- as.character(m)
      factor(vals, levels = unique(vals[!is.na(vals)]))
    })
  }
  local_log_prob <- self$log_prob; local_transform <- self$transform
  local_fixed_prior_specs <- self$fixed_prior_specs
  local_code_model <- self$code$model
  
  # Pre-calculate shadow list for fixed-prior removal
  env_shadow_local <- list2env(local_data, parent = asNamespace("BayesRTMB"))
  idx_s <- 1
  for (nm in names(local_par_list)) {
    p <- local_par_list[[nm]]
    shadow_vec <- local_pl_full$names[idx_s:(idx_s + p$length - 1)]
    if (length(p$dim) > 1) dim(shadow_vec) <- p$dim
    env_shadow_local[[nm]] <- shadow_vec
    idx_s <- idx_s + p$length
  }
  if (!is.null(local_transform)) {
    par_only <- mget(names(local_par_list), envir = env_shadow_local, inherits = FALSE)
    tran_res_shadow <- tryCatch(suppressMessages(local_transform(local_data, par_only)), error = function(e) list())
    for (nm in names(tran_res_shadow)) env_shadow_local[[nm]] <- tran_res_shadow[[nm]]
  }
  local_shadow_list <- as.list(env_shadow_local)

  random_effs <- names(local_par_list)[sapply(local_par_list, function(x) isTRUE(x$random))]
  use_random <- if (laplace && length(random_effs) > 0) random_effs else NULL
  base_init <- self$prepare_init(init)

  # Prepare f_ad beforehand to avoid capturing private environment
  f_ad_global <- private$.build_f_ad(
    data_local = local_data, 
    par_list_local = local_par_list, 
    log_prob_local = local_log_prob, 
    transform_local = local_transform, 
    jacobian_target = "all", 
    adreport = FALSE,
    fixed_prior_specs_local = local_fixed_prior_specs,
    code_model_local = local_code_model,
    shadow_list_local = local_shadow_list
  )
  environment(f_ad_global) <- list2env(list(
    data_local = local_data,
    par_list_local = local_par_list,
    log_prob_local = local_log_prob,
    transform_local = local_transform,
    jacobian_target = "all",
    adreport = FALSE,
    fixed_prior_specs_local = local_fixed_prior_specs,
    code_model_local = local_code_model,
    shadow_list_local = local_shadow_list
  ), parent = asNamespace("BayesRTMB"))

  run_advi_worker <- function(c, f_ad, p_callback = NULL, p_interval = 0) {
    unc_init_list <- to_unconstrained(constrained_vector_to_list(base_init, local_par_list), local_par_list)
    unc_init_vec <- unlist(unc_init_list, use.names = FALSE)
    unc_init_list_new <- unconstrained_vector_to_list(unc_init_vec, local_par_list)

    # Use the f_ad provided from outside
    ad_obj <- tryCatch({
      RTMB::MakeADFun(func = f_ad, parameters = unc_init_list_new, random = use_random, map = local_map, silent = TRUE)
    }, error = function(e) stop(.rtmb_format_makeadfun_error(e$message, context = "MakeADFun in parallel worker"), call. = FALSE))

    if (!is.null(use_random)) {
      orig_fn <- ad_obj$fn; orig_gr <- ad_obj$gr; idx_fixed_mask <- ad_obj$env$lfixed(); n_fixed <- sum(idx_fixed_mask)
      ad_obj$fn <- function(x, ...) { if (length(x) != n_fixed) x <- x[idx_fixed_mask]; orig_fn(x, ...) }
      ad_obj$gr <- function(x, ...) { if (length(x) != n_fixed) { g <- orig_gr(x[idx_fixed_mask], ...); g_full <- rep(0, length(x)); g_full[idx_fixed_mask] <- g; return(g_full) }; return(orig_gr(x, ...)) }
    }

    if (!is.null(p_callback)) {
      p_callback(msg = paste0("est", c, " started..."), amt = 1)
    }

    callback_count <- 0L
    advi_progress <- NULL
    if (!is.null(p_callback) && p_interval > 0L) {
      advi_progress <- function(amount = 1, ...) {
        amount <- max(1L, as.integer(amount))
        callback_count <<- callback_count + amount
        iter_now <- min(iter, callback_count * p_interval)
        p_callback(msg = paste0("est", c, ": iter ", iter_now), amt = amount)
      }
    }

    res <- ADVI_method(model = ad_obj, par_list = local_par_list, pl_full = local_pl_full, iter = iter, tol_rel_obj = tol_rel_obj, window_size = window_size, num_samples = num_samples, alpha = alpha, laplace = laplace, print_freq = if (is.null(advi_progress)) print_freq else 0, method = method, update_progress = advi_progress, update_interval = p_interval)
    if (!is.null(p_callback)) {
      p_callback(msg = paste0("est", c, " done"), amt = 1)
    }
    return(res)
  }

  worker_env <- list2env(list(
    base_init = base_init,
    local_par_list = local_par_list,
    use_random = use_random,
    local_map = local_map,
    iter = iter,
    tol_rel_obj = tol_rel_obj,
    window_size = window_size,
    num_samples = num_samples,
    alpha = alpha,
    laplace = laplace,
    print_freq = print_freq,
    method = method,
    local_pl_full = local_pl_full
  ), parent = asNamespace("BayesRTMB"))
  environment(run_advi_worker) <- worker_env

  future_globals <- function(extra = list()) {
    if (isTRUE(globals)) return(TRUE)
    c(list(run_advi_worker = run_advi_worker, f_ad_global = f_ad_global), extra)
  }

  results_list <- list()
  update_interval <- max(1L, floor(iter / 5L))
  if (parallel && num_estimate > 1) {
    if (identical(progress_mode, "none")) {
      results_list <- withCallingHandlers({
        futures <- lapply(seq_len(num_estimate), function(c) {
          estimate_id <- c
          future::future({
            run_advi_worker(estimate_id, f_ad = f_ad_global, p_callback = NULL, p_interval = 0L)
          }, seed = TRUE, packages = "RTMB", globals = future_globals(list(estimate_id = estimate_id)))
        })
        lapply(futures, future::value)
      }, warning = function(w) { if (grepl("BayesRTMB", conditionMessage(w))) invokeRestart("muffleWarning") })
    } else {
      progress_dir <- .rtmb_progress_file_dir()
      on.exit(unlink(progress_dir, recursive = TRUE, force = TRUE), add = TRUE)
      progress_files <- file.path(progress_dir, paste0("est_", seq_len(num_estimate), ".txt"))
      .rtmb_progress_line("Preparing parallel VB workers...", progress_mode)
      futures <- vector("list", num_estimate)
      line_counts <- integer(length(progress_files))
      for (c in seq_len(num_estimate)) {
        futures[[c]] <- local({
          estimate_id <- c
          progress_file <- progress_files[estimate_id]
          future::future({
          write_progress_file <- function(path, msg) {
            if (is.numeric(msg)) msg <- ""
            msg <- as.character(msg[1])
            if (!nzchar(msg)) return(invisible(FALSE))

            ok <- tryCatch({
              cat(msg, "\n", file = path, append = TRUE, sep = "")
              TRUE
            }, error = function(e) FALSE, warning = function(w) FALSE)
            invisible(isTRUE(ok))
          }
          write_progress <- function(msg = "", amt = 1, ...) {
            if (is.numeric(msg)) { amt <- msg; msg <- "" }
            write_progress_file(progress_file, msg)
          }
          withCallingHandlers({
            run_advi_worker(estimate_id, f_ad = f_ad_global, p_callback = write_progress, p_interval = update_interval)
          }, warning = function(w) {
            if (grepl("BayesRTMB", conditionMessage(w))) invokeRestart("muffleWarning")
          })
          }, seed = TRUE, packages = "RTMB", globals = future_globals(list(
            estimate_id = estimate_id,
            progress_file = progress_file,
            update_interval = update_interval
          )))
        })
        line_counts <- .rtmb_report_progress_files(progress_files, line_counts)
      }
      results_list <- .rtmb_collect_progress_futures(futures, progress_files, line_counts = line_counts)
    }
  } else {
    total_updates <- num_estimate * (2 + ceiling(iter / update_interval))
    meter <- .rtmb_progress_meter(total_updates, progress_mode, label = "VB estimation")
    on.exit(meter$finish(), add = TRUE)
    results_list <- lapply(seq_len(num_estimate), function(c) {
      run_advi_worker(c, f_ad = f_ad_global, p_callback = function(msg = "", amt = 1, ...) {
        if (is.numeric(msg)) { amt <- msg; msg <- "" }
        meter$advance(amt, msg = as.character(msg))
      }, p_interval = update_interval)
    })
    meter$finish()
  }

  P_fixed <- dim(results_list[[1]]$fit)[3] - 1; P_random <- if (!is.null(results_list[[1]]$random_fit)) dim(results_list[[1]]$random_fit)[3] else 0
  fit <- array(NA, dim = c(num_samples, num_estimate, P_fixed + 1)); dimnames(fit) <- list(iteration = NULL, chain = paste0("est", 1:num_estimate), variable = dimnames(results_list[[1]]$fit)[[3]])
  if (P_random > 0) { random_fit <- array(NA, dim = c(num_samples, num_estimate, P_random)); dimnames(random_fit) <- list(iteration = NULL, chain = paste0("est", 1:num_estimate), variable = dimnames(results_list[[1]]$random_fit)[[3]]) } else { random_fit <- NULL }

  elbo_history_list <- list(); elbo_final_vec <- numeric(num_estimate); rel_obj_vec <- numeric(num_estimate)
  for (c in 1:num_estimate) {
    res <- results_list[[c]]; fit[, c, ] <- res$fit[, 1, ]
    if (P_random > 0) random_fit[, c, ] <- res$random_fit[, 1, ]
    elbo_history_list[[c]] <- res$elbo_history - self$prior_correction
    elbo_final_vec[c] <- res$elbo_final - self$prior_correction; rel_obj_vec[c] <- res$rel_obj_final
  }

  if (!is.null(save_info)) {
    for (c in 1:num_estimate) {
      backup_file <- file.path(save_info$dir, paste0(save_info$name, "-", c, ".csv"))
      df_out <- if (!is.null(random_fit)) cbind(iteration = 1:num_samples, as.data.frame(fit[, c, ]), as.data.frame(random_fit[, c, ])) else cbind(iteration = 1:num_samples, as.data.frame(fit[, c, ]))
      write.csv(df_out, file = backup_file, row.names = FALSE)
    }
  }

  chain_score <- .rtmb_vb_chain_score(elbo_final_vec, rel_obj_vec, fit)
  if (!any(is.finite(chain_score))) {
    stop(
      "VB estimation finished, but no estimate had a finite ELBO, finite lp draw, or finite rel_obj. ",
      "Try smaller 'alpha', stronger priors, different initial values, or checking the model parameterization.",
      call. = FALSE
    )
  }
  best_chain <- which.max(chain_score)
  if (!is.finite(elbo_final_vec[best_chain])) {
    warning(
      sprintf(
        "All or some final VB ELBO values were non-finite; selected est%d using fallback diagnostics.",
        best_chain
      ),
      call. = FALSE
    )
  }
  posterior_mean <- numeric(length(self$pl_full$names)); names(posterior_mean) <- self$pl_full$names
  fixed_mean <- apply(fit[, best_chain, -1, drop = FALSE], 3, mean); posterior_mean[names(fixed_mean)] <- fixed_mean
  if (!is.null(random_fit)) { random_mean <- apply(random_fit[, best_chain, , drop = FALSE], 3, mean); posterior_mean[names(random_mean)] <- random_mean }

  cat("\nConvergence Diagnostics per estimate:\n")
  for (c in 1:num_estimate) {
    status <- if (!is.na(rel_obj_vec[c]) && rel_obj_vec[c] < tol_rel_obj) "Converged" else "Not Converged"
    best_marker <- if (c == best_chain) "  <-- BEST" else ""
    cat(sprintf("  est%d: ELBO = %10.2f, Final rel_obj = %.5f (%s)%s\n", c, elbo_final_vec[c], rel_obj_vec[c], status, best_marker))
  }
  
  res_obj <- VB_Fit$new(model = self, fit = fit, random_fit = random_fit, elbo_history = elbo_history_list, laplace = laplace, posterior_mean = posterior_mean, ELBO = elbo_final_vec, rel_obj_vals = rel_obj_vec, best_chain = best_chain, mu_history = results_list[[best_chain]]$mu_history)
  if (!is.null(self$transform)) res_obj$transformed_draws(self$transform, progress = progress_mode)
  if (!is.null(self$generate)) {
    res_obj$generated_quantities(self$code$generate, progress = progress_mode)
  }
  return(res_obj)
}
