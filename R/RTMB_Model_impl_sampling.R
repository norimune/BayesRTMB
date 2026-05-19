#' @noRd
.sample_impl <- function(self, private, sampling, warmup, chains, thin, seed, delta, 
                         max_treedepth, parallel, laplace, init, init_jitter, 
                         save_csv, map, fixed) {
  if (!is.null(fixed)) {
    return(private$.dispatch_fixed(.method_to_call = "sample", sampling = sampling, warmup = warmup, chains = chains, thin = thin,
      seed = seed, delta = delta, max_treedepth = max_treedepth,
      parallel = parallel, laplace = laplace, init = init,
      init_jitter = init_jitter, save_csv = save_csv, map = map, fixed = fixed))
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
    na_vars <- data_na_summary(self$data)
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
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE) ||
        !requireNamespace("progressr", quietly = TRUE)) {
      stop(
        "parallel = TRUE for sample() requires the suggested packages ",
        "'future', 'future.apply', and 'progressr'. ",
        "Please install them or use parallel = FALSE.",
        call. = FALSE
      )
    }
    if (inherits(future::plan(), "sequential")) {
      if (.Platform$OS.type == "unix") future::plan(future::multicore, workers = chains)
      else future::plan(future::multisession, workers = chains)
    }
    cat(paste0("Starting parallel sampling (chains = ", chains, ")...\n"))
  } else {
    cat(paste0("Starting sequential sampling (chains = ", chains, ")...\n"))
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

  run_chain <- function(c, f_ad, p_callback = NULL) {
    library(BayesRTMB)
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
    }, error = function(e) stop("Failed to setup MakeADFun in parallel worker.\n[Error]: ", e$message, call. = FALSE))
    ad_obj <- wrap_mcmc_pd_errors(ad_obj)

    init_lp <- tryCatch(-ad_obj$fn(ad_obj$par), error = function(e) NA_real_)
    if (!is.finite(init_lp)) {
      stop_nonfinite_lp(sprintf("MCMC initialization for chain %d", c))
    }

    res <- NUTS_method(model = ad_obj, sampling = sampling, warmup = warmup, delta = delta, max_treedepth = max_treedepth, chain = c, update_progress = p_callback, laplace = laplace, save_info = save_info)
    if (is.null(res$lp) || all(!is.finite(res$lp))) {
      stop_nonfinite_lp(sprintf("MCMC sampling for chain %d", c))
    }

    P_all_true <- length(local_pl_full$names)
    iter_total <- sampling + warmup
    para_final <- array(NA, dim = c(iter_total, P_all_true))

    for (i in 1:iter_total) {
      x_in <- as.numeric(res$para_fixed[i, ])
      if (laplace && length(ad_obj$env$random) > 0) {
        ad_obj$fn(x_in)
        para_list <- ad_obj$env$parList()
      } else {
        para_list <- ad_obj$env$parList(x = x_in)
      }
      con_list <- to_constrained(para_list, local_par_list)
      para_final[i, ] <- unlist(con_list, use.names = FALSE)
    }
    res$para <- para_final
    res$para_full <- NULL
    return(res)
  }

  results_list <- list()
  if (parallel) {
    iter <- sampling + warmup
    total_updates <- chains * (1 + floor(iter / 100))
    progressr::with_progress({
      p <- progressr::progressor(steps = total_updates)
      results_list <- withCallingHandlers({
        future.apply::future_lapply(1:chains, function(c) {
          run_chain(c, f_ad = f_ad_global, p_callback = function(msg = "", amt = 1, ...) {
            if (is.numeric(msg)) { amt <- msg; msg <- "" }
            p(amount = amt, message = as.character(msg))
          })
        }, future.seed = TRUE, future.packages = c("RTMB", "BayesRTMB"), future.globals = TRUE)
      }, warning = function(w) { if (grepl("package:BayesRTMB", conditionMessage(w))) invokeRestart("muffleWarning") })
    })
  } else {
    results_list <- lapply(1:chains, function(c) run_chain(c, f_ad = f_ad_global, p_callback = NULL))
  }

  mcmc_index <- seq(from = (warmup + 1), to = (warmup + sampling), by = thin)
  accept_mat <- array(NA, dim = c(length(mcmc_index), chains))
  td_mat <- array(NA, dim = c(length(mcmc_index), chains))
  eps_vec <- numeric(chains)
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
    for (j in 1:P_fixed) fit[, c, j + 1] <- res$para[mcmc_index, fixed_idx[j]]
    if (P_random > 0) { for (j in 1:P_random) random_fit[, c, j] <- res$para[mcmc_index, random_idx[j]] }
    accept_mat[, c] <- res$accept[mcmc_index]; td_mat[, c] <- res$treedepth[mcmc_index]; eps_vec[c] <- res$eps
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

  posterior_mean <- numeric(length(self$pl_full$names)); names(posterior_mean) <- self$pl_full$names
  fixed_mean <- apply(fit[, , -1, drop = FALSE], 3, mean); posterior_mean[names(fixed_mean)] <- fixed_mean
  if (!is.null(random_fit)) { random_mean <- apply(random_fit, 3, mean); posterior_mean[names(random_mean)] <- random_mean }

  if (!is.null(save_info)) {
    for (c in 1:chains) {
      backup_file <- file.path(save_info$dir, paste0(save_info$name, "-", c, ".csv"))
      df_metrics <- data.frame(iteration = mcmc_index, accept = accept_mat[, c], treedepth = td_mat[, c], eps = eps_vec[c])
      df_out <- if (!is.null(random_fit)) cbind(df_metrics, as.data.frame(fit[, c, ]), as.data.frame(random_fit[, c, ])) else cbind(df_metrics, as.data.frame(fit[, c, ]))
      write.csv(df_out, file = backup_file, row.names = FALSE)
    }
  }

  res_obj <- MCMC_Fit$new(
    model = self, fit = fit, random_fit = random_fit,
    eps = eps_chains, accept = accept_chains, treedepth = treedepth_chains,
    laplace = laplace, posterior_mean = posterior_mean,
    max_treedepth = max_treedepth,
    pd_error_count = pd_error_counts
  )
  if (!is.null(self$transform)) res_obj$transformed_draws(self$transform)
  if (!is.null(self$generate)) res_obj$generated_quantities(self$code$generate)
  return(res_obj)
}

#' @noRd
.variational_impl <- function(self, private, iter, tol_rel_obj, window_size, num_samples, 
                              num_estimate, alpha, laplace, print_freq, 
                              method = c("meanfield", "fullrank", "hybrid"), 
                              parallel, seed, init, save_csv, map, fixed) {
  if (!is.null(fixed)) {
    return(private$.dispatch_fixed(.method_to_call = "variational", iter = iter, tol_rel_obj = tol_rel_obj, window_size = window_size,
      num_samples = num_samples, num_estimate = num_estimate, alpha = alpha,
      laplace = laplace, print_freq = print_freq, method = method,
      parallel = parallel, seed = seed, init = init, save_csv = save_csv, map = map, fixed = fixed))
  }

  set.seed(seed); method <- match.arg(method)
  if (!is.null(save_csv)) {
    save_name <- if (!is.null(save_csv$name)) save_csv$name else "model_vb"
    save_dir <- if (!is.null(save_csv$dir)) save_csv$dir else "BayesRTMB_vb"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    save_info <- list(name = save_name, dir = save_dir)
  } else { save_info <- NULL }

  if (parallel && num_estimate > 1) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop(
        "parallel = TRUE for variational() requires the suggested packages ",
        "'future' and 'future.apply'. Please install them or use parallel = FALSE.",
        call. = FALSE
      )
    }
    if (inherits(future::plan(), "sequential")) {
      if (.Platform$OS.type == "unix") future::plan(future::multicore, workers = num_estimate)
      else future::plan(future::multisession, workers = num_estimate)
    }
    cat(paste0("Starting parallel VB estimation (num_estimate = ", num_estimate, ")...\n"))
  } else { cat(paste0("Starting sequential VB estimation (num_estimate = ", num_estimate, ")...\n")) }

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

  run_advi_worker <- function(c, f_ad, p_callback = NULL, p_interval = 0) {
    library(BayesRTMB)
    unc_init_list <- to_unconstrained(constrained_vector_to_list(base_init, local_par_list), local_par_list)
    unc_init_vec <- unlist(unc_init_list, use.names = FALSE)
    unc_init_list_new <- unconstrained_vector_to_list(unc_init_vec, local_par_list)

    # Use the f_ad provided from outside
    ad_obj <- tryCatch({
      RTMB::MakeADFun(func = f_ad, parameters = unc_init_list_new, random = use_random, map = local_map, silent = TRUE)
    }, error = function(e) stop("Failed to setup MakeADFun in parallel worker.\n[Error]: ", e$message, call. = FALSE))

    if (!is.null(use_random)) {
      orig_fn <- ad_obj$fn; orig_gr <- ad_obj$gr; idx_fixed_mask <- ad_obj$env$lfixed(); n_fixed <- sum(idx_fixed_mask)
      ad_obj$fn <- function(x, ...) { if (length(x) != n_fixed) x <- x[idx_fixed_mask]; orig_fn(x, ...) }
      ad_obj$gr <- function(x, ...) { if (length(x) != n_fixed) { g <- orig_gr(x[idx_fixed_mask], ...); g_full <- rep(0, length(x)); g_full[idx_fixed_mask] <- g; return(g_full) }; return(orig_gr(x, ...)) }
    }

    res <- ADVI_method(model = ad_obj, par_list = local_par_list, pl_full = local_pl_full, iter = iter, tol_rel_obj = tol_rel_obj, window_size = window_size, num_samples = num_samples, alpha = alpha, laplace = laplace, print_freq = if (is.null(p_callback)) print_freq else 0, method = method, update_progress = p_callback, update_interval = p_interval)
    return(res)
  }

  results_list <- list()
  if (parallel && num_estimate > 1) {
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::handlers(global = TRUE); update_interval <- max(1, floor(iter / 100)); total_steps <- ceiling(iter / update_interval) * num_estimate
      results_list <- progressr::with_progress({
        p <- progressr::progressor(steps = total_steps)
        withCallingHandlers({
          future.apply::future_lapply(1:num_estimate, function(c) {
            run_advi_worker(c = c, f_ad = f_ad_global, p_callback = function(amount = 1) p(amount = amount), p_interval = update_interval)
          }, future.seed = TRUE, future.packages = c("RTMB", "BayesRTMB"), future.globals = TRUE)
        }, warning = function(w) { if (grepl("BayesRTMB", conditionMessage(w))) invokeRestart("muffleWarning") })
      })
    } else {
      results_list <- future.apply::future_lapply(1:num_estimate, function(c) run_advi_worker(c, f_ad = f_ad_global, p_callback = NULL, p_interval = 0), future.seed = TRUE, future.packages = c("RTMB", "BayesRTMB"), future.globals = TRUE)
    }
  } else { results_list <- lapply(1:num_estimate, function(c) { if (print_freq > 0) cat(sprintf("\n--- Starting VB estimation: est%d ---\n", c)); run_advi_worker(c, f_ad = f_ad_global, p_callback = NULL, p_interval = 0) }) }

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

  best_chain <- which.max(elbo_final_vec)
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
  if (!is.null(self$transform)) { cat("Calculating transformed parameters...\n"); res_obj$transformed_draws(self$transform) }
  if (!is.null(self$generate)) {
    # --- Post-calculation of generated quantities ---
    # Execute the model's generate block for each posterior sample
    cat("Calculating generated quantities...\n"); res_obj$generated_quantities(self$code$generate) 
  }
  return(res_obj)
}
