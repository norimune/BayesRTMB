.diagnostic_row <- function(check, status, message) {
  data.frame(
    check = check,
    status = status,
    message = message,
    stringsAsFactors = FALSE
  )
}

.recommendation_row <- function(check, priority, message) {
  data.frame(
    check = check,
    priority = priority,
    message = message,
    stringsAsFactors = FALSE
  )
}

.diagnostic_status <- function(checks) {
  if (is.null(checks) || nrow(checks) == 0L) return("ok")
  if (any(checks$status == "problem", na.rm = TRUE)) "problem"
  else if (any(checks$status == "warning", na.rm = TRUE)) "warning"
  else "ok"
}

.empty_recommendations <- function() {
  data.frame(
    check = character(),
    priority = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

.combine_recommendations <- function(recommendations) {
  recommendations <- Filter(Negate(is.null), recommendations)
  if (length(recommendations) == 0L) return(.empty_recommendations())
  out <- do.call(rbind, recommendations)
  priorities <- c(high = 1L, medium = 2L, low = 3L)
  ord <- priorities[out$priority]
  ord[is.na(ord)] <- 99L
  out[order(ord), , drop = FALSE]
}

.new_diagnose_result <- function(fit_type, checks, details = list(),
                                 recommendations = NULL) {
  recommendations <- recommendations %||% .empty_recommendations()
  out <- list(
    fit_type = fit_type,
    status = .diagnostic_status(checks),
    checks = checks,
    recommendations = recommendations,
    details = details
  )
  class(out) <- "diagnose_BayesRTMB"
  out
}

.has_generated_log_lik <- function(fit) {
  if (!is.null(fit$generate) && !is.null(fit$generate$log_lik)) return(TRUE)
  if (!is.null(fit$generate_fit)) {
    gnames <- dimnames(fit$generate_fit)[[3L]]
    if (any(gnames == "log_lik" | grepl("^log_lik\\[", gnames))) return(TRUE)
  }
  if (!is.null(fit$se_samples) && !is.null(fit$se_samples$gq) &&
      !is.null(fit$se_samples$gq$log_lik)) return(TRUE)
  FALSE
}

.summarize_se_tables <- function(...) {
  tabs <- list(...)
  tabs <- tabs[vapply(tabs, is.data.frame, logical(1))]
  if (length(tabs) == 0L) return(NULL)
  do.call(rbind, tabs)
}

.selected_optimizer_status <- function(fit) {
  hist <- fit$opt_history
  if (is.data.frame(hist) && nrow(hist) > 0L) {
    if ("selected" %in% names(hist) && any(hist$selected %in% TRUE)) {
      idx <- which(hist$selected %in% TRUE)[1L]
    } else if ("objective" %in% names(hist)) {
      idx <- suppressWarnings(which.min(hist$objective))
      if (length(idx) == 0L || is.na(idx)) idx <- 1L
    } else {
      idx <- 1L
    }

    status <- if ("status" %in% names(hist)) hist$status[idx][1L] else NA_character_
    code <- if ("code" %in% names(hist)) hist$code[idx][1L] else fit$convergence
    message <- if ("message" %in% names(hist)) hist$message[idx][1L] else NA_character_

    if (is.na(status) || !nzchar(status)) {
      status <- if (!is.na(code) && code == 0L) "converged" else "not converged"
    }

    return(list(status = status, code = code, message = message))
  }

  status <- if (identical(fit$convergence, 0L) || identical(fit$convergence, 0)) "converged" else "not converged"
  list(status = status, code = fit$convergence, message = NA_character_)
}

.optimizer_status_message <- function(opt) {
  status <- opt$status %||% "not converged"
  status <- sub("\\s*\\([0-9]+\\)\\s*$", "", status)
  msg <- opt$message %||% NA_character_
  code <- opt$code %||% NA

  if (identical(status, "converged")) return("optimizer converged")
  if (!identical(status, "not converged")) {
    return(paste("optimizer ended with", status))
  }

  detail <- paste0(
    "optimizer not converged",
    if (!is.na(code)) paste0(" (code: ", code, ")") else "",
    if (!is.na(msg) && nzchar(msg)) paste0("; message: ", msg) else ""
  )
  detail
}

.diagnostic_check_status <- function(checks, check) {
  if (is.null(checks) || !is.data.frame(checks) || nrow(checks) == 0L) return(NA_character_)
  idx <- which(checks$check == check)
  if (length(idx) == 0L) return(NA_character_)
  checks$status[idx[1L]]
}

.diagnostic_is_flagged <- function(checks, check) {
  status <- .diagnostic_check_status(checks, check)
  !is.na(status) && status %in% c("warning", "problem")
}

.diagnose_map_recommendations <- function(checks) {
  recs <- list()
  if (.diagnostic_is_flagged(checks, "convergence")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "convergence", "high",
      "The optimizer did not converge cleanly. Try different initial values, increasing num_estimate, or adjusting optimizer control settings."
    )
  }
  if (.diagnostic_is_flagged(checks, "finite_objective")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "finite_objective", "high",
      "The objective is not finite. Check data for NA/Inf values, invalid link-function inputs, boundary values, and overly weak or incompatible priors."
    )
  }
  if (.diagnostic_is_flagged(checks, "hessian") ||
      .diagnostic_is_flagged(checks, "se_fallback") ||
      .diagnostic_is_flagged(checks, "vcov_unc") ||
      .diagnostic_is_flagged(checks, "standard_errors")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "standard_errors", "medium",
      "Hessian-based uncertainty may be unreliable. Consider se_method = 'sampling', rescaling variables, checking boundary estimates, or simplifying/reparameterizing the model."
    )
  }
  .combine_recommendations(recs)
}

.diagnose_classic_recommendations <- function(checks) {
  recs <- list()
  if (.diagnostic_is_flagged(checks, "convergence")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "convergence", "high",
      "The optimizer did not converge cleanly. Try different initial values, increasing num_estimate, or adjusting optimizer control settings."
    )
  }
  if (.diagnostic_is_flagged(checks, "finite_logLik")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "finite_logLik", "high",
      "The log-likelihood is not finite. Check data for NA/Inf values, invalid link-function inputs, and boundary values."
    )
  }
  if (.diagnostic_is_flagged(checks, "standard_errors")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "standard_errors", "medium",
      "Some standard errors are non-finite. Check for boundary estimates, weak identification, or variables on very different scales."
    )
  }
  .combine_recommendations(recs)
}

.diagnose_vb_recommendations <- function(checks) {
  recs <- list()
  if (.diagnostic_is_flagged(checks, "ELBO")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "ELBO", "high",
      "The final ELBO contains non-finite values. Try different initial values, stronger priors, or checking the data and model code for NA/NaN-producing operations."
    )
  }
  if (.diagnostic_is_flagged(checks, "relative_objective")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "relative_objective", "medium",
      "The variational objective has not stabilized. Consider increasing iter or num_estimate, or using MCMC for a more reliable posterior check."
    )
  }
  if (.diagnostic_is_flagged(checks, "best_chain")) {
    recs[[length(recs) + 1L]] <- .recommendation_row(
      "best_chain", "medium",
      "No reliable best variational run was selected. Inspect ELBO traces, increase num_estimate, or retry with different initial values."
    )
  }
  .combine_recommendations(recs)
}

diagnose_map_fit <- function(fit, recommend = TRUE, ...) {
  checks <- list()
  opt <- .selected_optimizer_status(fit)
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "convergence",
    if (identical(opt$status, "converged")) "ok" else "warning",
    .optimizer_status_message(opt)
  )
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "finite_objective",
    if (is.finite(fit$objective)) "ok" else "problem",
    if (is.finite(fit$objective)) "objective is finite" else "objective is not finite"
  )

  if (identical(fit$ci_method, "none") || is.null(fit$sd_rep)) {
    checks[[length(checks) + 1L]] <- .diagnostic_row("hessian", "skipped", "SE calculation was not requested")
  } else {
    pd <- fit$sd_rep$pdHess
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "hessian",
      if (isTRUE(pd)) "ok" else "warning",
      if (isTRUE(pd)) "pdHess is TRUE" else "pdHess is not TRUE; Wald SE may be unreliable"
    )
  }

  if (!is.null(fit$fallback_needed)) {
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "se_fallback",
      if (isTRUE(fit$fallback_needed)) "warning" else "ok",
      if (isTRUE(fit$fallback_needed)) "Hessian/SE fallback was used" else "no Hessian/SE fallback recorded"
    )
  }

  if (!is.null(fit$vcov_unc)) {
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "vcov_unc",
      if (all(is.finite(fit$vcov_unc))) "ok" else "warning",
      if (all(is.finite(fit$vcov_unc))) "unconstrained covariance is finite" else "unconstrained covariance contains non-finite values"
    )
  }

  se_tab <- .summarize_se_tables(fit$df_fixed, fit$random_effects, fit$df_transform, fit$df_generate)
  if (!is.null(se_tab) && "Std. Error" %in% names(se_tab)) {
    se_vals <- se_tab[["Std. Error"]]
    bad <- !is.na(se_vals) & !is.finite(se_vals)
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "standard_errors",
      if (any(bad)) "warning" else "ok",
      if (any(bad)) paste(sum(bad), "standard errors are non-finite") else "reported standard errors are finite or NA by design"
    )
  }

  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "generated_log_lik",
    if (.has_generated_log_lik(fit)) "ok" else "skipped",
    if (.has_generated_log_lik(fit)) "pointwise log_lik is available for WAIC" else "pointwise log_lik is not available"
  )

  checks_df <- do.call(rbind, checks)
  recommendations <- if (isTRUE(recommend)) .diagnose_map_recommendations(checks_df) else .empty_recommendations()

  .new_diagnose_result("MAP_Fit", checks_df, list(
    marginal_vars = fit$marginal_vars,
    laplace_random_vars = fit$laplace_random_vars,
    ci_method = fit$ci_method
  ), recommendations)
}

diagnose_mcmc_fit <- function(fit, rhat_warning = 1.01, rhat_problem = 1.05,
                              ess_warning = 400, ess_problem = 100,
                              divergent_problem_rate = 0.01,
                              divergent_strong_problem_rate = 0.05,
                              recommend = TRUE, ...) {
  checks <- list()
  min_accept <- NA_real_
  max_td <- NA_real_
  max_allowed <- NA_real_
  n_divergent <- NA_integer_
  n_transition <- NA_integer_
  divergent_rate <- NA_real_
  max_cond <- NA_real_
  max_rhat <- NA_real_
  min_ess <- NA_real_
  pd_error_total <- NA_integer_
  metric_auto_diag_reason <- character(0)
  ndraws <- dim(fit$fit)[1L]
  nchains <- dim(fit$fit)[2L]
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "draw_count",
    if (ndraws * nchains < ess_problem) "warning" else "ok",
    paste("retained draws:", ndraws * nchains, "from", nchains, "chain(s)")
  )

  lp_idx <- which(dimnames(fit$fit)[[3L]] == "lp")
  if (length(lp_idx) == 1L) {
    lp <- fit$fit[, , lp_idx, drop = TRUE]
    nonfinite <- sum(!is.finite(lp))
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "finite_lp",
      if (nonfinite == 0L) "ok" else "warning",
      if (nonfinite == 0L) "all retained lp values are finite" else paste(nonfinite, "retained lp values are non-finite")
    )
  }

  if (!is.null(fit$accept)) {
    min_accept <- min(fit$accept, na.rm = TRUE)
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "acceptance",
      if (min_accept < 0.4) "problem" else if (min_accept < 0.6) "warning" else "ok",
      sprintf("minimum mean chain acceptance is %.3f", min_accept)
    )
  }

  if (!is.null(fit$treedepth)) {
    max_td <- max(fit$treedepth, na.rm = TRUE)
    max_allowed <- fit$max_treedepth %||% NA
    td_status <- if (is.finite(max_allowed) && max_td >= max_allowed) "warning" else "ok"
    td_msg <- if (is.finite(max_allowed)) {
      sprintf("maximum treedepth is %s (limit %s)", max_td, max_allowed)
    } else {
      sprintf("maximum treedepth is %s", max_td)
    }
    checks[[length(checks) + 1L]] <- .diagnostic_row("treedepth", td_status, td_msg)
  }

  if (!is.null(fit$divergent)) {
    n_divergent <- sum(fit$divergent, na.rm = TRUE)
    n_transition <- sum(!is.na(fit$divergent))
    divergent_rate <- if (n_transition > 0L) n_divergent / n_transition else NA_real_
    divergent_status <- if (n_divergent == 0L) {
      "ok"
    } else if (is.finite(divergent_rate) && divergent_rate >= divergent_problem_rate) {
      "problem"
    } else {
      "warning"
    }
    divergent_msg <- if (is.finite(divergent_rate)) {
      sprintf(
        "divergent transitions: %d of %d (%.2f%%)",
        n_divergent, n_transition, 100 * divergent_rate
      )
    } else {
      paste("divergent transitions:", n_divergent)
    }
    if (is.finite(divergent_rate) && divergent_rate >= divergent_strong_problem_rate) {
      divergent_msg <- paste0(divergent_msg, "; high divergence rate")
    }
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "divergent",
      divergent_status,
      divergent_msg
    )
  }

  if (!is.null(fit$n_leapfrog)) {
    max_leapfrog <- suppressWarnings(max(fit$n_leapfrog, na.rm = TRUE))
    mean_leapfrog <- suppressWarnings(mean(fit$n_leapfrog, na.rm = TRUE))
    if (is.finite(max_leapfrog) && is.finite(mean_leapfrog)) {
      checks[[length(checks) + 1L]] <- .diagnostic_row(
        "leapfrog",
        "ok",
        sprintf("mean leapfrog steps %.1f; maximum %.0f", mean_leapfrog, max_leapfrog)
      )
    }
  }

  if (!is.null(fit$metric_type)) {
    metric_requested <- if ("metric_requested" %in% names(fit)) fit$metric_requested else fit$metric_type
    metric_effective <- if ("metric_effective" %in% names(fit)) fit$metric_effective else NULL
    metric_label <- if (identical(metric_requested, "auto")) {
      paste("NUTS metric: auto ->", fit$metric_type)
    } else {
      paste("NUTS metric:", fit$metric_type)
    }
    if (!is.null(metric_effective) && length(unique(metric_effective)) > 1L) {
      metric_label <- paste0(
        metric_label,
        " (",
        paste(names(table(metric_effective)), as.integer(table(metric_effective)),
              sep = "=", collapse = ", "),
        ")"
      )
    }
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "metric",
      "ok",
      paste(
        metric_label,
        "initial:", fit$metric_init %||% "unknown",
        "adaptation:", fit$metric_adaptation %||% "unknown"
      )
    )
  }
  metric_auto <- if ("metric_auto" %in% names(fit)) fit$metric_auto else NULL
  if (!is.null(metric_auto) && length(metric_auto) > 0L) {
    auto_entries <- Filter(Negate(is.null), metric_auto)
    if (length(auto_entries) > 0L) {
      auto_effective <- vapply(auto_entries, function(x) x$effective %||% NA_character_, character(1))
      auto_reasons <- unique(vapply(auto_entries, function(x) x$reason %||% NA_character_, character(1)))
      auto_reasons <- auto_reasons[!is.na(auto_reasons)]
      msg <- paste0(
        "auto selected ",
        paste(names(table(auto_effective)), as.integer(table(auto_effective)),
              sep = "=", collapse = ", ")
      )
      if (length(auto_reasons) > 0L) {
        msg <- paste0(msg, "; ", paste(utils::head(auto_reasons, 2L), collapse = " | "))
      }
      checks[[length(checks) + 1L]] <- .diagnostic_row("metric_auto", "ok", msg)
      if (any(auto_effective == "diag", na.rm = TRUE)) {
        metric_auto_diag_reason <- utils::head(auto_reasons, 2L)
      }
    }
  }
  if (!is.null(fit$warmup_diagnostics) && length(fit$warmup_diagnostics) > 0L) {
    warmup_df <- do.call(rbind, fit$warmup_diagnostics)
    if (is.data.frame(warmup_df) && nrow(warmup_df) > 0L) {
      mean_accept_w <- suppressWarnings(mean(warmup_df$accept, na.rm = TRUE))
      final_eps_w <- suppressWarnings(stats::median(vapply(fit$warmup_diagnostics, function(x) {
        if (!is.data.frame(x) || nrow(x) == 0L) return(NA_real_)
        utils::tail(x$eps, 1L)
      }, numeric(1)), na.rm = TRUE))
      max_cond_w <- suppressWarnings(max(warmup_df$metric_condition, na.rm = TRUE))
      n_metric_updates <- if ("metric_updated" %in% names(warmup_df)) {
        sum(warmup_df$metric_updated, na.rm = TRUE)
      } else {
        NA_integer_
      }
      msg <- sprintf("mean warmup acceptance %.3f; median final warmup eps %.3g", mean_accept_w, final_eps_w)
      if (is.finite(max_cond_w)) {
        msg <- paste0(msg, sprintf("; max metric condition %.2e", max_cond_w))
      }
      if (is.finite(n_metric_updates)) {
        msg <- paste0(msg, sprintf("; metric updates %d", n_metric_updates))
      }
      checks[[length(checks) + 1L]] <- .diagnostic_row("warmup", "ok", msg)
    }
  }
  if (length(fit$metric) > 0L && any(vapply(fit$metric, function(m) {
    (is.list(m) && identical(m$type, "hybrid")) || is.matrix(m)
  }, logical(1)))) {
    cond_vals <- vapply(fit$metric, function(m) {
      if (is.list(m) && identical(m$type, "hybrid")) {
        eig <- numeric(0)
        if (length(m$dense_idx) > 0L) {
          eig <- c(eig, tryCatch(eigen(m$dense, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_))
        }
        if (length(m$diag_idx) > 0L) eig <- c(eig, m$diag)
      } else if (is.matrix(m)) {
        eig <- tryCatch(eigen(m, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
      } else {
        return(NA_real_)
      }
      eig <- eig[is.finite(eig) & eig > 0]
      if (length(eig) == 0L) return(NA_real_)
      max(eig) / min(eig)
    }, numeric(1))
    max_cond <- suppressWarnings(max(cond_vals, na.rm = TRUE))
    if (is.finite(max_cond)) {
      checks[[length(checks) + 1L]] <- .diagnostic_row(
        "metric_condition",
        if (max_cond > 1e8) "warning" else "ok",
        sprintf("maximum %s metric condition number %.2e", fit$metric_type, max_cond)
      )
    }
  }

  summ <- tryCatch(fit$summary(max_rows = NULL, inc_random = FALSE, inc_transform = FALSE, inc_generate = FALSE), error = function(e) NULL)
  if (is.data.frame(summ)) {
    if ("rhat" %in% names(summ)) {
      max_rhat <- suppressWarnings(max(summ$rhat, na.rm = TRUE))
      if (is.finite(max_rhat)) {
        checks[[length(checks) + 1L]] <- .diagnostic_row(
          "rhat",
          if (max_rhat > rhat_problem) "problem" else if (max_rhat > rhat_warning) "warning" else "ok",
          sprintf("maximum R-hat is %.3f", max_rhat)
        )
      }
    }
    ess_cols <- intersect(c("ess_bulk", "ess_tail"), names(summ))
    if (length(ess_cols) > 0L) {
      min_ess <- suppressWarnings(min(as.matrix(summ[ess_cols]), na.rm = TRUE))
      if (is.finite(min_ess)) {
        checks[[length(checks) + 1L]] <- .diagnostic_row(
          "ess",
          if (min_ess < ess_problem) "problem" else if (min_ess < ess_warning) "warning" else "ok",
          sprintf("minimum reported ESS is %.1f", min_ess)
        )
      }
    }
  }

  if (!is.null(fit$pd_error_count)) {
    pd_error_total <- sum(fit$pd_error_count, na.rm = TRUE)
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "pd_errors",
      if (pd_error_total > 0L) "warning" else "ok",
      paste("positive-definite/singularity fallback count:", pd_error_total)
    )
  }

  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "generated_log_lik",
    if (.has_generated_log_lik(fit)) "ok" else "skipped",
    if (.has_generated_log_lik(fit)) "pointwise log_lik is available for WAIC" else "pointwise log_lik is not available"
  )

  checks_df <- do.call(rbind, checks)
  recs <- list()
  if (isTRUE(recommend)) {
    if (is.finite(min_accept) && min_accept < 0.6) {
      recs[[length(recs) + 1L]] <- .recommendation_row(
        "acceptance",
        if (min_accept < 0.4) "high" else "medium",
        "Mean acceptance is low. Try increasing delta, and if this persists, inspect model geometry and parameterization."
      )
    }
    if (is.finite(max_td) && is.finite(max_allowed) && max_td >= max_allowed) {
      recs[[length(recs) + 1L]] <- .recommendation_row(
        "treedepth", "medium",
        "Some transitions reached max_treedepth. First check divergences and parameterization; if the fit is otherwise healthy, consider increasing max_treedepth."
      )
    }
    if (is.finite(n_divergent) && n_divergent > 0L) {
      recs[[length(recs) + 1L]] <- .recommendation_row(
        "divergent",
        if (is.finite(divergent_rate) && divergent_rate >= divergent_problem_rate) "high" else "medium",
        "Divergences were detected. Try increasing delta, for example delta = 0.9 or 0.95; if divergences remain, inspect parameterization and divergent draws."
      )
    }
    if (is.finite(max_cond) && max_cond > 1e8) {
      recs[[length(recs) + 1L]] <- .recommendation_row(
        "metric_condition", "medium",
        "The mass-matrix condition number is very large. Consider rescaling predictors/parameters, using metric = 'diag', or simplifying the model."
      )
    }
    if (is.finite(max_rhat) && max_rhat > rhat_warning) {
      recs[[length(recs) + 1L]] <- .recommendation_row(
        "rhat",
        if (max_rhat > rhat_problem) "high" else "medium",
        "R-hat suggests incomplete mixing. Consider running more iterations, using more dispersed initial values, or checking for multimodality/label switching."
      )
    }
    if (is.finite(min_ess) && min_ess < ess_warning) {
      recs[[length(recs) + 1L]] <- .recommendation_row(
        "ess",
        if (min_ess < ess_problem) "high" else "medium",
        "Some effective sample sizes are low. Consider increasing sampling iterations or improving parameterization/metric choice."
      )
    }
    if (is.finite(pd_error_total) && pd_error_total > 0L) {
      recs[[length(recs) + 1L]] <- .recommendation_row(
        "pd_errors", "medium",
        "Positive-definite or singularity fallbacks occurred. Consider metric = 'diag', rescaling variables, or checking whether the posterior geometry is weakly identified."
      )
    }
    if (length(metric_auto_diag_reason) > 0L &&
        (.diagnostic_is_flagged(checks_df, "ess") ||
         .diagnostic_is_flagged(checks_df, "rhat"))) {
      recs[[length(recs) + 1L]] <- .recommendation_row(
        "metric_auto", "low",
        paste0(
          "metric = 'auto' selected a diagonal metric for at least one chain",
          if (length(metric_auto_diag_reason) > 0L) paste0(" (", paste(metric_auto_diag_reason, collapse = " | "), ")") else "",
          ". If mixing remains poor, compare with metric = 'hybrid' or inspect model scaling."
        )
      )
    }
  }
  recommendations <- .combine_recommendations(recs)

  .new_diagnose_result("MCMC_Fit", checks_df, list(
    laplace = fit$laplace,
    metric = fit$metric_type,
    metric_requested = if ("metric_requested" %in% names(fit)) fit$metric_requested else fit$metric_type,
    metric_effective = if ("metric_effective" %in% names(fit)) fit$metric_effective else NULL,
    metric_auto = if ("metric_auto" %in% names(fit)) fit$metric_auto else NULL,
    metric_init = fit$metric_init,
    metric_adaptation = fit$metric_adaptation,
    nuts_variant = fit$nuts_variant
  ), recommendations)
}

diagnose_vb_fit <- function(fit, rel_obj_warning = 1e-3, rel_obj_problem = 1e-2,
                            recommend = TRUE, ...) {
  checks <- list()
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "ELBO",
    if (all(is.finite(fit$ELBO))) "ok" else "problem",
    if (all(is.finite(fit$ELBO))) "final ELBO values are finite" else "ELBO contains non-finite values"
  )
  if (!is.null(fit$rel_obj_vals)) {
    best_rel <- suppressWarnings(min(fit$rel_obj_vals, na.rm = TRUE))
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "relative_objective",
      if (best_rel > rel_obj_problem) "problem" else if (best_rel > rel_obj_warning) "warning" else "ok",
      sprintf("best final relative objective value is %.3g", best_rel)
    )
  }
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "best_chain",
    if (!is.null(fit$best_chain) && is.finite(fit$best_chain)) "ok" else "warning",
    paste("best chain:", fit$best_chain %||% NA)
  )
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "generated_log_lik",
    if (.has_generated_log_lik(fit)) "ok" else "skipped",
    if (.has_generated_log_lik(fit)) "pointwise log_lik is available for WAIC" else "pointwise log_lik is not available"
  )
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "variational_note",
    "skipped",
    "ADVI is an approximate posterior; uncertainty and WAIC may be optimistic relative to MCMC"
  )
  checks_df <- do.call(rbind, checks)
  recommendations <- if (isTRUE(recommend)) .diagnose_vb_recommendations(checks_df) else .empty_recommendations()

  .new_diagnose_result("VB_Fit", checks_df, list(laplace = fit$laplace), recommendations)
}

diagnose_classic_fit <- function(fit, recommend = TRUE, ...) {
  checks <- list()
  opt <- .selected_optimizer_status(fit)
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "convergence",
    if (identical(opt$status, "converged")) "ok" else "warning",
    .optimizer_status_message(opt)
  )
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "finite_logLik",
    if (is.finite(as.numeric(fit$logLik()))) "ok" else "warning",
    if (is.finite(as.numeric(fit$logLik()))) "logLik is finite" else "logLik is not finite or unavailable"
  )
  se_tab <- .summarize_se_tables(fit$fit, fit$random_effects, fit$df_transform, fit$df_generate)
  if (!is.null(se_tab) && "Std. Error" %in% names(se_tab)) {
    se_vals <- se_tab[["Std. Error"]]
    bad <- !is.na(se_vals) & !is.finite(se_vals)
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "standard_errors",
      if (any(bad)) "warning" else "ok",
      if (any(bad)) paste(sum(bad), "standard errors are non-finite") else "reported standard errors are finite or NA by design"
    )
  }
  checks[[length(checks) + 1L]] <- .diagnostic_row("WAIC", "skipped", "WAIC is not defined for Classic_Fit")
  checks_df <- do.call(rbind, checks)
  recommendations <- if (isTRUE(recommend)) .diagnose_classic_recommendations(checks_df) else .empty_recommendations()

  .new_diagnose_result("Classic_Fit", checks_df, list(df_method = fit$df_method), recommendations)
}

#' @export
print.diagnose_BayesRTMB <- function(x, ...) {
  cat("BayesRTMB diagnostics\n")
  cat("Fit:", x$fit_type, "\n")
  cat("Status:", x$status, "\n\n")
  checks <- x$checks
  if (!is.null(checks) && nrow(checks) > 0L) {
    for (i in seq_len(nrow(checks))) {
      cat(sprintf("[%s] %s: %s\n", checks$status[i], checks$check[i], checks$message[i]))
    }
  }
  recommendations <- x$recommendations
  if (!is.null(recommendations) && nrow(recommendations) > 0L) {
    cat("\nRecommendations:\n")
    for (i in seq_len(nrow(recommendations))) {
      cat(sprintf("* [%s] %s: %s\n",
                  recommendations$priority[i],
                  recommendations$check[i],
                  recommendations$message[i]))
    }
  }
  invisible(x)
}
