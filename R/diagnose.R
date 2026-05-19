.diagnostic_row <- function(check, status, message) {
  data.frame(
    check = check,
    status = status,
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

.new_diagnose_result <- function(fit_type, checks, details = list()) {
  out <- list(
    fit_type = fit_type,
    status = .diagnostic_status(checks),
    checks = checks,
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

diagnose_map_fit <- function(fit, ...) {
  checks <- list()
  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "convergence",
    if (identical(fit$convergence, 0L) || identical(fit$convergence, 0)) "ok" else "warning",
    if (identical(fit$convergence, 0L) || identical(fit$convergence, 0)) "optimizer converged" else paste("optimizer convergence code:", fit$convergence)
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

  .new_diagnose_result("MAP_Fit", do.call(rbind, checks), list(
    marginal_vars = fit$marginal_vars,
    laplace_random_vars = fit$laplace_random_vars,
    ci_method = fit$ci_method
  ))
}

diagnose_mcmc_fit <- function(fit, rhat_warning = 1.01, rhat_problem = 1.05,
                              ess_warning = 400, ess_problem = 100, ...) {
  checks <- list()
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
    checks[[length(checks) + 1L]] <- .diagnostic_row(
      "pd_errors",
      if (sum(fit$pd_error_count, na.rm = TRUE) > 0L) "warning" else "ok",
      paste("positive-definite/singularity fallback count:", sum(fit$pd_error_count, na.rm = TRUE))
    )
  }

  checks[[length(checks) + 1L]] <- .diagnostic_row(
    "generated_log_lik",
    if (.has_generated_log_lik(fit)) "ok" else "skipped",
    if (.has_generated_log_lik(fit)) "pointwise log_lik is available for WAIC" else "pointwise log_lik is not available"
  )

  .new_diagnose_result("MCMC_Fit", do.call(rbind, checks), list(laplace = fit$laplace))
}

diagnose_vb_fit <- function(fit, rel_obj_warning = 1e-3, rel_obj_problem = 1e-2, ...) {
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
  .new_diagnose_result("VB_Fit", do.call(rbind, checks), list(laplace = fit$laplace))
}

diagnose_classic_fit <- function(fit, ...) {
  checks <- list()
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
  .new_diagnose_result("Classic_Fit", do.call(rbind, checks), list(df_method = fit$df_method))
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
  invisible(x)
}
