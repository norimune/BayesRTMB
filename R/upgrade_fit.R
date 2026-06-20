#' Upgrade a saved fit object to the current BayesRTMB class definitions
#'
#' `upgrade_fit()` rebuilds an existing fit object with the currently loaded
#' BayesRTMB R6 class definitions. This is useful for saved fit objects created
#' by older package versions when newer methods have been added to the fit
#' classes.
#'
#' The function does not refit the model or recompute estimates. It only copies
#' the stored estimation results into a new fit object. When `upgrade_model` is
#' `TRUE`, the embedded `RTMB_Model` object is also rebuilt from the model
#' components stored in the fit, so newer model methods are available as well.
#'
#' @param fit A BayesRTMB fit object. Supported classes are MCMC, VB, MAP, and
#'   classic fit objects.
#' @param upgrade_model Logical; if `TRUE`, rebuild the embedded model object
#'   from its stored data, parameter list, and model functions. If rebuilding
#'   fails, the original model is kept and a warning is issued.
#'
#' @return A fit object of the same type, using the currently loaded class
#'   definitions.
#' @export
upgrade_fit <- function(fit, upgrade_model = TRUE) {
  if (!inherits(fit, "RTMB_Fit_Base") &&
      !inherits(fit, "mcmc_fit") &&
      !inherits(fit, "advi_fit") &&
      !inherits(fit, "vb_fit") &&
      !inherits(fit, "map_fit") &&
      !inherits(fit, "Classic_Fit")) {
    stop("`fit` must be a BayesRTMB fit object.", call. = FALSE)
  }

  model <- .upgrade_fit_model(.upgrade_field(fit, "model"), upgrade_model)

  if (inherits(fit, "mcmc_fit")) {
    out <- .upgrade_mcmc_fit(fit, model)
  } else if (inherits(fit, "advi_fit") || inherits(fit, "vb_fit")) {
    out <- .upgrade_vb_fit(fit, model)
  } else if (inherits(fit, "map_fit")) {
    out <- .upgrade_map_fit(fit, model)
  } else if (inherits(fit, "Classic_Fit")) {
    out <- .upgrade_classic_fit(fit, model)
  } else {
    stop("Unsupported fit class: ", paste(class(fit), collapse = ", "), call. = FALSE)
  }

  .upgrade_restore_classes(out, fit)
}

.upgrade_field <- function(x, name, default = NULL) {
  if (is.null(x)) return(default)
  tryCatch(x[[name]], error = function(e) default)
}

.upgrade_set_field <- function(x, name, value) {
  tryCatch({
    x[[name]] <- value
  }, error = function(e) NULL)
  invisible(x)
}

.upgrade_restore_classes <- function(out, old) {
  class(out) <- unique(c(setdiff(class(old), class(out)), class(out)))
  out
}

.upgrade_fit_model <- function(model, upgrade_model) {
  if (!isTRUE(upgrade_model) || is.null(model) || !inherits(model, "RTMB_Model")) {
    return(model)
  }

  rebuilt <- tryCatch({
    old_silent <- getOption("BayesRTMB.silent", FALSE)
    options(BayesRTMB.silent = TRUE)
    on.exit(options(BayesRTMB.silent = old_silent), add = TRUE)

    new_model <- RTMB_Model$new(
      data = .upgrade_field(model, "data"),
      par_list = .upgrade_field(model, "par_list"),
      log_prob = .upgrade_field(model, "log_prob"),
      transform = .upgrade_field(model, "transform"),
      generate = .upgrade_field(model, "generate"),
      par_names = .upgrade_field(model, "par_names"),
      init = .upgrade_field(model, "init"),
      view = .upgrade_field(model, "view"),
      code = .upgrade_field(model, "code"),
      gr_test = FALSE
    )

    for (field in c(
      "formula", "raw_data", "family", "type", "extra", "contrasts",
      "requested_contrasts", "prior_correction", "fixed_prior_specs", "map"
    )) {
      .upgrade_set_field(new_model, field, .upgrade_field(model, field))
    }

    new_model
  }, error = function(e) {
    warning(
      "Could not upgrade the embedded model; reusing the original model. ",
      "Reason: ", conditionMessage(e),
      call. = FALSE
    )
    model
  })

  rebuilt
}

.upgrade_mcmc_fit <- function(fit, model) {
  out <- MCMC_Fit$new(
    model = model,
    fit = .upgrade_field(fit, "fit"),
    random_fit = .upgrade_field(fit, "random_fit"),
    eps = .upgrade_field(fit, "eps"),
    accept = .upgrade_field(fit, "accept"),
    treedepth = .upgrade_field(fit, "treedepth"),
    laplace = .upgrade_field(fit, "laplace"),
    posterior_mean = .upgrade_field(fit, "posterior_mean"),
    max_treedepth = .upgrade_field(fit, "max_treedepth"),
    pd_error_count = .upgrade_field(fit, "pd_error_count"),
    n_leapfrog = .upgrade_field(fit, "n_leapfrog"),
    divergent = .upgrade_field(fit, "divergent"),
    energy = .upgrade_field(fit, "energy"),
    metric = .upgrade_field(fit, "metric"),
    metric_type = .upgrade_field(fit, "metric_type"),
    metric_init = .upgrade_field(fit, "metric_init"),
    metric_requested = .upgrade_field(fit, "metric_requested"),
    metric_effective = .upgrade_field(fit, "metric_effective"),
    metric_auto = .upgrade_field(fit, "metric_auto"),
    metric_adaptation = .upgrade_field(fit, "metric_adaptation"),
    nuts_variant = .upgrade_field(fit, "nuts_variant"),
    warmup_diagnostics = .upgrade_field(fit, "warmup_diagnostics")
  )

  for (field in c(
    "transform_fit", "transform_dims", "generate_fit", "generate_dims",
    "log_ml", "comparison_fit"
  )) {
    .upgrade_set_field(out, field, .upgrade_field(fit, field))
  }

  out
}

.upgrade_vb_fit <- function(fit, model) {
  out <- VB_Fit$new(
    model = model,
    fit = .upgrade_field(fit, "fit"),
    random_fit = .upgrade_field(fit, "random_fit"),
    elbo_history = .upgrade_field(fit, "elbo_history"),
    laplace = .upgrade_field(fit, "laplace"),
    posterior_mean = .upgrade_field(fit, "posterior_mean"),
    ELBO = .upgrade_field(fit, "ELBO"),
    rel_obj_vals = .upgrade_field(fit, "rel_obj_vals"),
    best_chain = .upgrade_field(fit, "best_chain"),
    mu_history = .upgrade_field(fit, "mu_history")
  )

  for (field in c(
    "transform_fit", "transform_dims", "generate_fit", "generate_dims"
  )) {
    .upgrade_set_field(out, field, .upgrade_field(fit, field))
  }

  out
}

.upgrade_map_fit <- function(fit, model) {
  MAP_Fit$new(
    model = model,
    par_vec = .upgrade_field(fit, "par_vec"),
    par = .upgrade_field(fit, "par"),
    objective = .upgrade_field(fit, "objective"),
    log_ml = .upgrade_field(fit, "log_ml"),
    convergence = .upgrade_field(fit, "convergence"),
    sd_rep = .upgrade_field(fit, "sd_rep"),
    df_fixed = .upgrade_field(fit, "df_fixed"),
    random_effects = .upgrade_field(fit, "random_effects"),
    df_transform = .upgrade_field(fit, "df_transform"),
    df_generate = .upgrade_field(fit, "df_generate"),
    opt_history = .upgrade_field(fit, "opt_history"),
    transform = .upgrade_field(fit, "transform"),
    generate = .upgrade_field(fit, "generate"),
    se_samples = .upgrade_field(fit, "se_samples"),
    par_unc = .upgrade_field(fit, "par_unc"),
    ci_method = .upgrade_field(fit, "ci_method", "wald"),
    laplace = .upgrade_field(fit, "laplace", TRUE),
    map = .upgrade_field(fit, "map"),
    vcov_unc = .upgrade_field(fit, "vcov_unc"),
    marginal_vars = .upgrade_field(fit, "marginal_vars"),
    laplace_random_vars = .upgrade_field(fit, "laplace_random_vars"),
    idx_fix_active = .upgrade_field(fit, "idx_fix_active"),
    show_df = .upgrade_field(fit, "show_df", TRUE),
    view = .upgrade_field(fit, "view"),
    fallback_needed = .upgrade_field(fit, "fallback_needed")
  )
}

.upgrade_classic_fit <- function(fit, model) {
  out <- Classic_Fit$new(
    model = model,
    par_vec = .upgrade_field(fit, "par_vec"),
    par = .upgrade_field(fit, "par"),
    objective = .upgrade_field(fit, "objective"),
    log_lik = .upgrade_field(fit, "log_lik"),
    restricted_log_lik = .upgrade_field(fit, "restricted_log_lik"),
    convergence = .upgrade_field(fit, "convergence"),
    sd_rep = .upgrade_field(fit, "sd_rep"),
    df_fixed = .upgrade_field(fit, "df_fixed"),
    random_effects = .upgrade_field(fit, "random_effects"),
    df_transform = .upgrade_field(fit, "df_transform"),
    df_generate = .upgrade_field(fit, "df_generate"),
    opt_history = .upgrade_field(fit, "opt_history"),
    transform = .upgrade_field(fit, "transform"),
    generate = .upgrade_field(fit, "generate"),
    se_samples = .upgrade_field(fit, "se_samples"),
    par_unc = .upgrade_field(fit, "par_unc"),
    vcov_unc = .upgrade_field(fit, "vcov_unc"),
    ci_method = .upgrade_field(fit, "ci_method", "wald"),
    laplace = .upgrade_field(fit, "laplace", TRUE),
    map = .upgrade_field(fit, "map"),
    test_results = .upgrade_field(fit, "test_results", list()),
    view = .upgrade_field(fit, "view"),
    fit = .upgrade_field(fit, "fit"),
    vcov = .upgrade_field(fit, "vcov"),
    df_method = .upgrade_field(fit, "df_method", "auto"),
    idx_fix_active = .upgrade_field(fit, "idx_fix_active"),
    show_df = .upgrade_field(fit, "show_df", TRUE),
    rss = .upgrade_field(fit, "rss"),
    df_residual = .upgrade_field(fit, "df_residual")
  )

  for (field in c("se_method", "cluster", "robust_type", "bootstrap_results")) {
    .upgrade_set_field(out, field, .upgrade_field(fit, field))
  }

  out
}
