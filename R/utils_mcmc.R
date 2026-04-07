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
  test_dat <- read.csv(test_file, header = TRUE)
  n_samples <- nrow(test_dat)

  orig_pl <- model$par_list
  random_flags <- sapply(orig_pl, function(x) isTRUE(x$random))

  if (laplace && any(random_flags)) {
    pl_fixed  <- BayesRTMB:::parse_parameters(orig_pl[!random_flags])
    pl_random <- BayesRTMB:::parse_parameters(orig_pl[random_flags])
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

    dat <- read.csv(file_path, header = TRUE)
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
