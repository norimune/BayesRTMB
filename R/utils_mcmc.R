#' Check MCMC progress from saved files
#'
#' @param name Base name of the saved files.
#' @param dir Directory where progress files are saved. Default is "BayesRTMB_mcmc".
#' @param chains Number of chains to check.
#' @export
check_progress <- function(name, dir = "BayesRTMB_mcmc", chains = 4) {
  cat(sprintf("--- Current MCMC Progress ('%s') ---\n", name))
  for (c in 1:chains) {
    prog_file <- file.path(dir, paste0(name, "_progress_", c, ".txt"))
    if (file.exists(prog_file)) {
      iter <- readLines(prog_file, warn = FALSE)
      cat(sprintf("Chain %d: Iteration %s\n", c, iter[1]))
    } else {
      cat(sprintf("Chain %d: Waiting or file not found\n", c))
    }
  }
  cat("---------------------------------------\n")
}

#' Restore MCMC Fit from CSV
#'
#' @param model An RTMB_Model object.
#' @param name Base name of the saved CSVs.
#' @param dir Directory where CSVs are saved. Default is "BayesRTMB_mcmc".
#' @param chains Number of chains. Default is 4.
#' @param sampling Number of sampling iterations. Default is 1000.
#' @param warmup Number of warmup iterations. Default is 1000.
#' @param thin Thinning interval. Default is 1.
#' @param laplace Logical; whether Laplace approximation was used. Default is FALSE.
#' @return An MCMC_Fit object.
#' @export
read_mcmc_csv <- function(model, name, dir = "BayesRTMB_mcmc", chains = 4,
                          sampling = 1000, warmup = 1000, thin = 1, laplace = FALSE) {

  iter <- sampling + warmup
  mcmc_index <- seq(from = (warmup + 1), to = iter, by = thin)

  orig_pl <- model$par_list
  random_flags <- sapply(orig_pl, function(x) isTRUE(x$random))

  if (laplace && any(random_flags)) {
    pl_fixed  <- BayesRTMB:::parse_parameters(orig_pl[!random_flags])
    pl_random <- BayesRTMB:::parse_parameters(orig_pl[random_flags])
    fixed_idx  <- which(model$pl_full$names %in% pl_fixed$names)
    random_idx <- which(model$pl_full$names %in% pl_random$names)
  } else {
    pl_fixed   <- model$pl_full
    pl_random  <- NULL
    fixed_idx  <- 1:length(model$pl_full$names)
    random_idx <- integer(0)
  }

  P_fixed <- length(pl_fixed$names)
  P_random <- if(!is.null(pl_random)) length(pl_random$names) else 0

  fit <- array(NA, dim=c(length(mcmc_index), chains, P_fixed + 1))
  dimnames(fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = c("lp", pl_fixed$names))

  if (P_random > 0) {
    random_fit <- array(NA, dim=c(length(mcmc_index), chains, P_random))
    dimnames(random_fit) <- list(iteration = NULL, chain = paste0("chain", 1:chains), variable = pl_random$names)
  } else {
    random_fit <- NULL
  }

  accept_mat <- array(NA, dim=c(length(mcmc_index), chains))
  td_mat <- array(NA, dim=c(length(mcmc_index), chains))
  eps_vec <- numeric(chains)

  # 制約空間に戻すためのADオブジェクトを生成
  ad_setup <- model$build_ad_obj(init = model$pl_full$init, laplace = laplace, jacobian_target = "all")
  ad_obj <- ad_setup$ad_obj

  for (c in 1:chains) {
    file_path <- file.path(dir, paste0(name, "-", c, ".csv"))
    if (!file.exists(file_path)) stop(paste("File not found:", file_path))

    dat <- read.csv(file_path, header = TRUE)
    if (nrow(dat) < iter) warning(sprintf("Chain %d has %d iterations, expected %d", c, nrow(dat), iter))

    lp_all <- dat$lp
    accept_all <- dat$accept
    td_all <- dat$treedepth
    eps_all <- dat$eps
    para_unc <- as.matrix(dat[, 6:ncol(dat)])

    para_final <- array(NA, dim = c(length(mcmc_index), length(model$pl_full$names)))

    cat(sprintf("Reconstructing chain %d...\n", c))
    pb <- txtProgressBar(min = 0, max = length(mcmc_index), style = 3)

    # MCMCサンプリング後と同様に、無制約パラメータを制約パラメータに戻す
    for (i in seq_along(mcmc_index)) {
      orig_i <- mcmc_index[i]
      x_in <- para_unc[orig_i, ]

      if (laplace && length(ad_obj$env$random) > 0) {
        ad_obj$fn(x_in)
        para_list <- ad_obj$env$parList()
      } else {
        para_list <- ad_obj$env$parList(x = x_in)
      }

      con_list <- BayesRTMB:::to_constrained(para_list, model$par_list)
      para_final[i, ] <- unlist(con_list, use.names = FALSE)
      setTxtProgressBar(pb, i)
    }
    close(pb)

    fit[, c, 1] <- lp_all[mcmc_index]
    for(j in 1:P_fixed) fit[, c, j+1] <- para_final[, fixed_idx[j]]
    if (P_random > 0) {
      for (j in 1:P_random) random_fit[, c, j] <- para_final[, random_idx[j]]
    }

    accept_mat[, c] <- accept_all[mcmc_index]
    td_mat[, c] <- td_all[mcmc_index]
    eps_vec[c] <- eps_all[nrow(dat)]
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
