#' MCMC fit object
#'
#' An R6 class storing posterior samples and related information
#' from MCMC estimation.
#'
#' @param model An `RTMB_Model` object used for estimation.
#' @param fit Posterior draws for model parameters.
#' @param random_fit Posterior draws for random effects, if available.
#' @param transform_fit Posterior draws for transformed parameters, if available.
#' @param transform_dims Dimension information for transformed parameters.
#' @param generate_fit Posterior draws for generated quantities, if available.
#' @param generate_dims Dimension information for generated quantities.
#' @param eps Step size used by the sampler.
#' @param accept Acceptance statistics from sampling.
#' @param treedepth Tree depth used in HMC/NUTS sampling.
#' @param laplace Logical; whether Laplace approximation was used.
#' @param posterior_mean Posterior mean estimates.
#' @param pars Names of parameters to extract or summarize.
#' @param chains Chains to include.
#' @param inc_random Logical; whether to include random effects.
#' @param inc_transform Logical; whether to include transformed parameters.
#' @param inc_generate Logical; whether to include generated quantities.
#' @param max_rows Maximum number of rows to print in summaries.
#' @param digits Number of digits to print.
#' @param method Method name.
#' @param seed Random seed.
#' @param tran_fn A function for transformed parameters.
#' @param gq_fn A function for generated quantities.
#' @param target Target variable(s) for rotation or relabeling.
#' @param type Rotation type.
#' @param linked_straight Linked variables transformed in the same direction.
#' @param linked_inverse Linked variables transformed in the opposite direction.
#' @param overwrite Logical; whether to overwrite stored draws.
#' @param linked Linked variables.
#' @param loadings Factor loading variables.
#' @param scores Factor score variables.
#' @param linked_loadings Linked loading variables.
#' @param ... Additional arguments.
#'
#' @field model An `RTMB_Model` object used for estimation.
#' @field fit Posterior draws for model parameters.
#' @field random_fit Posterior draws for random effects.
#' @field transform_fit Posterior draws for transformed parameters.
#' @field transform_dims Dimension information for transformed parameters.
#' @field generate_fit Posterior draws for generated quantities.
#' @field generate_dims Dimension information for generated quantities.
#' @field eps Step size used by the sampler.
#' @field accept Acceptance statistics from sampling.
#' @field treedepth Tree depth used in HMC/NUTS sampling.
#' @field laplace Logical; whether Laplace approximation was used.
#' @field posterior_mean Posterior mean estimates.
#' @field log_ml Numeric value storing the calculated log marginal likelihood from bridge sampling.
#' @field null_fit An \code{MCMC_Fit} object containing the fitted null model. This is automatically cached when calculating a Bayes factor using a target string.
#'
MCMC_Fit <- R6::R6Class(
  classname = "mcmc_fit",
  inherit = RTMB_Fit_Base,

  public = list(
    # --- Fields ---
    model          = NULL, # Reference to RTMB_Model instance
    fit            = NULL,
    random_fit     = NULL,
    transform_fit  = NULL, # Store transformed quantities
    generate_fit   = NULL, # Store generated quantities
    transform_dims = NULL, # Store dimension info for transformed quantities
    generate_dims  = NULL, # Store dimension info for generated quantities
    eps            = NULL,
    accept         = NULL,
    treedepth      = NULL,
    laplace        = NULL,
    posterior_mean = NULL,
    log_ml          = NULL,
    null_fit       = NULL,

    # 1. Constructor
    #' @description Get point estimate for a target parameter (internal use).
    #' @param target Target parameter name.
    #' @return Matrix or array of point estimate.
    get_point_estimate = function(target) {
      target_draws <- self$draws(pars = target, inc_transform = TRUE, inc_generate = TRUE)
      if (dim(target_draws)[3] == 0) stop("Parameter not found: ", target)

      lp_draws <- self$draws(pars = "lp", inc_transform = FALSE, inc_generate = FALSE)
      max_idx <- which(lp_draws == max(lp_draws, na.rm = TRUE), arr.ind = TRUE)
      best_iter <- max_idx[1, 1]
      best_chain <- max_idx[1, 2]

      t_info <- self$model$par_list[[target]]
      if (!is.null(t_info)) {
        t_dim <- t_info$dim
      } else if (!is.null(self$transform_dims[[target]])) {
        t_dim <- self$transform_dims[[target]]
      } else if (!is.null(self$generate_dims[[target]])) {
        t_dim <- self$generate_dims[[target]]
      } else {
        t_dim <- dim(target_draws)[3]
      }

      target_map_flat <- target_draws[best_iter, best_chain, ]
      target_map <- target_map_flat
      if (length(t_dim) > 1) dim(target_map) <- t_dim

      return(target_map)
    },

    #' @description Create a new `MCMC_Fit` object.
    initialize = function(model, fit, random_fit, eps, accept, treedepth, laplace, posterior_mean) {
      self$model <- model
      self$fit <- fit
      self$random_fit <- random_fit
      self$eps <- eps
      self$accept <- accept
      self$treedepth <- treedepth
      self$laplace <- laplace
      self$posterior_mean <- posterior_mean
      self$transform_fit <- NULL
      self$transform_dims <- list()
      self$generate_fit <- NULL
      self$generate_dims <- list()
      # Ensure S3 dispatch works for base class methods
      class(self) <- c(class(self), "RTMB_Fit_Base")
    },

    #' @description Print a brief summary of the fitted object.
    #' @return The object itself, invisibly.
    print = function(...) {
      print(self$summary(...))
      invisible(self)
    },

    #' @description Extract posterior draws for selected parameters.
    #' @param pars Character or numeric vector specifying the names or indices of parameters to extract. If NULL, all available parameters are extracted.
    #' @param chains Numeric vector specifying the chains to extract. If NULL, draws from all chains are returned.
    #' @param best_chains Integer; number of best chains to retain based on mean log-posterior (lp).
    #' @param inc_random Logical; whether to include random effects in the output. Default is FALSE.
    #' @param inc_transform Logical; whether to include transformed parameters in the output. Default is TRUE.
    #' @param inc_generate Logical; whether to include generated quantities in the output. Default is TRUE.
    #' @return Posterior draws.
    draws = function(pars = NULL, chains = NULL, best_chains = NULL,
                     inc_random = FALSE, inc_transform = TRUE, inc_generate = TRUE) {

      # 1. Determine chains to use
      total_chains <- dim(self$fit)[2]
      available_chains <- 1:total_chains

      if (!is.null(chains)) {
        available_chains <- intersect(available_chains, chains)
        if (length(available_chains) == 0) {
          stop("The specified chains were not found.", call. = FALSE)
        }
      }

      if (!is.null(best_chains)) {
        lp_idx <- which(dimnames(self$fit)[[3]] == "lp")
        if (length(lp_idx) == 0) lp_idx <- 1

        # Calculate mean of lp for available_chains
        lp_means <- apply(self$fit[, available_chains, lp_idx, drop = FALSE], 2, mean, na.rm = TRUE)
        n_best <- min(best_chains, length(available_chains))
        ordered_idx <- order(lp_means, decreasing = TRUE)
        selected_rel_idx <- ordered_idx[1:n_best]
        available_chains <- available_chains[selected_rel_idx]
      }

      available_chains <- sort(available_chains)

      # 2. Extract and bind arrays
      out_array <- self$fit[, available_chains, , drop = FALSE]

      if (inc_random && !is.null(self$random_fit)) {
        P1 <- dim(out_array)[3]; P2 <- dim(self$random_fit)[3]
        I <- dim(out_array)[1]; C <- length(available_chains)
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$random_fit[, available_chains, , drop = FALSE]
        dimnames(new_out) <- list(
          iteration = dimnames(out_array)[[1]],
          chain = dimnames(out_array)[[2]],
          variable = c(dimnames(out_array)[[3]], dimnames(self$random_fit)[[3]])
        )
        out_array <- new_out
      }

      if (inc_transform && !is.null(self$transform_fit)) {
        P1 <- dim(out_array)[3]; P2 <- dim(self$transform_fit)[3]
        I <- dim(out_array)[1]; C <- length(available_chains)
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$transform_fit[, available_chains, , drop = FALSE]
        dimnames(new_out) <- list(
          iteration = dimnames(out_array)[[1]],
          chain = dimnames(out_array)[[2]],
          variable = c(dimnames(out_array)[[3]], dimnames(self$transform_fit)[[3]])
        )
        out_array <- new_out
      }

      if (inc_generate && !is.null(self$generate_fit)) {
        P1 <- dim(out_array)[3]; P2 <- dim(self$generate_fit)[3]
        I <- dim(out_array)[1]; C <- length(available_chains)
        new_out <- array(NA, dim = c(I, C, P1 + P2))
        new_out[,,1:P1] <- out_array
        new_out[,,(P1+1):(P1+P2)] <- self$generate_fit[, available_chains, , drop = FALSE]
        dimnames(new_out) <- list(
          iteration = dimnames(out_array)[[1]],
          chain = dimnames(out_array)[[2]],
          variable = c(dimnames(out_array)[[3]], dimnames(self$generate_fit)[[3]])
        )
        out_array <- new_out
      }

      P <- dim(out_array)[3]
      param_names <- dimnames(out_array)[[3]]
      if (is.null(param_names)) param_names <- paste0("V", 1:P)

      target_idx <- 1:P

      if (!is.null(pars)) {
        if (is.numeric(pars)) {
          valid_idx <- pars[pars >= 1 & pars <= P]
          if (length(valid_idx) == 0) {
            stop("The index specified in 'pars' was not found.", call. = FALSE)
          }
          target_idx <- valid_idx

        } else if (is.character(pars)) {
          base_names <- gsub("\\[.*\\]$", "", param_names)
          matched <- which(param_names %in% pars | base_names %in% pars)
          if (length(matched) == 0) {
            stop("The variable name specified in 'pars' was not found.", call. = FALSE)
          }
          target_idx <- matched

        } else {
          stop("'pars' must be either numeric or character.", call. = FALSE)
        }
      }

      return(out_array[, , target_idx, drop = FALSE])
    },

    #' @description Summarize posterior draws.
    #' @param pars Character or numeric vector specifying the names or indices of parameters to summarize. If NULL, all available parameters are summarized.
    #' @param chains Numeric vector specifying the chains to extract. If NULL, draws from all chains are used.
    #' @param best_chains Integer; number of best chains to retain based on mean log-posterior (lp).
    #' @param max_rows Integer; maximum number of rows to print in the summary table. Default is 10.
    #' @param digits Integer; number of decimal places to print. Default is 2.
    #' @param inc_random Logical; whether to include random effects in the summary. Default is FALSE.
    #' @param inc_transform Logical; whether to include transformed parameters in the summary. Default is TRUE.
    #' @param inc_generate Logical; whether to include generated quantities in the summary. Default is TRUE.
    #' @return A summary object.
    summary = function(pars = NULL, chains = NULL, best_chains = NULL,
                       max_rows = 10, digits = 2,
                       inc_random = FALSE, inc_transform = TRUE, inc_generate = TRUE){
      draws_array <- self$draws(pars = pars,
                                chains = chains,
                                best_chains = best_chains,
                                inc_random = inc_random,
                                inc_transform = inc_transform,
                                inc_generate = inc_generate)

      P <- dim(draws_array)[3]
      param_names <- dimnames(draws_array)[[3]]

      target_idx <- 1:P

      # --- Priority sorting by lp and model$view ---
      if (length(target_idx) > 0) {
        current_names <- param_names[target_idx]
        base_names <- gsub("\\[.*\\]$", "", current_names)

        # Always prioritize lp
        target_views <- c("lp")
        if (!is.null(self$model$view)) {
          target_views <- c(target_views, self$model$view)
        }

        priority_sub_idx <- integer(0)
        for (v in target_views) {
          match_idx <- which(current_names == v | base_names == v)
          priority_sub_idx <- c(priority_sub_idx, match_idx)
        }
        priority_sub_idx <- unique(priority_sub_idx)
        other_sub_idx <- setdiff(seq_along(current_names), priority_sub_idx)
        target_idx <- target_idx[c(priority_sub_idx, other_sub_idx)]
      }

      if (!is.null(max_rows)) {
        limit <- min(length(target_idx), as.integer(max_rows))
        target_idx <- target_idx[1:limit]
      }

      res_list_sum <- vector("list", length(target_idx))

      for (i in seq_along(target_idx)) {
        p <- target_idx[i]
        mat_p <- as.matrix(draws_array[, , p])
        vec_p <- as.vector(mat_p)
        valid_vec <- vec_p[is.finite(vec_p)]

        if (length(valid_vec) == 0) {
          res_list_sum[[i]] <-
            data.frame(
              variable = param_names[p],
              mean = NA,
              sd = NA,
              map = NA,
              q2.5 = NA,
              q97.5 = NA,
              ess_bulk = NA,
              ess_tail = NA,
              rhat = NA,
              stringsAsFactors = FALSE
            )
          next
        }

        sd_val <- sd(valid_vec)
        if (is.na(sd_val) || sd_val < 1e-10) {
          map_val <- valid_vec[1]; q95 <- c(valid_vec[1], valid_vec[1])
          rhat_val <- NA; ebulk_val <- NA; etail_val <- NA
        } else {
          map_val   <- map_est(valid_vec)
          q95       <- quantile95(valid_vec)
          rhat_val  <- r_hat(mat_p)
          ebulk_val <- ess_bulk(mat_p)
          etail_val <- ess_tail95(mat_p)
        }

        res_list_sum[[i]] <- data.frame(
          variable = param_names[p],
          mean     = mean(valid_vec),
          sd       = sd_val,
          map      = map_val,
          q2.5     = unname(q95[1]),
          q97.5    = unname(q95[2]),
          ess_bulk = ebulk_val,
          ess_tail = etail_val,
          rhat     = rhat_val,
          stringsAsFactors = FALSE
        )
      }
      res_df <- do.call(rbind, res_list_sum)
      # Replace extremely small values close to zero with exact 0 in numeric columns of the dataframe
      num_cols <- sapply(res_df, is.numeric)
      res_df[num_cols] <- lapply(res_df[num_cols], function(x) {
        x[abs(x) < 1e-12 & !is.na(x)] <- 0
        return(x)
      })
      class(res_df) <- c("summary_BayesRTMB", "data.frame")
      attr(res_df, "digits") <- digits
      return(res_df)
    },

    #' @description Transform posterior draws to the unconstrained scale.
    #' @return Posterior draws on the unconstrained scale.
    unconstrain_draws = function() {
      I <- dim(self$fit)[1]
      C <- dim(self$fit)[2]
      orig_pl <- self$model$par_list
      random_flags <- sapply(orig_pl, function(x) isTRUE(x$random))

      target_par_list <- if (self$laplace && any(random_flags)) orig_pl[!random_flags] else orig_pl
      draws_ob <- self$fit[, , -1, drop = FALSE]

      # --- Modified: Calculate the number of effective parameters considering map ---
      map_list <- self$model$map
      if (!is.null(map_list)) {
        Q <- 0
        for (name in names(target_par_list)) {
          L <- target_par_list[[name]]$unc_length
          if (L > 0) {
            if (!is.null(map_list[[name]])) {
              Q <- Q + sum(!is.na(map_list[[name]]))
            } else {
              Q <- Q + L
            }
          }
        }
      } else {
        Q <- sum(sapply(target_par_list, function(x) x$unc_length))
      }

      draws_uc <- array(NA, dim = c(I, C, Q))

      for (i in 1:I) {
        for (j in 1:C) {
          con_vec <- draws_ob[i, j, ]
          con_list <- constrained_vector_to_list(con_vec, target_par_list)
          unc_list <- to_unconstrained(con_list, target_par_list)

          # --- Modified: Extract only non-fixed parameters if map is specified ---
          if (!is.null(map_list)) {
            unc_vec <- c()
            for (name in names(target_par_list)) {
              val <- unc_list[[name]]
              if (length(val) > 0) {
                if (!is.null(map_list[[name]])) {
                  unc_vec <- c(unc_vec, val[!is.na(map_list[[name]])])
                } else {
                  unc_vec <- c(unc_vec, val)
                }
              }
            }
            draws_uc[i, j, ] <- unc_vec
          } else {
            draws_uc[i, j, ] <- unlist(unc_list, use.names = FALSE)
          }
        }
      }
      return(matrix(draws_uc, nrow = I * C, ncol = Q))
    },

    #' @description Evaluate log-probability values.
    #' @return Numeric vector of log-probability values.
    log_prob = function() {
      ad_setup <- self$model$build_ad_obj(init = self$posterior_mean, laplace = self$laplace, jacobian_target = "all")
      ad_obj <- ad_setup$ad_obj

      fn <- function(q) {
        val <- tryCatch({ -ad_obj$fn(q) }, error = function(e) -Inf)
        if (is.na(val) || is.nan(val)) return(-Inf)
        return(val)
      }
      return(fn)
    },

    #' @description Estimate the marginal likelihood by bridge sampling.
    #' @param method Character; the method to use for bridge sampling (e.g., "warp3", "normal"). Default is "warp3".
    #' @param use_neff Logical; whether to use the effective sample size (ESS) to adjust for autocorrelation. Default is TRUE.
    #' @param seed Integer; random seed for reproducibility. Default is NULL.
    #' @param max_iter Integer; maximum number of iterations for the estimation algorithm. Default is 100.
    #' @return Bridge sampling result.
    bridgesampling = function(method = "normal", use_neff = TRUE, seed = NULL, max_iter = 100) {

      if (!is.null(seed)) set.seed(seed)

      draws_uc <- self$unconstrain_draws()
      log_prob_fn <- self$log_prob()

      N_total <- nrow(draws_uc)
      Q <- ncol(draws_uc)

      index1 <- seq(1, N_total - 1, by = 2)
      index2 <- seq(2, N_total, by = 2)
      z_post <- draws_uc[index1, , drop = FALSE]
      z_fit  <- draws_uc[index2, , drop = FALSE]

      M1 <- nrow(z_post)
      M2 <- nrow(z_fit)
      M  <- M1 + M2

      # Pre-computation: parameter estimation of proposal distribution
      meanz <- apply(z_fit, 2, mean)
      covz <- cov(z_fit)

      # Calculate log-posterior for posterior samples
      lp_post <- apply(z_post, 1, log_prob_fn)

      if (isTRUE(use_neff)) {
        n_eff <- tryCatch({
          posterior::ess_basic(lp_post)
        }, error = function(e) M1)
        if (is.na(n_eff) || n_eff <= 0) n_eff <- M1
      } else {
        n_eff <- M1
      }

      # Optimal weights based on Meng & Wong
      S1 <- n_eff / (n_eff + M2)
      S2 <- M2 / (n_eff + M2)

      if (method == "normal") {
        z_propose <- MASS::mvrnorm(M2, meanz, covz)
        log_propose <- function(z) mvtnorm::dmvnorm(z, meanz, covz, log = TRUE)

        # Reuse pre-calculated lp_post
        log_L1 <- lp_post - log_propose(z_post)
        log_L2 <- apply(z_propose, 1, log_prob_fn) - log_propose(z_propose)

        valid_log_L1 <- log_L1[is.finite(log_L1)]
        if (length(valid_log_L1) == 0) stop("No valid posterior samples found.")
        log_Lm <- median(valid_log_L1)

        L1 <- exp(log_L1 - log_Lm)
        L2 <- exp(log_L2 - log_Lm)
        L1[is.na(L1) | is.nan(L1)] <- 0
        L2[is.na(L2) | is.nan(L2)] <- 0

        ml_t <- exp(log_prob_fn(meanz) - log_Lm)
        if (is.na(ml_t) || is.nan(ml_t) || ml_t == 0) ml_t <- 1e-10

      } else if (method == "warp3") {
        L_chol <- t(chol(covz + diag(0.0001, Q)))
        log_det_L <- sum(log(diag(L_chol)))

        eta_post <- t(solve(L_chol, t(sweep(z_post, 2, meanz, "-"))))
        eta_prop <- matrix(rnorm(M2 * Q), nrow = M2, ncol = Q)

        log_propose_std <- function(eta) {
          sum(dnorm(eta, mean = 0, sd = 1, log = TRUE))
        }

        log_target_warp3 <- function(eta) {
          z_plus  <- meanz + as.vector(L_chol %*% eta)
          z_minus <- meanz - as.vector(L_chol %*% eta)

          lp_plus  <- log_prob_fn(z_plus)
          lp_minus <- log_prob_fn(z_minus)

          max_lp <- max(lp_plus, lp_minus)
          if (is.infinite(max_lp)) return(-Inf)

          log_mix <-
            log(0.5) + max_lp + log(exp(lp_plus - max_lp) + exp(lp_minus - max_lp))
          return(log_mix + log_det_L)
        }

        log_L1 <- apply(eta_post, 1, log_target_warp3) - apply(eta_post, 1, log_propose_std)
        log_L2 <- apply(eta_prop, 1, log_target_warp3) - apply(eta_prop, 1, log_propose_std)

        valid_log_L1 <- log_L1[is.finite(log_L1)]
        if (length(valid_log_L1) == 0) stop("No valid posterior samples found.")
        log_Lm <- median(valid_log_L1)

        L1 <- exp(log_L1 - log_Lm)
        L2 <- exp(log_L2 - log_Lm)
        L1[is.na(L1) | is.nan(L1)] <- 0
        L2[is.na(L2) | is.nan(L2)] <- 0

        ml_t <- exp(log_target_warp3(rep(0, Q)) - log_Lm)
        if (is.na(ml_t) || is.nan(ml_t) || ml_t == 0) ml_t <- 1e-10

      } else {
        stop("The 'method' argument must be either 'normal' or 'warp3'.")
      }

      Trial <- max_iter
      for(t in 1:Trial){
        ml_t_old <- ml_t
        bunbo <- 0
        bunshi <- 0
        for(m in 1:M1) bunbo <- bunbo + 1/(S1*L1[m] + S2*ml_t)
        for(m in 1:M2) bunshi <- bunshi + L2[m] / (S1*L2[m] + S2*ml_t)

        ml_t <- (bunshi/M2) / (bunbo/M1)

        if (is.na(ml_t) || is.nan(ml_t) || ml_t <= 0) {
          warning("ml_t became 0 or NA during calculation.")
          break
        }
        if(abs(log(ml_t) - log(ml_t_old)) < 0.000001){
          break
        }
      }

      logml.bs <- as.numeric(log(ml_t) + log_Lm)

      correction <- self$model$prior_correction
      if (!is.null(correction) && correction != 0) {
        logml.bs <- logml.bs - correction
      }

      # 2. Calculation of approximate standard error
      res <- logml.bs

      if (method == "normal") {
        f1_vec <- ml_t / (S1 * L1 + S2 * ml_t)
        f2_vec <- L2 / (S1 * L2 + S2 * ml_t)

        V1 <- var(f1_vec)
        V2 <- var(f2_vec)
        E1 <- mean(f1_vec)
        E2 <- mean(f2_vec)

        re2 <- V1 / (n_eff * E1^2) + V2 / (M2 * E2^2)
        error_logml <- sqrt(re2)

        cat(sprintf("Bridge Sampling Converged: LogML = %.3f (Error = %.4f, ESS = %.1f)\n", logml.bs, error_logml, n_eff))
        attr(res, "error") <- error_logml

      } else {
        # Return NA if error approximation is difficult, such as in warp3
        cat(sprintf("Bridge Sampling Converged: LogML = %.3f (ESS = %.1f)\n", logml.bs, n_eff))
        attr(res, "error") <- NA_real_
      }

      attr(res, "ess") <- n_eff
      self$log_ml <- res

      return(res)
    },

    #' @description Calculate the Bayes Factor against a null model or another fit object.
    #' @param null_model Either a character string specifying the null target (e.g., "rho ~ uniform(-1, 1)") or another MCMC_Fit object.
    #' @param bs_method Character; the method to use for bridge sampling ("normal" or "warp3"). Default is "normal".
    #' @param error_threshold Numeric; threshold for the approximate error warning. Default is 0.2.
    #' @param ... Additional arguments passed to the sample() method when fitting a null model (e.g., \code{chains = 4}, \code{sampling = 4000}).
    bayes_factor = function(null_model, bs_method = "normal", error_threshold = 0.2, ...) {

      # 1. Calculate and assign marginal likelihood if not already calculated
      if (is.null(self$log_ml)) {
        cat("Calculating marginal likelihood for the full model...\n")
        self$log_ml <- self$bridgesampling(method = bs_method)
      }
      log_ml1 <- self$log_ml

      log_ml2 <- NULL

      # 2. Branch by the type of comparison target (null_model)
      if (is.character(null_model) && length(null_model) == 1) {
        cat(sprintf("\n--- Preparing and Sampling Null Model (%s) ---\n", null_model))

        mdl_null <- self$model$null_model(target = null_model)

        # --- Modified: Inherit full model settings ---
        sample_args <- list(...) # List of arguments specified by the user

        # If not specified, obtain iterations and chains from own fit array dimensions to supplement
        if (is.null(sample_args$sampling)) {
          sample_args$sampling <- dim(self$fit)[1]
        }
        if (is.null(sample_args$chains)) {
          sample_args$chains <- dim(self$fit)[2]
        }

        # Sample null model using supplemented argument list
        fit_null <- do.call(mdl_null$sample, sample_args)

        cat("\n--- Calculating marginal likelihood for the null model ---\n")
        fit_null$log_ml <- fit_null$bridgesampling(method = bs_method)
        log_ml2 <- fit_null$log_ml

        self$null_fit <- fit_null

      } else if (inherits(null_model, "R6") && "bridgesampling" %in% names(null_model)) {
        if (is.null(null_model$log_ml)) {
          cat("Calculating marginal likelihood for the comparison model...\n")
          null_model$log_ml <- null_model$bridgesampling(method = bs_method)
        }
        log_ml2 <- null_model$log_ml

      } else {
        stop("The 'null_model' argument must be a null_model string or another MCMC_Fit object.")
      }

      # 3. Calculation of Bayes factor and evaluation of error
      val1 <- as.numeric(log_ml1)
      val2 <- as.numeric(log_ml2)
      log_bf <- val1 - val2
      bf <- exp(log_bf)

      err1 <- attr(log_ml1, "error")
      err2 <- attr(log_ml2, "error")
      has_error <- !is.null(err1) && !is.null(err2) && !is.na(err1) && !is.na(err2)

      if (has_error) {
        log_bf_err <- sqrt(err1^2 + err2^2)
        if (log_bf_err > error_threshold) {
          warning(
            sprintf(
              "The estimation error of the log Bayes factor (%.3f) exceeds the threshold (%g). Interpretation of the results may be unstable.\n",
              log_bf_err, error_threshold
            ),
            "Consider increasing the number of MCMC samples or ESS to improve precision.\n",
            call. = FALSE, immediate. = TRUE
          )
        }
      } else {
        log_bf_err <- NA_real_
      }

      # 4. Assignment of interpretation
      if (bf > 100) evidence <- "Decisive evidence for Model 1"
      else if (bf > 10) evidence <- "Strong evidence for Model 1"
      else if (bf > 3) evidence <- "Substantial evidence for Model 1"
      else if (bf > 1) evidence <- "Anecdotal evidence for Model 1"
      else if (bf == 1) evidence <- "No evidence"
      else if (bf >= 1/3) evidence <- "Anecdotal evidence for Model 2"
      else if (bf >= 1/10) evidence <- "Substantial evidence for Model 2"
      else if (bf >= 1/100) evidence <- "Strong evidence for Model 2"
      else evidence <- "Decisive evidence for Model 2"

      res <- list(
        BF12 = bf,
        log_BF12 = log_bf,
        log_BF_error = log_bf_err,
        interpretation = evidence
      )

      class(res) <- "bayes_factor"
      return(res)
    },


    #' @description Compute transformed parameters from posterior draws.
    #' @return Transformed parameter draws.
    transformed_draws = function(tran_fn = NULL) {
      all_draws <- self$draws(
        inc_random = TRUE,
        inc_transform = FALSE,
        inc_generate = FALSE
      )

      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2]

      if (is.null(tran_fn)) tran_fn <- self$model$transform

      wrapper_tran_fn <- function(dat, param) {
        res <- list()

        if (!is.null(tran_fn)) {
          user_res <- tran_fn(dat, param)

          if (is.null(user_res)) {
            user_res <- list()
          } else if (!is.list(user_res)) {
            stop("'transformed_parameters' must return a list.", call. = FALSE)
          }

          dup_names <- intersect(names(res), names(user_res))
          if (length(dup_names) > 0) {
            stop(
              sprintf(
                "Transformed parameter names are duplicated: %s",
                paste(dup_names, collapse = ", ")
              ),
              call. = FALSE
            )
          }

          res <- c(res, user_res)
        }

        return(res)
      }

      test_para <- constrained_vector_to_list(all_draws[1, 1, -1], self$model$par_list)
      test_tran <- wrapper_tran_fn(self$model$data, test_para)

      if (length(test_tran) == 0) return(invisible(self))

      tran_names <- character(0)
      self$transform_dims <- list()

      for (name in names(test_tran)) {
        val <- test_tran[[name]]
        len <- length(val)
        dim_val <- dim(val)
        if (is.null(dim_val)) dim_val <- len
        self$transform_dims[[name]] <- dim_val

        names_def <- self$model$par_names[[name]]
        flat_nms <- generate_flat_names(name, dim_val, names_def)
        tran_names <- c(tran_names, flat_nms)
      }

      tran_array <- array(NA, dim = c(iter, chains, length(tran_names)))
      dimnames(tran_array) <- list(
        iteration = NULL,
        chain = paste0("chain", seq_len(chains)),
        variable = tran_names
      )

      cat("Calculating transformed parameters...\n")

      total_steps <- iter * chains
      pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
      counter <- 0

      update_interval <- max(1, floor(total_steps / 100))

      for (c in seq_len(chains)) {
        for (i in seq_len(iter)) {
          p_list <- constrained_vector_to_list(all_draws[i, c, -1], self$model$par_list)
          res <- wrapper_tran_fn(self$model$data, p_list)
          tran_array[i, c, ] <- unlist(res, use.names = FALSE)

          counter <- counter + 1
          if (counter %% update_interval == 0) {
            setTxtProgressBar(pb, counter)
          }
        }
      }
      setTxtProgressBar(pb, total_steps)
      close(pb)

      self$transform_fit <- tran_array
      return(invisible(self))
    },

    #' @description Compute generated quantities from posterior draws.
    #' @param code An `rtmb_code({ ... })` or `{ ... }` block containing the logic
    #' to be calculated using posterior samples.
    #' @return The `MCMC_Fit` object itself (invisibly).
    #' Results are appended to the `generate_fit` field.
    generated_quantities = function(code) {
      raw_code <- substitute(code)

      if (is.name(raw_code)) {
        evaluated <- tryCatch(eval(raw_code, envir = parent.frame()), error = function(e) NULL)
        if (is.language(evaluated) || is.call(evaluated)) {
          code <- evaluated
          raw_code <- evaluated
        }
      }

      if (is.call(raw_code) && identical(raw_code[[1]], as.name("rtmb_code"))) {
        parsed_code <- eval(raw_code, envir = parent.frame())
        if (!"generate" %in% names(parsed_code)) {
          stop("There is no 'generate' block in rtmb_code().")
        }
        gen_ast <- parsed_code$generate
      } else if (is.call(raw_code) && identical(raw_code[[1]], as.name("{"))) {
        gen_ast <- raw_code
      } else if (is.call(code) && identical(code[[1]], as.name("{"))) {
        gen_ast <- code
      } else if (is.list(code) && "generate" %in% names(code)) {
        gen_ast <- code$generate
      } else {
        stop("'code' must be specified in the format rtmb_code(generate = { ... }) or { ... }.")
      }

      cat("Running generated_quantities...\n")

      gen_fn <- eval(bquote(transform_code(.(gen_ast))))
      environment(gen_fn) <- parent.env(globalenv())

      all_draws <- self$draws(
        inc_random = TRUE,
        inc_transform = TRUE,
        inc_generate = FALSE
      )

      iter   <- dim(all_draws)[1]
      chains <- dim(all_draws)[2]

      test_p_list <- constrained_vector_to_list(all_draws[1, 1, -1], self$model$par_list)
      if (!is.null(self$model$transform)) {
        test_tran <- self$model$transform(self$model$data, test_p_list)
        if (is.list(test_tran)) test_p_list <- c(test_p_list, test_tran)
      }

      if (!is.null(self$generate_fit)) {
        existing_gq <- dimnames(self$generate_fit)[[3]]
        for (g_name in names(self$generate_dims)) {
          g_dim <- self$generate_dims[[g_name]]
          flat_nms <- generate_flat_names(g_name, g_dim, self$model$par_names[[g_name]])
          if (all(flat_nms %in% existing_gq)) {
            val <- self$generate_fit[1, 1, flat_nms]
            if (length(g_dim) > 1) dim(val) <- g_dim
            test_p_list[[g_name]] <- val
          }
        }
      }

      test_gq <- gen_fn(self$model$data, test_p_list)

      if (is.null(test_gq) || length(test_gq) == 0) {
        cat("No generated quantities returned.\n")
        return(invisible(self))
      }

      gq_names <- character(0)
      if (is.null(self$generate_dims)) self$generate_dims <- list()

      for (name in names(test_gq)) {
        val <- test_gq[[name]]
        dim_val <- dim(val)
        if (is.null(dim_val)) dim_val <- length(val)
        self$generate_dims[[name]] <- dim_val

        names_def <- self$model$par_names[[name]]
        flat_nms <- generate_flat_names(name, dim_val, names_def)
        gq_names <- c(gq_names, flat_nms)
      }

      # (前略)
      new_gq_array <- array(NA, dim = c(iter, chains, length(gq_names)))
      dimnames(new_gq_array) <- list(
        iteration = NULL,
        chain = paste0("chain", seq_len(chains)),
        variable = gq_names
      )

      total_steps <- iter * chains
      pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
      counter <- 0

      update_interval <- max(1, floor(total_steps / 100))

      for (c in seq_len(chains)) {
        for (i in seq_len(iter)) {
          p_list <- constrained_vector_to_list(all_draws[i, c, -1], self$model$par_list)
          if (!is.null(self$model$transform)) {
            tran_res <- self$model$transform(self$model$data, p_list)
            if (is.list(tran_res)) p_list <- c(p_list, tran_res)
          }
          if (!is.null(self$generate_fit)) {
            existing_gq <- dimnames(self$generate_fit)[[3]]
            for (g_name in names(self$generate_dims)) {
              g_dim <- self$generate_dims[[g_name]]
              flat_nms <- generate_flat_names(g_name, g_dim, self$model$par_names[[g_name]])
              if (all(flat_nms %in% existing_gq)) {
                val <- self$generate_fit[i, c, flat_nms]
                if (length(g_dim) > 1) dim(val) <- g_dim
                p_list[[g_name]] <- val
              }
            }
          }
          res <- gen_fn(self$model$data, p_list)
          new_gq_array[i, c, ] <- unlist(res, use.names = FALSE)

          counter <- counter + 1

          if (counter %% update_interval == 0) {
            setTxtProgressBar(pb, counter)
          }
        }
      }
      setTxtProgressBar(pb, total_steps)
      close(pb)

      if (is.null(self$generate_fit)) {
        self$generate_fit <- new_gq_array
      } else {
        old_gq <- self$generate_fit
        I <- dim(old_gq)[1]; C <- dim(old_gq)[2]
        P1 <- dim(old_gq)[3]; P2 <- dim(new_gq_array)[3]
        merged_gq <- array(NA, dim = c(I, C, P1 + P2))
        merged_gq[,,1:P1] <- old_gq
        merged_gq[,,(P1+1):(P1+P2)] <- new_gq_array
        dimnames(merged_gq) <- list(dimnames(old_gq)[[1]], dimnames(old_gq)[[2]], c(dimnames(old_gq)[[3]], dimnames(new_gq_array)[[3]]))
        self$generate_fit <- merged_gq
      }

      cat("Generated quantities added to samples.\n")
      invisible(self)
    },

    #' @description Resolve label switching in posterior draws.
    #' @param target Character string specifying the target variable to base the relabeling on.
    #' @param linked Character vector of variable names to be relabeled in the same order as the target. Default is NULL.
    #' @param overwrite Logical; whether to overwrite the stored draws in the current object. Default is TRUE.
    #' @param scalar_fns A named list of functions to apply to scalar variables for relabeling. Default is an empty list.
    #' @return Relabeled draws or updated object.
    resolve_switching = function(target, linked = NULL, overwrite = TRUE, scalar_fns = list()) {
      cat(sprintf("Resolving label switching based on '%s'...\n", target))

      obj <- if (isTRUE(overwrite)) self else self$clone(deep = TRUE)

      f_arr <- obj$fit
      r_arr <- obj$random_fit
      t_arr <- obj$transform_fit
      g_arr <- obj$generate_fit

      v_names_f <- dimnames(f_arr)[[3]]
      v_names_r <- if (!is.null(r_arr)) dimnames(r_arr)[[3]] else character(0)
      v_names_t <- if (!is.null(t_arr)) dimnames(t_arr)[[3]] else character(0)
      v_names_g <- if (!is.null(g_arr)) dimnames(g_arr)[[3]] else character(0)

      get_var_info <- function(vname) {
        pattern <- paste0("^", vname, "\\[")

        idx_f <- grep(pattern, v_names_f)
        if (length(idx_f) > 0) return(list(loc = "fixed", idx = idx_f, dim = obj$model$par_list[[vname]]$dim))
        idx_r <- grep(pattern, v_names_r)
        if (length(idx_r) > 0) return(list(loc = "random", idx = idx_r, dim = obj$model$par_list[[vname]]$dim))
        idx_t <- grep(pattern, v_names_t)
        if (length(idx_t) > 0) return(list(loc = "tran", idx = idx_t, dim = obj$transform_dims[[vname]]))
        idx_g <- grep(pattern, v_names_g)
        if (length(idx_g) > 0) return(list(loc = "gq", idx = idx_g, dim = obj$generate_dims[[vname]]))

        idx_f0 <- which(v_names_f == vname)
        if (length(idx_f0) == 1) return(list(loc = "fixed_scalar", idx = idx_f0, dim = 1))
        idx_r0 <- which(v_names_r == vname)
        if (length(idx_r0) == 1) return(list(loc = "random_scalar", idx = idx_r0, dim = 1))
        idx_t0 <- which(v_names_t == vname)
        if (length(idx_t0) == 1) return(list(loc = "tran_scalar", idx = idx_t0, dim = 1))
        idx_g0 <- which(v_names_g == vname)
        if (length(idx_g0) == 1) return(list(loc = "gq_scalar", idx = idx_g0, dim = 1))

        stop(paste0("Label switching failed: Variable '", vname, "' not found."))
      }

      t_info <- get_var_info(target)

      if (grepl("scalar", t_info$loc)) stop(paste0("Target variable '", target, "' must be a vector, not a scalar."))
      if (length(t_info$dim) != 1) stop(paste0("Target variable '", target, "' must be a vector."))

      K <- t_info$dim[1]
      if (length(t_info$idx) != K) stop(paste0("Target variable '", target, "' has inconsistent dimension information."))

      linked_info_list <- list()
      if (!is.null(linked)) {
        for (lvar in linked) {
          l_info <- get_var_info(lvar)

          if (grepl("scalar", l_info$loc)) {
            linked_info_list[[lvar]] <- l_info
          } else {
            if (length(l_info$dim) != 1) {
              warning(sprintf("Variable '%s' is not a vector. Skipping.", lvar))
            } else if (l_info$dim[1] != K) {
              warning(sprintf("Length of '%s' (%d) != target '%s' (%d). Skipping.", lvar, l_info$dim[1], target, K))
            } else {
              linked_info_list[[lvar]] <- l_info
            }
          }
        }
      }

      iter_total   <- dim(f_arr)[1]
      chains_total <- dim(f_arr)[2]

      get_values <- function(arr, i, c, idx) { arr[i, c, idx] }
      set_values <- function(arr, i, c, idx, value) {
        arr[i, c, idx] <- value
        arr
      }

      for (c in seq_len(chains_total)) {
        for (i in seq_len(iter_total)) {

          if (t_info$loc == "fixed") t_vals <- get_values(f_arr, i, c, t_info$idx)
          else if (t_info$loc == "random") t_vals <- get_values(r_arr, i, c, t_info$idx)
          else if (t_info$loc == "tran") t_vals <- get_values(t_arr, i, c, t_info$idx)
          else if (t_info$loc == "gq") t_vals <- get_values(g_arr, i, c, t_info$idx)
          else stop("Target must be a vector variable.")

          ord <- order(t_vals)

          if (any(ord != seq_len(K))) {
            if (t_info$loc == "fixed") f_arr <- set_values(f_arr, i, c, t_info$idx, t_vals[ord])
            else if (t_info$loc == "random") r_arr <- set_values(r_arr, i, c, t_info$idx, t_vals[ord])
            else if (t_info$loc == "tran") t_arr <- set_values(t_arr, i, c, t_info$idx, t_vals[ord])
            else if (t_info$loc == "gq") g_arr <- set_values(g_arr, i, c, t_info$idx, t_vals[ord])

            for (lvar in names(linked_info_list)) {
              l_info <- linked_info_list[[lvar]]

              if (l_info$loc == "fixed") f_arr <- set_values(f_arr, i, c, l_info$idx, f_arr[i, c, l_info$idx][ord])
              else if (l_info$loc == "random") r_arr <- set_values(r_arr, i, c, l_info$idx, r_arr[i, c, l_info$idx][ord])
              else if (l_info$loc == "tran") t_arr <- set_values(t_arr, i, c, l_info$idx, t_arr[i, c, l_info$idx][ord])
              else if (l_info$loc == "gq") g_arr <- set_values(g_arr, i, c, l_info$idx, g_arr[i, c, l_info$idx][ord])
              else if (K == 2) {
                s_fn <- function(x) 1 - x
                if (!is.null(scalar_fns) && lvar %in% names(scalar_fns)) {
                  s_fn <- scalar_fns[[lvar]]
                }

                if (l_info$loc == "fixed_scalar") f_arr[i, c, l_info$idx] <- s_fn(f_arr[i, c, l_info$idx])
                else if (l_info$loc == "random_scalar") r_arr[i, c, l_info$idx] <- s_fn(r_arr[i, c, l_info$idx])
                else if (l_info$loc == "tran_scalar") t_arr[i, c, l_info$idx] <- s_fn(t_arr[i, c, l_info$idx])
                else if (l_info$loc == "gq_scalar") g_arr[i, c, l_info$idx] <- s_fn(g_arr[i, c, l_info$idx])
              }
            }
          }
        }
      }

      obj$fit <- f_arr
      obj$random_fit <- r_arr
      obj$transform_fit <- t_arr
      obj$generate_fit <- g_arr

      fixed_mean_new <- apply(obj$fit[, , -1, drop = FALSE], 3, mean)
      new_posterior_mean <- obj$posterior_mean
      new_posterior_mean[names(fixed_mean_new)] <- fixed_mean_new

      if (!is.null(obj$random_fit)) {
        random_mean_new <- apply(obj$random_fit, 3, mean)
        new_posterior_mean[names(random_mean_new)] <- random_mean_new
      }

      obj$posterior_mean <- new_posterior_mean

      if (isTRUE(overwrite)) invisible(self) else obj
    }

  )
)
