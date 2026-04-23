#' Automatic Differentiation Variational Inference (ADVI)
#'
#' @param model An RTMB objective function object (`ad_obj`).
#' @param par_list A list defining the structure of parameters to be estimated.
#' @param pl_full A list defining the full structure of parameters including random effects.
#' @param iter Integer; fixed number of iterations for the optimization. Default is 3000.
#' @param tol_rel_obj Numeric; relative tolerance for the ELBO to check convergence. Default is 0.001.
#' @param window_size Integer; size of the moving window to calculate the median ELBO. Default is 100.
#' @param num_samples Integer; number of posterior draws to generate after optimization. Default is 1000.
#' @param alpha Numeric; learning rate (step size) for the Adam optimizer. Default is 0.01.
#' @param laplace Logical; whether Laplace approximation is used. Default is FALSE.
#' @param print_freq Integer; frequency of printing progress to the console. Set to 0 to disable. Default is 500.
#' @param method Character; type of variational distribution. One of "meanfield", "fullrank", or "hybrid". Default is "meanfield".
#' @param update_progress Optional function to update a progress bar.
#' @param update_interval Integer; interval for updating the progress bar. Default is 100.
#' @return A list containing `fit`, `random_fit`, `elbo_history`, `elbo_final`, `rel_obj_final`, and `converged`.
ADVI_method <- function(model, par_list, pl_full,
                        iter = 3000, tol_rel_obj = 0.001,
                        window_size = 100, num_samples = 1000,
                        alpha = 0.01, laplace = FALSE,
                        print_freq = 500,
                        method = c("meanfield", "fullrank", "hybrid"),
                        update_progress = NULL, update_interval = 100) {

  method <- match.arg(method)

  P <- length(model$par)
  mu <- model$par
  par_names_rtmb <- names(model$par)

  # fixed and random effect index
  random_flags <- sapply(par_list, function(x) isTRUE(x$random))
  random_bases <- names(par_list)[random_flags]

  idx_random <- which(par_names_rtmb %in% random_bases)
  idx_fixed  <- which(!(par_names_rtmb %in% random_bases))

  P_fixed <- length(idx_fixed)
  P_random <- length(idx_random)

  if (method == "hybrid") {
    if (P_random == 0) {
      method <- "fullrank"
    } else if (P_fixed == 0) {
      method <- "meanfield"
    }
  }

  beta1 <- 0.9
  beta2 <- 0.999
  epsilon <- 1e-8
  entropy_const <- (P / 2) * log(2 * pi * exp(1))

  if (print_freq > 0) cat("Starting ADVI optimization with Adam...\n")

  if (method == "hybrid") {
    L_diag <- rep(-2, P_fixed)
    L_off <- rep(0, P_fixed * (P_fixed - 1) / 2)
    omega <- rep(-2, P_random)

    m_mu <- rep(0, P); v_mu <- rep(0, P)
    m_diag <- rep(0, P_fixed); v_diag <- rep(0, P_fixed)
    m_off <- rep(0, length(L_off)); v_off <- rep(0, length(L_off))
    m_omega <- rep(0, P_random); v_omega <- rep(0, P_random)
  } else if (method == "fullrank") {
    L_diag <- rep(-2, P)
    L_off <- rep(0, P * (P - 1) / 2)
    m_mu <- rep(0, P); v_mu <- rep(0, P)
    m_diag <- rep(0, P); v_diag <- rep(0, P)
    m_off <- rep(0, length(L_off)); v_off <- rep(0, length(L_off))
  } else { # meanfield
    omega <- rep(-2, P)
    m_mu <- rep(0, P); v_mu <- rep(0, P)
    m_omega <- rep(0, P); v_omega <- rep(0, P)
  }

  elbo_history <- numeric(iter)
  converged <- FALSE
  rel_obj_final <- NA

  mu_history <- matrix(NA, nrow = window_size, ncol = P)

  for (t in 1:iter) {
    # --- update progress bar ---
    if (!is.null(update_progress) && t %% update_interval == 0) {
      update_progress(1)
    }

    eps <- rnorm(P)

    if (method == "hybrid") {
      eps_fixed <- eps[idx_fixed]
      eps_random <- eps[idx_random]

      L <- matrix(0, P_fixed, P_fixed)
      if (length(L_off) > 0) L[lower.tri(L)] <- L_off
      diag(L) <- exp(L_diag)
      theta_fixed <- mu[idx_fixed] + as.vector(L %*% eps_fixed)

      sigma_random <- exp(omega)
      theta_random <- mu[idx_random] + sigma_random * eps_random

      theta <- numeric(P)
      theta[idx_fixed] <- theta_fixed
      theta[idx_random] <- theta_random
    } else if (method == "fullrank") {
      L <- matrix(0, P, P)
      if (length(L_off) > 0) L[lower.tri(L)] <- L_off
      diag(L) <- exp(L_diag)
      theta <- mu + as.vector(L %*% eps)
    } else {
      sigma <- exp(omega)
      theta <- mu + sigma * eps
    }

    fn_val <- -model$fn(theta)
    gr_val <- as.vector(-model$gr(theta))

    if (!is.finite(fn_val) || any(!is.finite(gr_val))) {
      elbo_history[t] <- if(t > 1) elbo_history[t-1] else 0
      next
    }

    if (method == "hybrid") {
      elbo_history[t] <- fn_val + sum(L_diag) + sum(omega) + entropy_const

      grad_mu <- gr_val

      # fixed effect (Full-rank)
      gr_fixed <- gr_val[idx_fixed]
      grad_L_mat <- outer(gr_fixed, eps_fixed)
      grad_diag <- diag(grad_L_mat) * exp(L_diag) + 1
      grad_off <- grad_L_mat[lower.tri(grad_L_mat)]

      # random effect (Mean-field)
      gr_random <- gr_val[idx_random]
      grad_omega <- gr_random * (eps_random * sigma_random) + 1

      m_mu <- beta1 * m_mu + (1 - beta1) * grad_mu
      v_mu <- beta2 * v_mu + (1 - beta2) * (grad_mu^2)
      mu <- mu + alpha * (m_mu / (1 - beta1^t)) / (sqrt(v_mu / (1 - beta2^t)) + epsilon)

      m_diag <- beta1 * m_diag + (1 - beta1) * grad_diag
      v_diag <- beta2 * v_diag + (1 - beta2) * (grad_diag^2)
      L_diag <- L_diag + alpha * (m_diag / (1 - beta1^t)) / (sqrt(v_diag / (1 - beta2^t)) + epsilon)

      if (length(L_off) > 0) {
        m_off <- beta1 * m_off + (1 - beta1) * grad_off
        v_off <- beta2 * v_off + (1 - beta2) * (grad_off^2)
        L_off <- L_off + alpha * (m_off / (1 - beta1^t)) / (sqrt(v_off / (1 - beta2^t)) + epsilon)
      }

      if (length(omega) > 0) {
        m_omega <- beta1 * m_omega + (1 - beta1) * grad_omega
        v_omega <- beta2 * v_omega + (1 - beta2) * (grad_omega^2)
        omega <- omega + alpha * (m_omega / (1 - beta1^t)) / (sqrt(v_omega / (1 - beta2^t)) + epsilon)
      }
    } else if (method == "fullrank") {
      elbo_history[t] <- fn_val + sum(L_diag) + entropy_const
      grad_mu <- gr_val
      grad_L_mat <- outer(gr_val, eps)
      grad_diag <- diag(grad_L_mat) * exp(L_diag) + 1
      grad_off <- grad_L_mat[lower.tri(grad_L_mat)]

      m_mu <- beta1 * m_mu + (1 - beta1) * grad_mu
      v_mu <- beta2 * v_mu + (1 - beta2) * (grad_mu^2)
      mu <- mu + alpha * (m_mu / (1 - beta1^t)) / (sqrt(v_mu / (1 - beta2^t)) + epsilon)

      m_diag <- beta1 * m_diag + (1 - beta1) * grad_diag
      v_diag <- beta2 * v_diag + (1 - beta2) * (grad_diag^2)
      L_diag <- L_diag + alpha * (m_diag / (1 - beta1^t)) / (sqrt(v_diag / (1 - beta2^t)) + epsilon)

      if (length(L_off) > 0) {
        m_off <- beta1 * m_off + (1 - beta1) * grad_off
        v_off <- beta2 * v_off + (1 - beta2) * (grad_off^2)
        L_off <- L_off + alpha * (m_off / (1 - beta1^t)) / (sqrt(v_off / (1 - beta2^t)) + epsilon)
      }
    } else {
      elbo_history[t] <- fn_val + sum(omega) + entropy_const
      grad_mu <- gr_val
      grad_omega <- gr_val * (eps * sigma) + 1

      m_mu <- beta1 * m_mu + (1 - beta1) * grad_mu
      v_mu <- beta2 * v_mu + (1 - beta2) * (grad_mu^2)
      mu <- mu + alpha * (m_mu / (1 - beta1^t)) / (sqrt(v_mu / (1 - beta2^t)) + epsilon)

      m_omega <- beta1 * m_omega + (1 - beta1) * grad_omega
      v_omega <- beta2 * v_omega + (1 - beta2) * (grad_omega^2)
      omega <- omega + alpha * (m_omega / (1 - beta1^t)) / (sqrt(v_omega / (1 - beta2^t)) + epsilon)
    }

    idx <- (t - 1) %% window_size + 1
    mu_history[idx, ] <- mu

    if (print_freq > 0 && t %% print_freq == 0) {
      cat(sprintf("Iter %d: Approx ELBO = %.2f\n", t, elbo_history[t]))
    }
  }

  # convergence index
  check_start <- 2 * window_size
  if (iter > check_start) {
    med_prev <- median(elbo_history[(iter - 2 * window_size + 1):(iter - window_size)])
    med_curr <- median(elbo_history[(iter - window_size + 1):iter])
    rel_obj_final <- abs(med_curr - med_prev) / (abs(med_prev) + 1e-8)
    if (rel_obj_final < tol_rel_obj) {
      converged <- TRUE
    }
  }

  if (print_freq > 0) cat("Generating posterior samples from variational distribution...\n")

  fit_matrix <- matrix(NA, nrow = num_samples, ncol = P)

  if (method == "hybrid") {
    L <- matrix(0, P_fixed, P_fixed)
    if (length(L_off) > 0) L[lower.tri(L)] <- L_off
    diag(L) <- exp(L_diag)
    sigma_random <- exp(omega)

    for (i in 1:num_samples) {
      eps_samp <- rnorm(P)
      eps_fixed <- eps_samp[idx_fixed]
      eps_random <- eps_samp[idx_random]

      theta_fixed <- mu[idx_fixed] + as.vector(L %*% eps_fixed)
      theta_random <- mu[idx_random] + sigma_random * eps_random

      theta_samp <- numeric(P)
      theta_samp[idx_fixed] <- theta_fixed
      theta_samp[idx_random] <- theta_random
      fit_matrix[i, ] <- theta_samp
    }
  } else if (method == "fullrank") {
    L <- matrix(0, P, P)
    if (length(L_off) > 0) L[lower.tri(L)] <- L_off
    diag(L) <- exp(L_diag)
    for (i in 1:num_samples) {
      fit_matrix[i, ] <- mu + as.vector(L %*% rnorm(P))
    }
  } else {
    sigma <- exp(omega)
    for (i in 1:num_samples) {
      fit_matrix[i, ] <- mu + sigma * rnorm(P)
    }
  }

  base_names_full <- sub("\\[.*\\]", "", pl_full$names)
  random_flags <- sapply(par_list, function(x) isTRUE(x$random))

  if (laplace && any(random_flags)) {
    base_names_fixed <- names(par_list)[!random_flags]
    base_names_random <- names(par_list)[random_flags]
  } else {
    base_names_fixed <- names(par_list)
    base_names_random <- character(0)
  }

  fixed_names <- pl_full$names[base_names_full %in% base_names_fixed]
  random_names <- pl_full$names[base_names_full %in% base_names_random]

  fixed_idx  <- which(pl_full$names %in% fixed_names)
  random_idx <- if (length(random_names) > 0) which(pl_full$names %in% random_names) else integer(0)

  P_all_true <- length(pl_full$names)
  para_final <- array(NA, dim = c(num_samples, P_all_true))
  lp_final <- numeric(num_samples)

  for (i in 1:num_samples) {
    zeta_sample <- fit_matrix[i, ]
    lp_final[i] <- -model$fn(zeta_sample)

    if (laplace && length(model$env$random) > 0) {
      para_list_res <- model$env$parList(x = model$env$last.par)
    } else {
      para_list_res <- model$env$parList(x = zeta_sample)
    }

    con_list <- to_constrained(para_list_res, par_list)
    para_final[i, ] <- unlist(con_list, use.names = FALSE)
  }

  if (method == "hybrid") {
    elbo_final <- mean(lp_final) + sum(L_diag) + sum(omega) + entropy_const
  } else if (method == "fullrank") {
    elbo_final <- mean(lp_final) + sum(L_diag) + entropy_const
  } else {
    elbo_final <- mean(lp_final) + sum(omega) + entropy_const
  }

  fit <- array(NA, dim = c(num_samples, 1, length(fixed_names) + 1))
  dimnames(fit) <- list(iteration = NULL, chain = "est1", variable = c("lp", fixed_names))

  random_fit <- NULL
  if (length(random_names) > 0) {
    random_fit <- array(NA, dim = c(num_samples, 1, length(random_names)))
    dimnames(random_fit) <- list(iteration = NULL, chain = "est1", variable = random_names)
  }

  for (i in 1:num_samples) {
    fit[i, 1, 1] <- lp_final[i]
    fit[i, 1, -1] <- para_final[i, fixed_idx]
    if (!is.null(random_fit)) {
      random_fit[i, 1, ] <- para_final[i, random_idx]
    }
  }

  valid_hist_idx <- which(mu_hist_filled & apply(mu_history, 1, function(z) all(is.finite(z))))

  mu_hist_constrained <- matrix(
    NA,
    nrow = length(valid_hist_idx),
    ncol = length(fixed_idx) + length(random_idx)
  )
  colnames(mu_hist_constrained) <- c(fixed_names, random_names)

  for (j in seq_along(valid_hist_idx)) {
    i <- valid_hist_idx[j]
    zeta_hist <- mu_history[i, ]

    model$fn(zeta_hist)
    if (laplace && length(model$env$random) > 0) {
      para_list_res <- model$env$parList(x = model$env$last.par)
    } else {
      para_list_res <- model$env$parList(x = zeta_hist)
    }

    con_list <- to_constrained(para_list_res, par_list)
    para_hist_all <- unlist(con_list, use.names = FALSE)

    mu_hist_constrained[j, 1:length(fixed_idx)] <- para_hist_all[fixed_idx]
    if (length(random_idx) > 0) {
      mu_hist_constrained[j, (length(fixed_idx) + 1):ncol(mu_hist_constrained)] <- para_hist_all[random_idx]
    }
  }

  fit <- Re(fit)
  if (!is.null(random_fit)) random_fit <- Re(random_fit)

  return(list(
    fit           = fit,
    random_fit    = random_fit,
    elbo_history  = elbo_history,
    elbo_final    = elbo_final,
    rel_obj_final = rel_obj_final,
    converged     = converged,
    mu_history    = mu_hist_constrained
  ))
}
