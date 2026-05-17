.calculate_satterthwaite_df_impl <- function(self, private, ad_obj, idx_fix_active = NULL, L_u_total = NULL, opt_par = NULL, max_df = NULL, silent = FALSE, return_sensitivities = FALSE) {
  # if (!silent) cat("Estimating Satterthwaite degrees of freedom...\n")

  par <- if (!is.null(opt_par)) opt_par else ad_obj$par
  P <- length(par)
  total_len <- if (!is.null(L_u_total)) L_u_total else P
  df_full <- rep(Inf, total_len)

  # Step 1: Compute Hessian at optimum via Jacobian of analytical gradient
  # For marginal likelihood (Laplace), numerical differentiation on the gradient is most reliable

  H0 <- tryCatch(private$.simple_jacobian(ad_obj$gr, par), warning = function(w) NULL, error = function(e) NULL)

  if (is.null(H0) || any(!is.finite(H0))) {
    if (!silent) cat("SE warning: Hessian computation failed in Satterthwaite calculation. Using normal approximation.\n")
    return(df_full)
  }
  V <- tryCatch(solve(H0), error = function(e) {
    if (!silent) cat("SE warning: Hessian matrix was singular in Satterthwaite calculation; using MASS::ginv().\n")
    tryCatch(MASS::ginv(H0), error = function(e2) NULL)
  })
  if (is.null(V)) {
    if (!silent) cat("SE warning: Satterthwaite covariance calculation failed. Using normal approximation.\n")
    return(df_full)
  }

  # Step 2: Compute dV_ii/dtheta_k via central differences of the Hessian
  # Optimized: Compute dH/dtheta_k once per k and vectorize diag(V dH V)
  # Improved: Smaller eps base for better sensitivity to 3rd derivatives
  # Improved: Slightly larger eps base to smooth out 3rd derivative noise
  eps <- 1e-3 * pmax(abs(par), 0.1)
  grad_V_diag <- matrix(0, nrow = P, ncol = P) # grad_V_diag[k, i] = dV_ii / dtheta_k
  dH_list <- if (return_sensitivities) vector("list", P) else NULL

  for (k in 1:P) {
    par_plus <- par_minus <- par
    par_plus[k] <- par[k] + eps[k]
    par_minus[k] <- par[k] - eps[k]

    H_plus <- tryCatch(private$.simple_jacobian(ad_obj$gr, par_plus), warning = function(w) NULL, error = function(e) NULL)
    H_minus <- tryCatch(private$.simple_jacobian(ad_obj$gr, par_minus), warning = function(w) NULL, error = function(e) NULL)

    if (is.null(H_plus) || is.null(H_minus)) next

    # Third derivative approximation: dH/dtheta_k
    dH_k <- (H_plus - H_minus) / (2 * eps[k])
    if (return_sensitivities) dH_list[[k]] <- dH_k

    # dV/dtheta_k = -V %*% dH_k %*% V
    # Diagonal elements: diag(dV/dtheta_k) = -rowSums((V %*% dH_k) * V)
    grad_V_diag[k, ] <- -rowSums((V %*% dH_k) * V)
  }

  rel_tol <- 1e-10
  df_par <- rep(Inf, P)
  for (i in 1:P) {
    grad_vi <- grad_V_diag[, i]
    if (all(abs(grad_vi) < 1e-30)) next

    # Var(V_ii) = grad(V_ii)^T %*% V %*% grad(V_ii)
    var_vi <- as.numeric(t(grad_vi) %*% V %*% grad_vi)
    
    # Numerical stability: If variance of V_ii is extremely small relative to V_ii^2, 
    # it means sensitivity is negligible, so DF should be Infinite.
    if (!is.na(var_vi) && is.finite(var_vi) && var_vi > rel_tol * (V[i, i]^2) && V[i, i] > 0) {
      df_par[i] <- 2 * (V[i, i]^2) / var_vi
    } else {
      df_par[i] <- Inf
    }
  }

  if (!is.null(max_df)) {
    df_par[df_par > max_df] <- Inf
  }
  df_par[!is.finite(df_par)] <- Inf
  df_par <- pmax(df_par, 2.1)

  # Map to full parameter vector
  if (!is.null(idx_fix_active) && length(idx_fix_active) == P) {
    df_full[idx_fix_active] <- df_par
  } else if (P == total_len) {
    df_full <- df_par
  }

  finite_dfs <- df_par[is.finite(df_par)]

  if (return_sensitivities) {
    return(list(df = df_full, sensitivities = dH_list, V = V))
  }
  return(df_full)
}

.calculate_reml_satterthwaite_df_impl <- function(self, private, ad_obj, opt_par, beta_idx, max_df = NULL, silent = FALSE, return_sensitivities = FALSE) {
  if (!silent) cat("Estimating Satterthwaite degrees of freedom for fixed effects (REML)...\n")

  P <- length(opt_par)
  n_beta <- length(beta_idx)

  # Helper to get covariance of joint Hessian for beta
  get_beta_cov = function(theta) {
    p_full <- ad_obj$env$last.par.best
    idx_fixed <- ad_obj$env$lfixed()
    p_full[idx_fixed] <- theta

    H <- tryCatch(ad_obj$env$spHess(p_full, random = TRUE), warning = function(w) NULL, error = function(e) NULL)
    if (is.null(H)) {
      return(NULL)
    }

    # Convert to dense matrix for solve to ensure extraction works reliably
    V_full <- tryCatch(solve(as.matrix(H)), warning = function(w) NULL, error = function(e) NULL)
    if (is.null(V_full)) {
      return(NULL)
    }

    as.matrix(V_full[beta_idx, beta_idx, drop = FALSE])
  }

  # Early return helper to ensure consistent return types
  fallback_reml_result <- function() {
    if (return_sensitivities) {
      return(list(
        df = rep(Inf, n_beta),
        se = rep(NA_real_, n_beta),
        sensitivities = vector("list", P), # Keep list structure
        V_theta = NULL,
        V_beta = NULL
      ))
    }
    list(
      df = rep(Inf, n_beta),
      se = rep(NA_real_, n_beta)
    )
  }

  V_beta_0 <- tryCatch(get_beta_cov(opt_par), error = function(e) {
    return(NULL)
  })
  if (is.null(V_beta_0)) {
    return(fallback_reml_result())
  }
  V_beta_diag_0 <- diag(V_beta_0)

  # Covariance of variance components (theta)
  H_theta <- tryCatch(private$.simple_jacobian(ad_obj$gr, opt_par), warning = function(w) NULL, error = function(e) NULL)
  if (is.null(H_theta)) {
    return(fallback_reml_result())
  }

  V_theta <- tryCatch(solve(H_theta), error = function(e) {
    if (!silent) cat("SE warning: Hessian matrix was singular in REML Satterthwaite calculation; using MASS::ginv().\n")
    tryCatch(MASS::ginv(H_theta), error = function(e2) NULL)
  })
  if (is.null(V_theta)) {
    if (!silent) cat("SE warning: REML Satterthwaite covariance calculation failed. Using infinite degrees of freedom.\n")
    return(fallback_reml_result())
  }

  # Gradient of V_beta with respect to theta
  eps <- 1e-5 * pmax(abs(opt_par), 0.1)
  grad_V_beta_diag <- matrix(0, nrow = P, ncol = n_beta)
  dV_beta_list <- if (return_sensitivities) vector("list", P) else NULL

  for (k in 1:P) {
    p_plus <- p_minus <- opt_par
    p_plus[k] <- opt_par[k] + eps[k]
    p_minus[k] <- opt_par[k] - eps[k]

    v_plus <- tryCatch(get_beta_cov(p_plus), warning = function(w) NULL, error = function(e) NULL)
    v_minus <- tryCatch(get_beta_cov(p_minus), warning = function(w) NULL, error = function(e) NULL)

    if (is.null(v_plus) || is.null(v_minus)) next
    dV_k <- (v_plus - v_minus) / (2 * eps[k])
    if (return_sensitivities) dV_beta_list[[k]] <- dV_k
    grad_V_beta_diag[k, ] <- diag(dV_k)
  }

  dfs <- rep(Inf, n_beta)
  rel_tol <- 1e-10
  for (i in 1:n_beta) {
    grad_vi <- grad_V_beta_diag[, i]
    var_vi <- as.numeric(t(grad_vi) %*% V_theta %*% grad_vi)
    v_ii <- V_beta_diag_0[i]
    
    if (!is.na(var_vi) && is.finite(var_vi) && var_vi > rel_tol * (v_ii^2) && v_ii > 0) {
      dfs[i] <- 2 * (v_ii^2) / var_vi
    } else {
      dfs[i] <- Inf
    }
  }

  if (!is.null(max_df)) {
    dfs[dfs > max_df] <- Inf
  }
  dfs[!is.finite(dfs)] <- Inf
  dfs <- pmax(dfs, 2.1)
  if (return_sensitivities) {
    return(list(df = dfs, se = sqrt(V_beta_diag_0), sensitivities = dV_beta_list, V_theta = V_theta, V_beta = V_beta_0))
  }
  return(list(df = dfs, se = sqrt(V_beta_diag_0)))
}
