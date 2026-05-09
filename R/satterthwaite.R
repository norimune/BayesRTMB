#' Calculate Satterthwaite degrees of freedom for fixed effects
#'
#' @param model An \code{RTMB_Model} object.
#' @param ad_obj An RTMB AD object.
#' @param idx_fix_active Integer vector of active fixed parameter indices.
#' @param L_u_total Total number of unconstrained parameters.
#' @param opt_par Vector of optimized parameter values.
#' @return A vector of estimated degrees of freedom.
#' @keywords internal
.calculate_satterthwaite_df_impl <- function(model, ad_obj, idx_fix_active = NULL, L_u_total = NULL, opt_par = NULL) {
  P <- if (!is.null(L_u_total)) L_u_total else length(ad_obj$par)
  df_full <- rep(Inf, P)
  
  if (is.null(opt_par)) opt_par <- ad_obj$par

  # Robust numerical Hessian calculation
  H0 <- tryCatch(simple_jacobian(ad_obj$gr, opt_par), warning = function(w) NULL, error = function(e) NULL)
  if (is.null(H0)) return(df_full)

  h <- 1e-4
  p_active <- length(opt_par)
  df_par <- rep(Inf, p_active)

  for (i in 1:p_active) {
    h_i <- max(abs(opt_par[i]), 1e-4) * h
    par_plus <- par_minus <- opt_par
    par_plus[i] <- opt_par[i] + h_i
    par_minus[i] <- opt_par[i] - h_i

    H_plus <- tryCatch(simple_jacobian(ad_obj$gr, par_plus), warning = function(w) NULL, error = function(e) NULL)
    H_minus <- tryCatch(simple_jacobian(ad_obj$gr, par_minus), warning = function(w) NULL, error = function(e) NULL)

    if (!is.null(H_plus) && !is.null(H_minus)) {
      dH_dtheta <- (H_plus - H_minus) / (2 * h_i)
      var_theta <- tryCatch(solve(H0), error = function(e) NULL)
      if (!is.null(var_theta)) {
        g <- dH_dtheta[i, i]
        V_g <- 2 * sum(diag(var_theta %*% dH_dtheta %*% var_theta %*% dH_dtheta))
        df_val <- 2 * (g^2) / V_g
        df_par[i] <- max(df_val, 1e-6)
      }
    }
  }

  if (!is.null(idx_fix_active) && length(idx_fix_active) == p_active) {
    df_full[idx_fix_active] <- df_par
  } else {
    df_full <- df_par
  }

  return(df_full)
}

#' Calculate Satterthwaite degrees of freedom for REML models
#'
#' @param model An \code{RTMB_Model} object.
#' @param ad_obj An RTMB AD object.
#' @param opt_par Vector of optimized parameter values.
#' @param target_ran_idx Indices of target parameters within the random effect vector.
#' @return A vector of estimated degrees of freedom for the target parameters.
#' @keywords internal
.calculate_reml_satterthwaite_df_impl <- function(model, ad_obj, opt_par, target_ran_idx) {
  n_target <- length(target_ran_idx)
  df_out <- rep(Inf, n_target)

  # 0. 初期チェック: 関数と勾配が現在のパラメータで評価可能か
  f0 <- tryCatch(ad_obj$fn(opt_par), error = function(e) NULL)
  g0 <- tryCatch(ad_obj$gr(opt_par), error = function(e) NULL)
  
  if (is.null(f0) || is.null(g0) || any(is.na(g0))) {
    return(df_out)
  }

  # (A) 分散パラメータの分散共分散行列 (V_theta)
  # 数値微分によるヘッセ行列の推定
  H_theta <- tryCatch(simple_jacobian(ad_obj$gr, opt_par), warning = function(w) NULL, error = function(e) NULL)
  
  if (is.null(H_theta) || any(is.na(H_theta))) {
    # 失敗した場合は勾配の外部積で近似 (BHHH)
    H_theta <- outer(g0, g0)
  }
  
  # 正則化
  diag(H_theta) <- diag(H_theta) + 1e-8
  V_theta <- tryCatch(solve(H_theta), error = function(e) NULL)
  
  if (is.null(V_theta)) {
    V_theta <- tryCatch(MASS::ginv(H_theta), error = function(e) NULL)
  }
  
  if (is.null(V_theta)) return(df_out)

  h <- 1e-4
  P <- length(opt_par)

  # Function to get variance of fixed effects given theta
  get_fixed_var <- function(theta) {
    # Ensure internal state is updated to theta (finds mode of random effects)
    # Use tryCatch to prevent hard crash
    lp <- tryCatch(ad_obj$fn(theta), error = function(e) NA)
    if (is.na(lp)) return(rep(NA, n_target))
    
    # --- Plan A: Use sdreport (most robust, forces re-calc) ---
    sdr <- tryCatch(RTMB::sdreport(ad_obj, getJointPrecision = TRUE), error = function(e) NULL)
    
    if (!is.null(sdr) && !is.null(sdr$jointPrecision)) {
      H_joint <- sdr$jointPrecision
      diag(H_joint) <- diag(H_joint) + 1e-9
      V_joint <- tryCatch(solve(H_joint), error = function(e) NULL)
      
      if (!is.null(V_joint)) {
        # random indices in joint precision matrix
        idx_ran_in_joint <- ad_obj$env$random
        target_joint_idx <- idx_ran_in_joint[target_ran_idx]
        return(diag(V_joint)[target_joint_idx])
      }
    }
    
    # --- Fallback Plan D: spHess with manual internal state reset ---
    # Force evaluation at order 0 to refresh internal tape pointers
    tryCatch(ad_obj$env$f(ad_obj$env$last.par, order = 0), error = function(e) NULL)
    
    H_full <- tryCatch(ad_obj$env$spHess(random = TRUE), error = function(e) NULL)
    if (is.null(H_full)) return(rep(NA, n_target))
    
    diag(H_full) <- diag(H_full) + 1e-9
    V_full <- tryCatch(solve(H_full), error = function(e) NULL)
    if (is.null(V_full)) {
      eig_f <- tryCatch(eigen(as.matrix(H_full), symmetric = TRUE), error = function(e) NULL)
      if (is.null(eig_f)) return(rep(NA, n_target))
      eig_f$values <- pmax(eig_f$values, 1e-10)
      V_full <- eig_f$vectors %*% diag(1/eig_f$values) %*% t(eig_f$vectors)
    }
    return(diag(V_full)[target_ran_idx])
  }

  V_base <- get_fixed_var(opt_par)
  if (any(is.na(V_base))) {
    message("Satterthwaite: V_base calculation failed.")
    return(df_out)
  }

  message(sprintf("Satterthwaite: Calculating sensitivity for %d parameters...", P))
  # Matrix to store gradients: [n_target x P]
  grad_V <- matrix(0, nrow = n_target, ncol = P)
  
  for (i in 1:P) {
    h_i <- max(abs(opt_par[i]), 1e-4) * h
    par_plus <- par_minus <- opt_par
    par_plus[i] <- opt_par[i] + h_i
    par_minus[i] <- opt_par[i] - h_i

    V_plus <- get_fixed_var(par_plus)
    V_minus <- get_fixed_var(par_minus)
    
    if (all(!is.na(V_plus)) && all(!is.na(V_minus))) {
      grad_V[, i] <- (V_plus - V_minus) / (2 * h_i)
    } else {
      # Fallback to one-sided difference if one fails
      if (all(!is.na(V_plus))) {
        grad_V[, i] <- (V_plus - V_base) / h_i
      } else if (all(!is.na(V_minus))) {
        grad_V[, i] <- (V_base - V_minus) / h_i
      }
    }
  }

  for (j in 1:n_target) {
    grad_Vj <- grad_V[j, ]
    denom <- as.numeric(t(grad_Vj) %*% V_theta %*% grad_Vj)
    if (!is.na(denom) && denom > 1e-16) {
      df_val <- 2 * (V_base[j]^2) / denom
      df_out[j] <- max(df_val, 1e-6)
    }
  }

  return(df_out)
}
