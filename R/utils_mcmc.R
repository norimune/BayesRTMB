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

  # --- Fix: Add check.names = FALSE ---
  test_dat <- read.csv(test_file, header = TRUE, check.names = FALSE)
  n_samples <- nrow(test_dat)

  orig_pl <- model$par_list
  random_flags <- sapply(orig_pl, function(x) isTRUE(x$random))

  if (laplace && any(random_flags)) {
    pl_fixed  <- BayesRTMB:::parse_parameters(orig_pl[!random_flags], model$par_names)
    pl_random <- BayesRTMB:::parse_parameters(orig_pl[random_flags], model$par_names)
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

    # --- Fix: Add check.names = FALSE ---
    dat <- read.csv(file_path, header = TRUE, check.names = FALSE)

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
# Helper to assign class to summary results
as_summary_df <- function(df) {
  class(df) <- c("summary_BayesRTMB", "data.frame")
  return(df)
}

#' print for summary_BayesRTMB class
#' @method print summary_BayesRTMB
#' @param x An object of class summary_BayesRTMB
#' @param digits integer
#' @param ... Additional arguments.
#' @export
print.summary_BayesRTMB <- function(x, digits = NULL,...) {
  df <- x

  if (is.null(digits)) {
    digits <- attr(x, "digits")
    if (is.null(digits)) digits <- 2
  }

  # 1. Convert each column to string (specify display format)
  out_char <- as.data.frame(lapply(names(df), function(cn) {
    val <- df[[cn]]
    if (is.numeric(val)) {
      if (grepl("ess|iter|count", cn, ignore.case = TRUE)) {
        # ESS and similar metrics are integers (including NA handling)
        ifelse(is.na(val), "NA", sprintf("%.0f", val))
      } else {
        fmt <- paste0("%.", digits, "f")
        ifelse(is.na(val), "NA", sprintf(fmt, val))
      }
    } else {
      ifelse(is.na(val), "NA", as.character(val))
    }
  }), stringsAsFactors = FALSE)
  colnames(out_char) <- names(df)

  col_widths <- sapply(seq_along(names(out_char)), function(i) {
    max(nchar(names(out_char)[i]), nchar(out_char[, i]), na.rm = TRUE)
  })

  # 2. Calculate column widths
  col_widths <- sapply(seq_along(names(out_char)), function(i) {
    max(nchar(names(out_char)[i]), nchar(out_char[, i]))
  })

  # 3. Display header (all titles right-aligned)
  header_parts <- sapply(seq_along(col_widths), function(i) {
    sprintf(paste0("%", col_widths[i], "s"), names(out_char)[i])
  })
  cat(paste(header_parts, collapse = "  "), "\n")

  # 4. Display data (only the first column is left-aligned)
  for (r in seq_len(nrow(out_char))) {
    row_parts <- sapply(seq_along(col_widths), function(c) {
      if (c == 1) {
        # Left-align only the variable name data
        sprintf(paste0("%-", col_widths[c], "s"), out_char[r, c])
      } else {
        # Right-align numeric data
        sprintf(paste0("%", col_widths[c], "s"), out_char[r, c])
      }
    })
    cat(paste(row_parts, collapse = "  "), "\n")
  }

  return(invisible(x))
}

#' Calculate Conditional Effects
#' @param fit Model fit object.
#' @param effect Name of the variable to visualize the effect.
#' @param ... Additional arguments.
#' @export
conditional_effects <- function(fit, effect, ...) {
  UseMethod("conditional_effects")
}

#' Calculate conditional effects for MCMC fit objects
#' @method conditional_effects mcmc_fit
#' @param fit An object of class `MCMC_Fit`.
#' @param effect Name of the explanatory variable to visualize (e.g., "X1" or "X1:X2").
#' @param resolution Grid resolution to calculate for continuous variables (default is 100).
#' @param prob Probability for the credible interval (default is 0.95).
#' @param ... Additional arguments.
#' @export
conditional_effects.mcmc_fit <- function(fit, effect, resolution = 100, prob = 0.95, ...) {
  model_obj <- fit$model
  if (is.null(model_obj$formula) || is.null(model_obj$raw_data)) {
    stop("This model object does not contain a formula or the original data.")
  }

  form <- model_obj$formula
  raw_data <- model_obj$raw_data
  fam <- model_obj$family

  # Check for interaction and split
  eff_vars <- strsplit(effect, ":")[[1]]
  if (length(eff_vars) > 2) {
    stop("Interaction plots for 3 or more variables are not currently supported.")
  }
  eff1 <- eff_vars[1]
  eff2 <- if (length(eff_vars) == 2) eff_vars[2] else NULL

  # 1. Create base data frame (fixing other variables to mean/mode)
  base_data <- lapply(raw_data, function(x) {
    if (is.numeric(x)) {
      mean(x, na.rm = TRUE)
    } else if (is.factor(x) || is.character(x)) {
      tbl <- table(x)
      factor(names(tbl)[which.max(tbl)], levels = levels(as.factor(x)))
    } else {
      x[1]
    }
  })
  base_data <- as.data.frame(base_data)

  # 2. Create newdata moving only the effect variable
  val1 <- raw_data[[eff1]]
  if (is.null(val1)) stop(sprintf("Variable '%s' not found in the data.", eff1))

  is_numeric1 <- is.numeric(val1)
  if (is_numeric1) {
    seq1 <- seq(min(val1, na.rm = TRUE), max(val1, na.rm = TRUE), length.out = resolution)
  } else {
    seq1 <- sort(unique(val1))
  }

  if (!is.null(eff2)) {
    val2 <- raw_data[[eff2]]
    if (is.null(val2)) stop(sprintf("Variable '%s' not found in the data.", eff2))

    # If the second variable is continuous and has many unique values, restrict to 3 representative points (Mean-1SD, Mean, Mean+1SD)
    if (is.numeric(val2) && length(unique(val2)) > 5) {
      seq2 <- c(mean(val2, na.rm=TRUE) - sd(val2, na.rm=TRUE),
                mean(val2, na.rm=TRUE),
                mean(val2, na.rm=TRUE) + sd(val2, na.rm=TRUE))
      seq2 <- round(seq2, 2) # Round for better display
    } else {
      seq2 <- sort(unique(val2))
    }

    grid_data <- expand.grid(eff1 = seq1, eff2 = seq2)
    names(grid_data) <- c(eff1, eff2)

    newdata <- base_data[rep(1, nrow(grid_data)), , drop = FALSE]
    newdata[[eff1]] <- grid_data[[eff1]]
    newdata[[eff2]] <- grid_data[[eff2]]
  } else {
    newdata <- base_data[rep(1, length(seq1)), , drop = FALSE]
    newdata[[eff1]] <- seq1
  }

  # 3. Create design matrix
  rhs <- delete.response(terms(form))
  X_new <- model.matrix(rhs, data = newdata)

  # 4. Obtain posterior samples (fixed effects only)
  beta_samples <- fit$draws(pars = "beta", inc_random = FALSE, inc_transform = FALSE, inc_generate = FALSE)
  I <- dim(beta_samples)[1]
  C <- dim(beta_samples)[2]
  P <- dim(beta_samples)[3]
  beta_flat <- matrix(beta_samples, nrow = I * C, ncol = P)

  # 5. Calculate linear predictor
  eta <- X_new %*% t(beta_flat)

  # 6. Convert to expected value via inverse link function
  if (is.null(fam)) fam <- "gaussian"
  inv_link <- switch(fam,
                     "gaussian" = , "lognormal" = , "student_t" = function(x) x,
                     "poisson" = , "neg_binomial" = , "gamma" = exp,
                     "bernoulli" = , "binomial" = plogis,
                     function(x) x
  )
  mu <- inv_link(eta)

  # 7. Summarize posterior distribution
  alpha <- 1 - prob
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2

  res_df <- data.frame(
    estimate   = apply(mu, 1, mean),
    lower      = apply(mu, 1, quantile, probs = lower_q),
    upper      = apply(mu, 1, quantile, probs = upper_q)
  )
  res_df <- cbind(newdata[, eff_vars, drop = FALSE], res_df)

  # 8. Combine results into a list and assign class
  res <- list(data = res_df, effect_vars = eff_vars, is_numeric = is_numeric1)
  class(res) <- "ce_rtmb"
  return(res)
}

#' Plot method for ce_rtmb class (Base R)
#' @method plot ce_rtmb
#' @param x An object of class ce_rtmb
#' @param ... Additional arguments.
#' @export
plot.ce_rtmb <- function(x, ...) {
  df <- x$data
  eff_vars <- x$effect_vars
  eff1 <- eff_vars[1]
  has_interaction <- length(eff_vars) > 1

  x_val <- df[[eff1]]
  y_est <- df$estimate
  y_low <- df$lower
  y_up  <- df$upper

  if (!has_interaction) {
    # --- Without interaction ---
    col_line <- rgb(0, 0.45, 0.7)
    col_ribbon <- rgb(0, 0.45, 0.7, 0.2)

    if (x$is_numeric) {
      plot(x_val, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlab = eff1, ylab = "Predicted value", main = paste("Conditional effect of", eff1), ...)
      polygon(c(x_val, rev(x_val)), c(y_low, rev(y_up)), col = col_ribbon, border = NA)
      lines(x_val, y_est, col = col_line, lwd = 2)
    } else {
      x_num <- as.numeric(as.factor(x_val))
      plot(x_num, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlim = c(min(x_num) - 0.5, max(x_num) + 0.5), xaxt = "n",
           xlab = eff1, ylab = "Predicted value", main = paste("Conditional effect of", eff1), ...)
      axis(1, at = x_num, labels = as.character(x_val))
      segments(x0 = x_num, y0 = y_low, x1 = x_num, y1 = y_up, col = col_line, lwd = 2)
      points(x_num, y_est, col = col_line, pch = 16, cex = 1.5)
    }

  } else {
    # --- With interaction ---
    eff2 <- eff_vars[2]
    groups <- unique(df[[eff2]])
    n_groups <- length(groups)

    # Create color palette based on number of groups
    cols_line <- hcl.colors(n_groups, palette = "Dark 2")
    cols_ribbon <- sapply(cols_line, function(col) {
      rgb_val <- col2rgb(col) / 255
      rgb(rgb_val[1], rgb_val[2], rgb_val[3], 0.2) # Transparency 0.2
    })

    if (x$is_numeric) {
      plot(x_val, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlab = eff1, ylab = "Predicted value",
           main = paste("Conditional effect of", eff1, "by", eff2), ...)

      for (i in seq_along(groups)) {
        idx <- df[[eff2]] == groups[i]
        xv <- df[[eff1]][idx]
        polygon(c(xv, rev(xv)), c(df$lower[idx], rev(df$upper[idx])), col = cols_ribbon[i], border = NA)
        lines(xv, df$estimate[idx], col = cols_line[i], lwd = 2)
      }
    } else {
      x_fct <- as.factor(x_val)
      x_num <- as.numeric(x_fct)

      plot(x_num, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlim = c(min(x_num) - 0.5, max(x_num) + 0.5), xaxt = "n",
           xlab = eff1, ylab = "Predicted value",
           main = paste("Conditional effect of", eff1, "by", eff2), ...)
      axis(1, at = unique(x_num), labels = levels(x_fct))

      # Shift X-axis slightly to avoid overlapping error bars
      offset_step <- 0.1
      offsets <- seq(-offset_step * (n_groups-1)/2, offset_step * (n_groups-1)/2, length.out = n_groups)

      for (i in seq_along(groups)) {
        idx <- df[[eff2]] == groups[i]
        xv <- x_num[idx] + offsets[i]
        segments(x0 = xv, y0 = df$lower[idx], x1 = xv, y1 = df$upper[idx], col = cols_line[i], lwd = 2)
        points(xv, df$estimate[idx], col = cols_line[i], pch = 16, cex = 1.5)
      }
    }

    # Add legend
    legend("topright", title = eff2, legend = format(groups, digits = 3),
           col = cols_line, lty = 1, pch = ifelse(x$is_numeric, NA, 16), lwd = 2, bty = "n")
  }
  invisible(x)
}

#' Print method for ce_rtmb class (automatically calls plot)
#' @method print ce_rtmb
#' @param x An object of class ce_rtmb
#' @param ... Additional arguments.
#' @export
print.ce_rtmb <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}

#' Calculate Item Information Function
#' @param x An object of class RTMB_Fit_Base
#' @param ... Additional arguments.
#' @export
item_info <- function(x, ...) UseMethod("item_info")

#' Calculate Test Information Function
#' @param x An object of class RTMB_Fit_Base
#' @param ... Additional arguments.
#' @export
test_info <- function(x, ...) UseMethod("test_info")

#' Item Information Function for RTMB_Fit_Base
#' @method item_info RTMB_Fit_Base
#' @param x An object of class RTMB_Fit_Base.
#' @param theta_seq Sequence of trait values (ability) to evaluate.
#' @param items Index or item names to restrict the calculation to specific items (optional).
#' @param ... Additional arguments.
#' @export
item_info.RTMB_Fit_Base <- function(x, theta_seq = seq(-4, 4, length.out = 100), items = NULL, ...) {
  # MAP_Fit uses $par, MCMC/VB use EAP()
  est <- if (!is.null(x$par)) x$par else x$EAP()
  b <- est$b

  # Fix: Refer to the correct item names (par_names) saved during model construction instead of names(b)
  par_names_b <- x$model$par_names$b

  # Get data type and item names
  if (is.matrix(b)) {
    type <- "ordered"
    J <- nrow(b)
    K_cat <- ncol(b) + 1
    item_names <- if (is.list(par_names_b)) par_names_b[[1]] else paste0("Item", 1:J)
  } else {
    type <- "binary"
    J <- length(b)
    item_names <- if (!is.null(par_names_b)) par_names_b else paste0("Item", 1:J)
  }

  # Identify indices of target items
  if (!is.null(items)) {
    if (is.character(items)) {
      target_idx <- match(items, item_names)
      if (any(is.na(target_idx))) {
        warning("Specified item names not found: ", paste(items[is.na(target_idx)], collapse = ", "))
        target_idx <- target_idx[!is.na(target_idx)]
      }
    } else {
      target_idx <- items[items >= 1 & items <= J]
    }
    if (length(target_idx) == 0) stop("No valid items selected.")
  } else {
    target_idx <- 1:J
  }

  a <- if (!is.null(est$a)) est$a else rep(1.0, J)
  c_param <- if (!is.null(est$c)) est$c else rep(0.0, J)

  info_mat <- matrix(0, nrow = length(theta_seq), ncol = length(target_idx))
  colnames(info_mat) <- item_names[target_idx]

  for (idx in seq_along(target_idx)) {
    j <- target_idx[idx]
    if (type == "binary") {
      eta <- a[j] * (theta_seq - b[j])
      P <- c_param[j] + (1 - c_param[j]) * plogis(eta)
      Q <- 1 - P

      info_mat[, idx] <- (a[j]^2 * Q / P) * ((P - c_param[j]) / (1 - c_param[j]))^2

    } else if (type == "ordered") {
      P_star <- matrix(0, nrow = length(theta_seq), ncol = K_cat + 1)
      P_star[, 1] <- 1.0
      P_star[, K_cat + 1] <- 0.0

      for (k in 1:(K_cat - 1)) {
        P_star[, k + 1] <- plogis(b[j, k] - a[j] * theta_seq)
      }

      item_info_j <- rep(0, length(theta_seq))
      for (k in 1:K_cat) {
        Pk <- P_star[, k] - P_star[, k + 1]

        dP_star_k   <- -a[j] * P_star[, k] * (1 - P_star[, k])
        dP_star_kp1 <- -a[j] * P_star[, k + 1] * (1 - P_star[, k + 1])

        if (k == 1) dP_star_k <- rep(0, length(theta_seq))
        if (k == K_cat) dP_star_kp1 <- rep(0, length(theta_seq))

        dPk <- dP_star_k - dP_star_kp1

        Pk_safe <- pmax(Pk, 1e-10)
        item_info_j <- item_info_j + (dPk^2) / Pk_safe
      }
      info_mat[, idx] <- item_info_j
    }
  }

  res <- list(theta = theta_seq, info = info_mat)
  class(res) <- "rtmb_item_info"
  return(res)
}

#' @method item_info mcmc_fit
#' @export
item_info.mcmc_fit <- item_info.RTMB_Fit_Base
#' @method item_info advi_fit
#' @export
item_info.advi_fit <- item_info.RTMB_Fit_Base
#' @method item_info map_fit
#' @export
item_info.map_fit <- item_info.RTMB_Fit_Base

#' @method test_info RTMB_Fit_Base
#' @export
test_info.RTMB_Fit_Base <- function(x, theta_seq = seq(-4, 4, length.out = 100), ...) {
  ii <- item_info(x, theta_seq, ...)
  test_info_vec <- rowSums(ii$info)

  res <- list(theta = theta_seq, info = test_info_vec)
  class(res) <- "rtmb_test_info"
  return(res)
}

#' @method test_info mcmc_fit
#' @export
test_info.mcmc_fit <- test_info.RTMB_Fit_Base
#' @method test_info advi_fit
#' @export
test_info.advi_fit <- test_info.RTMB_Fit_Base
#' @method test_info map_fit
#' @export
test_info.map_fit <- test_info.RTMB_Fit_Base

# --- Plotting Methods ---

#' @method plot rtmb_item_info
#' @export
plot.rtmb_item_info <- function(x, legend = TRUE, ...) {
  matplot(x$theta, x$info, type = "l", lty = 1,
          xlab = expression(theta ~ "(Ability)"), ylab = "Information",
          main = "Item Information Functions", ...)

  # Add legend (if legend = TRUE)
  if (legend && ncol(x$info) > 0) {
    legend("topright", legend = colnames(x$info), col = 1:ncol(x$info), lty = 1, cex = 0.8)
  }
}

#' @method plot rtmb_test_info
#' @export
plot.rtmb_test_info <- function(x, ...) {
  plot(x$theta, x$info, type = "l", lwd = 2, col = "blue",
       xlab = expression(theta ~ "(Ability)"), ylab = "Information",
       main = "Test Information Function", ...)
}

#' @method print rtmb_item_info
#' @export
print.rtmb_item_info <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}

#' @method print rtmb_test_info
#' @export
print.rtmb_test_info <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}


#' Calculate Item Response Curve / Category Response Curve
#' @export
#' @param x An object of class
#' @param ... Additional arguments.
item_curve <- function(x, ...) UseMethod("item_curve")

#' Item Response Curve for RTMB_Fit_Base
#' @method item_curve RTMB_Fit_Base
#' @param x An object of class RTMB_Fit_Base.
#' @param theta_seq Sequence of trait values (ability) to evaluate.
#' @param items Index or item names to restrict the calculation to specific items (optional).
#' @param ... Additional arguments.
#' @export
item_curve.RTMB_Fit_Base <- function(x, theta_seq = seq(-4, 4, length.out = 100), items = NULL, ...) {
  est <- if (!is.null(x$par)) x$par else x$EAP()
  b <- est$b

  par_names_b <- x$model$par_names$b

  # Get data type and item names
  if (is.matrix(b)) {
    type <- "ordered"
    J <- nrow(b)
    K_cat <- ncol(b) + 1
    item_names <- if (is.list(par_names_b)) par_names_b[[1]] else paste0("Item", 1:J)
  } else {
    type <- "binary"
    J <- length(b)
    item_names <- if (!is.null(par_names_b)) par_names_b else paste0("Item", 1:J)
  }

  # Identify indices of target items
  if (!is.null(items)) {
    if (is.character(items)) {
      target_idx <- match(items, item_names)
      if (any(is.na(target_idx))) {
        warning("Specified item names not found: ", paste(items[is.na(target_idx)], collapse = ", "))
        target_idx <- target_idx[!is.na(target_idx)]
      }
    } else {
      target_idx <- items[items >= 1 & items <= J]
    }
    if (length(target_idx) == 0) stop("No valid items selected.")
  } else {
    target_idx <- 1:J
  }

  a <- if (!is.null(est$a)) est$a else rep(1.0, J)
  c_param <- if (!is.null(est$c)) est$c else rep(0.0, J)

  out_curves <- list()

  for (idx in seq_along(target_idx)) {
    j <- target_idx[idx]

    if (type == "binary") {
      # Probability of correct response for binary model
      eta <- a[j] * (theta_seq - b[j])
      P <- c_param[j] + (1 - c_param[j]) * plogis(eta)

      out_curves[[item_names[j]]] <- matrix(P, ncol = 1)
      colnames(out_curves[[item_names[j]]]) <- "P(Y=1)"

    } else if (type == "ordered") {
      # Category selection probability for ordered model (GRM)
      eta <- a[j] * theta_seq

      # Cumulative probability P(Y <= k)
      P_cum <- matrix(1.0, nrow = length(theta_seq), ncol = K_cat + 1)
      P_cum[, 1] <- 0.0 # P(Y <= 0) = 0

      for (k in 1:(K_cat - 1)) {
        P_cum[, k + 1] <- plogis(b[j, k] - eta)
      }
      # P_cum[, K_cat + 1] is initialized to 1.0 (P(Y <= K) = 1)

      # Probability of each category P(Y = k) = P(Y <= k) - P(Y <= k-1)
      P_cat <- matrix(0, nrow = length(theta_seq), ncol = K_cat)
      for (k in 1:K_cat) {
        P_cat[, k] <- P_cum[, k + 1] - P_cum[, k]
      }

      colnames(P_cat) <- paste0("Cat", 1:K_cat)
      out_curves[[item_names[j]]] <- P_cat
    }
  }

  res <- list(theta = theta_seq, curves = out_curves, type = type)
  class(res) <- "rtmb_item_curve"
  return(res)
}

#' @method item_curve mcmc_fit
#' @export
item_curve.mcmc_fit <- item_curve.RTMB_Fit_Base
#' @method item_curve advi_fit
#' @export
item_curve.advi_fit <- item_curve.RTMB_Fit_Base
#' @method item_curve map_fit
#' @export
item_curve.map_fit <- item_curve.RTMB_Fit_Base

#' @method plot rtmb_item_curve
#' @export
plot.rtmb_item_curve <- function(x, legend = TRUE, ...) {
  n_items <- length(x$curves)

  if (x$type == "binary") {
    # Binary model: Draw all items overlaid on a single plot
    P_mat <- do.call(cbind, x$curves)
    matplot(x$theta, P_mat, type = "l", lty = 1, ylim = c(0, 1),
            xlab = expression(theta ~ "(Ability)"), ylab = "Probability",
            main = "Item Response Curves", ...)

    if (legend && ncol(P_mat) > 0) {
      legend("topleft", legend = names(x$curves), col = 1:n_items, lty = 1, cex = 0.8)
    }

  } else if (x$type == "ordered") {
    # Ordered model: Draw category response curves in separate panels per item
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    if (n_items > 1) {
      nr <- ceiling(sqrt(n_items))
      nc <- ceiling(n_items / nr)
      par(mfrow = c(nr, nc))
    }

    for (item_name in names(x$curves)) {
      P_mat <- x$curves[[item_name]]
      matplot(x$theta, P_mat, type = "l", lty = 1:ncol(P_mat), ylim = c(0, 1),
              xlab = expression(theta), ylab = "Probability",
              main = paste("Category Curves:", item_name), ...)

      if (legend) {
        legend("topleft", legend = colnames(P_mat), col = 1:ncol(P_mat),
               lty = 1:ncol(P_mat), cex = 0.7)
      }
    }
  }
}

#' @method print rtmb_item_curve
#' @export
print.rtmb_item_curve <- function(x, ...) {
  plot(x, ...)
  invisible(x)
}
#' Sort and display factor loadings neatly
#'
#' @param loadings Matrix, data frame, or list of factor loadings (if list, the first element is used)
#' @param cutoff Absolute loadings below this value will be displayed as blank (default is 0.0)
#' @param round_digits Number of decimal places to display (default is 3)
#' @return Sorted loading matrix (returned invisibly, allowing assignment to a variable)
#' @export
sort_loadings <- function(loadings, cutoff = 0.0, round_digits = 3) {

  if (is.list(loadings) && !is.data.frame(loadings)) {
    if (length(loadings) == 0) stop("The input list is empty.")
    loadings <- loadings[[1]]
  }

  load_mat <- as.matrix(loadings)

  if (is.null(rownames(load_mat))) {
    rownames(load_mat) <- paste0("V", 1:nrow(load_mat))
  }

  max_factor <- apply(abs(load_mat), 1, which.max)
  max_loading <- apply(abs(load_mat), 1, max)

  df_order <- data.frame(
    var_name = rownames(load_mat),
    factor = max_factor,
    loading = max_loading,
    stringsAsFactors = FALSE
  )

  df_order <- df_order[order(df_order$factor, -df_order$loading), ]
  sorted_mat <- load_mat[df_order$var_name, , drop = FALSE]

  print_mat <- matrix("", nrow = nrow(sorted_mat), ncol = ncol(sorted_mat))
  rownames(print_mat) <- rownames(sorted_mat)
  colnames(print_mat) <- colnames(sorted_mat)

  # Create a string to identify zeros like "-.000"
  zero_str <- paste0("-.", paste(rep("0", round_digits), collapse = ""))

  for (i in 1:nrow(sorted_mat)) {
    for (j in 1:ncol(sorted_mat)) {
      val <- sorted_mat[i, j]
      if (abs(val) >= cutoff) {
        s <- sprintf(paste0("%.", round_digits, "f"), val)
        s <- sub("^0\\.", ".", s)
        s <- sub("^-0\\.", "-.", s)

        # Replace "-.000" with ".000"
        if (s == zero_str) {
          s <- sub("^-", "", s)
        }

        print_mat[i, j] <- s
      }
    }
  }

  print_df <- as.data.frame(print_mat, stringsAsFactors = FALSE)
  print(print_df, quote = FALSE, right = TRUE)

  return(invisible(sorted_mat))
}
#' Calculate Bayes factor from log marginal likelihoods
#'
#' @param logml1 Log marginal likelihood of Model 1 (e.g., target model)
#' @param logml2 Log marginal likelihood of Model 2 (e.g., reference/null model)
#' @return An object containing Bayes factor, log Bayes factor, estimation error, and interpretation
#' @export
bayes_factor <- function(logml1, logml2) {
  # Strip attributes (error, ess) and calculate as pure numbers
  val1 <- as.numeric(logml1)
  val2 <- as.numeric(logml2)

  # Calculate log Bayes factor and Bayes factor
  log_bf <- val1 - val2
  bf <- exp(log_bf)

  # Error propagation (calculate if logml1 and logml2 contain "error" attribute)
  err1 <- attr(logml1, "error")
  err2 <- attr(logml2, "error")

  has_error <- !is.null(err1) && !is.null(err2) && !is.na(err1) && !is.na(err2)

  if (has_error) {
    # Standard error assuming the two MCMC samplings are independent
    log_bf_err <- sqrt(err1^2 + err2^2)
  } else {
    log_bf_err <- NA_real_
  }

  # Interpretation of evidence strength based on Jeffreys (1961) / Kass & Raftery (1995)
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
}

#' Print method for bayes_factor objects
#' @param x An object of class bayes_factor.
#' @param digits Number of decimal places to display (default is 4).
#' @param ... Additional arguments.
#' @export
print.bayes_factor <- function(x, digits = 4, ...) {
  cat("Bayes Factor (BF12) :", round(x$BF12, digits), "\n")

  if (!is.na(x$log_BF_error)) {
    cat(sprintf("Log Bayes Factor    : %.4f (Approx. Error = %.4f)\n", x$log_BF12, x$log_BF_error))
  } else {
    cat(sprintf("Log Bayes Factor    : %.4f\n", x$log_BF12))
  }

  cat("Interpretation      :", x$interpretation, "\n")
  invisible(x)
}
