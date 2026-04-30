#' Calculate Conditional Effects
#'
#' @description
#' Calculate and visualize the predicted values (marginal effects) of a variable,
#' potentially conditional on the levels of another variable (interaction).
#'
#' @param fit Model fit object (e.g., `mcmc_fit`).
#' @param effect Name of the variable to visualize (e.g., "X1" or "X1:X2").
#' @param sd_multiplier Numeric. Multiplier for standard deviation when splitting continuous moderators (default is 1).
#' @param ... Additional arguments.
#'
#' @return An object of class `ce_rtmb` containing the predicted values and their credible intervals.
#'
#' @examples
#' \dontrun{
#'   fit <- rtmb_lm(mpg ~ wt * hp, data = mtcars)
#'   mcmc_fit <- fit$sample()
#'   ce <- conditional_effects(mcmc_fit, effect = "wt:hp")
#'   plot(ce)
#'   summary(ce)
#' }
#' @export
conditional_effects <- function(fit, effect, sd_multiplier = 1, ...) {
  UseMethod("conditional_effects")
}

#' Calculate conditional effects for MCMC fit objects
#' @method conditional_effects mcmc_fit
#' @param fit An object of class `MCMC_Fit`.
#' @param effect Name of the explanatory variable to visualize (e.g., "X1" or "X1:X2").
#' @param resolution Grid resolution to calculate for continuous variables (default is 100).
#' @param prob Probability for the credible interval (default is 0.95).
#' @param sd_multiplier Multiplier for SD for continuous moderators.
#' @param ... Additional arguments.
#'
#' @return A `ce_rtmb` object.
#' @export
conditional_effects.mcmc_fit <- function(fit, effect, resolution = 100, prob = 0.95, sd_multiplier = 1, ...) {
  model_obj <- fit$model
  if (is.null(model_obj$formula) || is.null(model_obj$raw_data)) {
    stop("This model object does not contain a formula or the original data.")
  }

  form <- model_obj$formula
  raw_data <- model_obj$raw_data
  fam <- model_obj$family
  X_mean <- model_obj$data$X_mean # For back-transforming Intercept_c if needed

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
      seq2 <- c(mean(val2, na.rm=TRUE) - sd_multiplier * sd(val2, na.rm=TRUE),
                mean(val2, na.rm=TRUE),
                mean(val2, na.rm=TRUE) + sd_multiplier * sd(val2, na.rm=TRUE))
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
  # Ensure all factor levels are known to model.matrix
  mf_raw <- model.frame(rhs, data = raw_data, na.action = na.pass)
  xlev <- .getXlevels(terms(rhs), mf_raw)
  X_new <- model.matrix(rhs, data = newdata, xlev = xlev)

  # 4. Obtain posterior samples (fixed effects)
  # Handle both wrapper functions (b, Intercept_c) and custom models (beta)
  all_pars <- dimnames(fit$draws(inc_transform=TRUE, inc_generate=FALSE))[[3]]
  
  if (any(grepl("^beta(\\[|$)", all_pars))) {
    beta_samples <- fit$draws(pars = "beta", inc_random = FALSE, inc_transform = FALSE, inc_generate = FALSE)
  } else if (any(grepl("^b(\\[|$)", all_pars))) {
    b_samps <- fit$draws(pars = "b", inc_random = FALSE, inc_transform = FALSE, inc_generate = FALSE)
    
    # Try to get Intercept (transformed) or Intercept_c
    if ("Intercept" %in% all_pars) {
      int_samps <- fit$draws(pars = "Intercept", inc_transform = TRUE)
    } else if ("Intercept_c" %in% all_pars) {
      # Back-calculate Intercept from Intercept_c if X_mean is available
      ic_samps <- fit$draws(pars = "Intercept_c")
      if (!is.null(X_mean)) {
        b_mat <- matrix(b_samps, nrow = dim(b_samps)[1] * dim(b_samps)[2], ncol = dim(b_samps)[3])
        int_val <- as.numeric(ic_samps) - (b_mat %*% as.numeric(X_mean))
        int_samps <- array(int_val, dim = c(dim(b_samps)[1], dim(b_samps)[2], 1))
      } else {
        int_samps <- ic_samps
      }
    } else {
      int_samps <- NULL
    }
    
    # Combine Intercept and b to match model.matrix output (Intercept first)
    if (!is.null(int_samps)) {
      I <- dim(b_samps)[1]; C <- dim(b_samps)[2]; P <- dim(b_samps)[3]
      beta_samples <- array(NA, dim = c(I, C, P + 1))
      beta_samples[,,1] <- int_samps
      beta_samples[,,2:(P+1)] <- b_samps
    } else {
      beta_samples <- b_samps
    }
  } else {
    stop("Could not find fixed effect parameters (b or beta) in the fit object.")
  }

  I <- dim(beta_samples)[1]
  C <- dim(beta_samples)[2]
  P <- dim(beta_samples)[3]
  beta_flat <- matrix(beta_samples, nrow = I * C, ncol = P)

  # Check dimension compatibility
  if (ncol(X_new) != ncol(beta_flat)) {
    # If mismatch, try to adjust X_new (e.g., if Intercept is missing in beta but present in X_new)
    if (ncol(X_new) == ncol(beta_flat) + 1 && colnames(X_new)[1] == "(Intercept)") {
      X_new <- X_new[, -1, drop = FALSE]
    } else if (ncol(beta_flat) == ncol(X_new) + 1) {
       # Maybe beta has intercept but X_new doesn't (rare but possible)
       X_new <- cbind(1, X_new)
    } else {
      stop(sprintf("Dimension mismatch: Design matrix has %d columns, but fixed effects have %d parameters.", ncol(X_new), ncol(beta_flat)))
    }
  }

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
      x_fct <- as.factor(x_val)
      x_num <- as.numeric(x_fct)
      lvls <- levels(x_fct)
      
      plot(x_num, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlim = c(0.5, length(lvls) + 0.5), xaxt = "n",
           xlab = eff1, ylab = "Predicted value", main = paste("Conditional effect of", eff1), ...)
      axis(1, at = 1:length(lvls), labels = lvls)
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
      lvls <- levels(x_fct)

      plot(x_num, y_est, type = "n", ylim = range(c(y_low, y_up)),
           xlim = c(0.5, length(lvls) + 0.5), xaxt = "n",
           xlab = eff1, ylab = "Predicted value",
           main = paste("Conditional effect of", eff1, "by", eff2), ...)
      axis(1, at = 1:length(lvls), labels = lvls)

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

#' Summary method for ce_rtmb class
#' @method summary ce_rtmb
#' @param object An object of class ce_rtmb
#' @param ... Additional arguments.
#' @export
summary.ce_rtmb <- function(object, ...) {
  return(object$data)
}


#' Calculate Simple Effects
#'
#' @description
#' Calculate the effect of a focal variable at different levels of a moderator.
#' For categorical focal variables, it calculates pairwise differences (contrasts).
#' For continuous focal variables, it calculates the slope (simple slopes).
#'
#' @param fit Model fit object (mcmc_fit).
#' @param effect Character string of the interaction (e.g., "A:B"). The first variable is the focal variable.
#' @param sd_multiplier Multiplier for SD for continuous moderators (default is 1).
#' @param ... Additional arguments.
#'
#' @return A `ce_simple` object (data frame) containing the estimated effects and their credible intervals.
#'
#' @examples
#' \dontrun{
#'   fit <- rtmb_lm(len ~ supp * dose, data = ToothGrowth)
#'   mcmc_fit <- fit$sample()
#'   # Effect of supplement at each dose level
#'   se <- simple_effects(mcmc_fit, effect = "supp:dose")
#'   print(se)
#' }
#' @export
simple_effects <- function(fit, effect, sd_multiplier = 1, ...) {
  UseMethod("simple_effects")
}

#' Simple effects for MCMC fit objects
#' @method simple_effects mcmc_fit
#' @param fit An object of class `MCMC_Fit`.
#' @param effect Interaction term (e.g., "A:B").
#' @param prob Probability for credible intervals.
#' @param sd_multiplier Multiplier for SD for continuous moderators.
#' @param ... Additional arguments.
#' @export
simple_effects.mcmc_fit <- function(fit, effect, prob = 0.95, sd_multiplier = 1, ...) {
  eff_vars <- strsplit(effect, ":")[[1]]
  if (length(eff_vars) != 2) {
    stop("simple_effects requires an interaction of exactly two variables (e.g., 'A:B').")
  }
  
  focal <- eff_vars[1]
  mod <- eff_vars[2]
  
  # 1. Reuse conditional_effects logic to get predicted values
  # We use a higher resolution for moderators if continuous
  ce <- conditional_effects(fit, effect = effect, resolution = 10, sd_multiplier = sd_multiplier, ...)
  df <- ce$data
  
  # Identify levels/values of the moderator
  mod_vals <- unique(df[[mod]])
  
  res_list <- list()
  
  if (ce$is_numeric) {
    # Simple Slopes for continuous focal variable
    eps <- 1e-4
    results <- data.frame()
    
    # [Internal Helper] Get predicted samples (shared with categorical logic)
    get_pred_samples <- function(fit, newdata) {
      model_obj <- fit$model
      form <- model_obj$formula
      raw_data <- model_obj$raw_data
      X_mean <- model_obj$data$X_mean
      rhs <- delete.response(terms(form))
      
      # Ensure all factor levels are known to model.matrix
      mf_raw <- model.frame(rhs, data = raw_data, na.action = na.pass)
      xlev <- .getXlevels(terms(rhs), mf_raw)
      X_new <- model.matrix(rhs, data = newdata, xlev = xlev)
      
      all_pars <- dimnames(fit$draws(inc_transform=TRUE))[[3]]
      if (any(grepl("^beta(\\[|$)", all_pars))) {
        beta_samples <- fit$draws(pars = "beta")
      } else {
        b_samps <- fit$draws(pars = "b")
        if ("Intercept" %in% all_pars) {
          int_samps <- fit$draws(pars = "Intercept", inc_transform = TRUE)
        } else if ("Intercept_c" %in% all_pars) {
          ic_samps <- fit$draws(pars = "Intercept_c")
          b_mat <- matrix(b_samps, nrow = dim(b_samps)[1] * dim(b_samps)[2], ncol = dim(b_samps)[3])
          int_val <- as.numeric(ic_samps) - (b_mat %*% as.numeric(X_mean))
          int_samps <- array(int_val, dim = c(dim(b_samps)[1], dim(b_samps)[2], 1))
        } else {
          int_samps <- NULL
        }
        if (!is.null(int_samps)) {
          I <- dim(b_samps)[1]; C <- dim(b_samps)[2]; P <- dim(b_samps)[3]
          beta_samples <- array(NA, dim = c(I, C, P + 1))
          beta_samples[,,1] <- int_samps
          beta_samples[,,2:(P+1)] <- b_samps
        } else {
          beta_samples <- b_samps
        }
      }
      
      I <- dim(beta_samples)[1]; C <- dim(beta_samples)[2]; P <- dim(beta_samples)[3]
      beta_flat <- matrix(beta_samples, nrow = I * C, ncol = P)
      
      if (ncol(X_new) == ncol(beta_flat) + 1 && colnames(X_new)[1] == "(Intercept)") {
        X_new <- X_new[, -1, drop = FALSE]
      }
      
      eta <- X_new %*% t(beta_flat)
      inv_link <- switch(model_obj$family %||% "gaussian",
                         "gaussian" = function(x) x, "poisson" = exp, "bernoulli" = plogis, function(x) x)
      return(inv_link(eta))
    }

    for (mv in mod_vals) {
      # Use the mean of the focal variable as the point of evaluation for the slope
      f_val <- mean(df[[focal]], na.rm=TRUE)
      
      # Create newdata for x and x + eps
      # We take one row of original data as template to keep other variables at their baseline
      template_data <- fit$model$raw_data[1, , drop=FALSE]
      # Set moderator value
      template_data[[mod]] <- mv
      
      sub_newdata1 <- template_data
      sub_newdata1[[focal]] <- f_val
      sub_newdata2 <- template_data
      sub_newdata2[[focal]] <- f_val + eps
      
      combined_newdata <- rbind(sub_newdata1, sub_newdata2)
      pred_samples <- get_pred_samples(fit, combined_newdata)
      
      # Slope = (y2 - y1) / eps
      slope_samples <- (pred_samples[2, ] - pred_samples[1, ]) / eps
      
      res_row <- data.frame(
        moderator = mod,
        mod_value = mv,
        term = paste("Slope of", focal),
        estimate = mean(slope_samples),
        lower = quantile(slope_samples, (1-prob)/2),
        upper = quantile(slope_samples, 1-(1-prob)/2)
      )
      results <- rbind(results, res_row)
    }
    
    names(results)[2] <- mod
    class(results) <- c("ce_simple", "data.frame")
    return(results)
    
  } else {
    # Pairwise differences for categorical focal variable
    focal_lvls <- unique(df[[focal]])
    if (length(focal_lvls) < 2) stop("Focal variable must have at least 2 levels.")
    
    # We need the posterior samples of the predicted values to calculate the differences
    # (re-running the core logic of conditional_effects but keeping samples)
    
    # [Internal Helper] Get predicted samples
    get_pred_samples <- function(fit, newdata) {
      model_obj <- fit$model
      form <- model_obj$formula
      X_mean <- model_obj$data$X_mean
      rhs <- delete.response(terms(form))
      X_new <- model.matrix(rhs, data = newdata)
      
      all_pars <- dimnames(fit$draws(inc_transform=TRUE))[[3]]
      if (any(grepl("^beta(\\[|$)", all_pars))) {
        beta_samples <- fit$draws(pars = "beta")
      } else {
        b_samps <- fit$draws(pars = "b")
        if ("Intercept" %in% all_pars) {
          int_samps <- fit$draws(pars = "Intercept", inc_transform = TRUE)
        } else if ("Intercept_c" %in% all_pars) {
          ic_samps <- fit$draws(pars = "Intercept_c")
          b_mat <- matrix(b_samps, nrow = dim(b_samps)[1] * dim(b_samps)[2], ncol = dim(b_samps)[3])
          int_val <- as.numeric(ic_samps) - (b_mat %*% as.numeric(X_mean))
          int_samps <- array(int_val, dim = c(dim(b_samps)[1], dim(b_samps)[2], 1))
        } else {
          int_samps <- NULL
        }
        if (!is.null(int_samps)) {
          I <- dim(b_samps)[1]; C <- dim(b_samps)[2]; P <- dim(b_samps)[3]
          beta_samples <- array(NA, dim = c(I, C, P + 1))
          beta_samples[,,1] <- int_samps
          beta_samples[,,2:(P+1)] <- b_samps
        } else {
          beta_samples <- b_samps
        }
      }
      
      I <- dim(beta_samples)[1]; C <- dim(beta_samples)[2]; P <- dim(beta_samples)[3]
      beta_flat <- matrix(beta_samples, nrow = I * C, ncol = P)
      
      if (ncol(X_new) == ncol(beta_flat) + 1 && colnames(X_new)[1] == "(Intercept)") {
        X_new <- X_new[, -1, drop = FALSE]
      }
      
      eta <- X_new %*% t(beta_flat)
      inv_link <- switch(model_obj$family %||% "gaussian",
                         "gaussian" = function(x) x, "poisson" = exp, "bernoulli" = plogis, function(x) x)
      return(inv_link(eta))
    }
    
    # Calculate differences for each level of moderator
    results <- data.frame()
    
    for (mv in mod_vals) {
      sub_newdata <- df[df[[mod]] == mv, names(df) %in% names(fit$model$raw_data)]
      # Ensure focal levels are in consistent order
      sub_newdata <- sub_newdata[order(sub_newdata[[focal]]), ]
      
      pred_samples <- get_pred_samples(fit, sub_newdata) # (N_levels x N_samples)
      
      # All pairwise differences
      lvls <- sub_newdata[[focal]]
      for (i in 1:(length(lvls)-1)) {
        for (j in (i+1):length(lvls)) {
          diff_samples <- pred_samples[j, ] - pred_samples[i, ]
          
          res_row <- data.frame(
            moderator = mod,
            mod_value = mv,
            term = paste(lvls[j], "-", lvls[i]),
            estimate = mean(diff_samples),
            lower = quantile(diff_samples, (1-prob)/2),
            upper = quantile(diff_samples, 1-(1-prob)/2)
          )
          results <- rbind(results, res_row)
        }
      }
    }
    
    names(results)[2] <- mod
    class(results) <- c("ce_simple", "data.frame")
    return(results)
  }
}

#' Print method for ce_simple
#' @method print ce_simple
#' @export
print.ce_simple <- function(x, ...) {
  cat("--- Simple Effects Analysis ---\n")
  print.data.frame(x, row.names = FALSE)
  cat("\n")
  invisible(x)
}


#' Calculate Item Information Function
#' @param x An object of class RTMB_Fit_Base
#' @param ... Additional arguments.
#' @examples
#' \dontrun{
#'   fit <- rtmb_irt(data = BigFive[, 1:5], model = "2pl")
#'   map_fit <- fit$optimize()
#'   ii <- item_info(map_fit)
#'   plot(ii)
#' }
#' @export
item_info <- function(x, ...) UseMethod("item_info")

#' Calculate Test Information Function
#' @param x An object of class RTMB_Fit_Base
#' @param ... Additional arguments.
#' @examples
#' \dontrun{
#'   fit <- rtmb_irt(data = BigFive[, 1:5], model = "2pl")
#'   map_fit <- fit$optimize()
#'   ti <- test_info(map_fit)
#'   plot(ti)
#' }
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
#' @examples
#' \dontrun{
#'   fit <- rtmb_irt(data = BigFive[, 1:5], model = "2pl")
#'   map_fit <- fit$optimize()
#'   ic <- item_curve(map_fit)
#'   plot(ic)
#' }
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
#' @param logml1 Log marginal likelihood of Model 1, or a fitted model object (e.g., `mcmc_fit`, `map_fit`, `advi_fit`).
#' @param logml2 Log marginal likelihood of Model 2, or a fitted model object.
#' @param error_threshold Numeric; threshold for the approximate error warning. Default is 0.2.
#' @return An object containing Bayes factor, log Bayes factor, estimation error, and interpretation
#' @examples
#' \dontrun{
#'   # Compare two models using Bayes Factor
#'   fit1 <- rtmb_lm(mpg ~ wt, data = mtcars)
#'   fit2 <- rtmb_lm(mpg ~ wt + hp, data = mtcars)
#'   mcmc1 <- fit1$sample(sampling = 500, warmup = 500)
#'   mcmc2 <- fit2$sample(sampling = 500, warmup = 500)
#'   bf <- bayes_factor(mcmc1, mcmc2)
#'   print(bf)
#' }
#' @export
bayes_factor <- function(logml1, logml2, error_threshold = 0.2) {
  # Helper function to extract or calculate log marginal likelihood
  get_lml <- function(x, name = "Model") {
    if (inherits(x, "mcmc_fit")) {
      if (is.null(x$log_ml)) {
        cat(sprintf("Calculating marginal likelihood for %s...\n", name))
        x$log_ml <- x$bridgesampling()
      }
      return(x$log_ml)
    } else if (inherits(x, "map_fit")) {
      return(x$log_ml)
    } else if (inherits(x, "advi_fit")) {
      # Use the best ELBO as an approximation of log marginal likelihood
      return(max(x$ELBO, na.rm = TRUE))
    } else if (is.numeric(x)) {
      return(x)
    } else {
      stop(sprintf("Unsupported object type for %s. Must be numeric or a fit object.", name))
    }
  }

  val1_obj <- get_lml(logml1, "Model 1")
  val2_obj <- get_lml(logml2, "Model 2")

  # Strip attributes (error, ess) and calculate as pure numbers
  val1 <- as.numeric(val1_obj)
  val2 <- as.numeric(val2_obj)

  # Calculate log Bayes factor and Bayes factor
  log_bf <- val1 - val2
  bf <- exp(log_bf)

  # Error propagation (calculate if logml1 and logml2 contain "error" attribute)
  err1 <- attr(val1_obj, "error")
  err2 <- attr(val2_obj, "error")

  has_error <- !is.null(err1) && !is.null(err2) && !is.na(err1) && !is.na(err2)

  if (has_error) {
    # Standard error assuming the two MCMC samplings are independent
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

  # Interpretation of evidence strength based on Jeffreys (1961) / Kass & Raftery (1995)
  if (is.na(bf) || is.nan(bf)) {
    evidence <- "Indeterminate"
  } else if (bf > 100) {
    evidence <- "Decisive evidence for Model 1"
  } else if (bf > 10) {
    evidence <- "Strong evidence for Model 1"
  } else if (bf > 3) {
    evidence <- "Substantial evidence for Model 1"
  } else if (bf > 1) {
    evidence <- "Anecdotal evidence for Model 1"
  } else if (bf == 1) {
    evidence <- "No evidence"
  } else if (bf >= 1/3) {
    evidence <- "Anecdotal evidence for Model 2"
  } else if (bf >= 1/10) {
    evidence <- "Substantial evidence for Model 2"
  } else if (bf >= 1/100) {
    evidence <- "Strong evidence for Model 2"
  } else {
    evidence <- "Decisive evidence for Model 2"
  }

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
