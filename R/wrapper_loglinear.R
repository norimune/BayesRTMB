#' RTMB-based Log-linear analysis (Poisson regression)
#'
#' @description
#' Performs Bayesian or Frequentist log-linear analysis (Poisson regression)
#' on a contingency table or raw data.
#'
#' @param formula A formula (e.g., `~ A + B + A:B`) or a contingency table.
#' @param data A data frame (required if `formula` is used).
#' @param prior An object of class "rtmb_prior" specifying the prior distribution.
#' @param y_range Optional theoretical minimum and maximum values of the count
#'   response. If supplied with the default flat prior, `prior_weak()` is used.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Additional arguments passed to `rtmb_glm()`.
#' @return An `RTMB_Model` object.
#' @example inst/examples/ex_loglinear.R
#' @export
rtmb_loglinear <- function(formula, data, prior = prior_flat(), y_range = NULL,
                           fixed = NULL, WAIC = FALSE, ...) {

  if (is.null(prior)) prior <- prior_flat()
  if (missing(data) || is.null(data)) {
    stop("'data' must be supplied as a data frame, table, or matrix.", call. = FALSE)
  }

  # 1. Data Preparation
  # Convert table or matrix to long data frame
  was_table <- FALSE
  if (inherits(data, c("table", "matrix"))) {
    data <- as.data.frame(as.table(data))
    was_table <- TRUE
  }

  # Identify variables and frequency column
  all_vars <- all.vars(formula)
  response_var <- NULL

  # Check if formula has a response
  if (length(formula) == 3) {
    response_var <- as.character(formula[[2]])
  }

  # If no response is specified
  if (is.null(response_var)) {
    if (was_table && "Freq" %in% names(data)) {
      # For tables/matrices, we automatically use the Freq column
      formula <- update(formula, Freq ~ .)
      response_var <- "Freq"
    } else {
      # For raw data frames, we aggregate
      # Check if variables exist in data
      missing_vars <- setdiff(all_vars, names(data))
      if (length(missing_vars) > 0) {
        stop(sprintf("Variables not found in data: %s", paste(missing_vars, collapse = ", ")))
      }

      # Aggregate raw observations into counts
      data_counts <- as.data.frame(table(data[all_vars]))
      # Update formula to include 'Freq' as response
      formula <- update(formula, Freq ~ .)
      data <- data_counts
      response_var <- "Freq"
    }
  }


  # Ensure response is numeric (counts)
  data[[response_var]] <- as.numeric(data[[response_var]])

  # Call the engine
  # We use rtmb_glmer which handles formulas, random effects, and priors.
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && identical(prior$type, "flat")) {
    prior <- prior_weak()
  }
  prior <- .validate_prior_type(
    prior,
    allowed = c("flat", "normal", "weak", "rhs", "ssp"),
    context = "rtmb_loglinear()"
  )

  # For log-linear models, use weak Poisson defaults when requested.
  if (inherits(prior, "rtmb_prior") && identical(prior$type, "weak")) {
    # Stan-style defaults for log-linear (Poisson)
    if (is.null(prior$sd_ratio)) prior$sd_ratio <- 5.0 # For Intercept (alpha_prior_sd)
    if (is.null(prior$max_beta)) prior$max_beta <- 2.5 # For coefficients
  }

  # Call the engine
  # We use rtmb_glmer which handles formulas, random effects, and priors.

  # 2. Determine centering and predictors (to keep generate block clean in print_code)
  use_centering <- prior$type %in% c("normal", "weak", "rhs", "ssp", "jzs")
  is_centered <- use_centering

  # Identify if there are predictors
  X_tmp <- stats::model.matrix(nobars(formula), data[1, , drop = FALSE])
  K <- ncol(X_tmp) - (if (attr(stats::terms(nobars(formula)), "intercept")) 1 else 0)

  if (is_centered) {
    int_name <- as.name("Intercept_c")
    x_name <- as.name("X_c")
  } else {
    int_name <- as.name("Intercept")
    x_name <- as.name("X")
  }

  if (K > 0) {
    gen_block <- bquote({
      eta <- .(int_name) + .(x_name) %*% b
      mu <- exp(eta)
      list(mu = mu)
    })
  } else {
    gen_block <- bquote({
      eta <- rep(.(int_name), N)
      mu <- exp(eta)
      list(mu = mu)
    })
  }

  res <- rtmb_glmer(
    formula = formula,
    data = data,
    family = "poisson",
    prior = prior,
    y_range = y_range,
    generate = gen_block, fixed = fixed, WAIC = WAIC,
    ...
  )

  # Tag as loglinear
  res$type <- "loglinear"
  res$extra$source <- "wrapper"
  # res$extra$prior_type is maintained from rtmb_glmer
  res$extra$marginal <- if (is_centered) "Intercept_c" else "Intercept"
  if (K > 0) res$extra$marginal <- c(res$extra$marginal, "b")
  res$extra$obs_Y <- data[[response_var]]
  fit <- res

  # Tag as rtmb_loglinear for potential custom methods
  class(fit) <- c("rtmb_loglinear", class(fit))

  return(fit)
}
