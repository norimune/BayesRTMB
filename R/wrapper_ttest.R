#' RTMB-based Bayesian two-sample t-test wrapper function
#'
#' @description
#' Performs a Bayesian or Frequentist two-sample t-test using RTMB.
#'
#' @param x Numeric vector of responses for group 1, a formula (e.g., `y ~ group`), or a column name (unquoted) if `data` is provided.
#' @param y Numeric vector of responses for group 2, or a column name (unquoted) if `data` is provided. Required if `x` is not a formula.
#' @param data Data frame containing the variables.
#' @param r Numeric; Cauchy prior scale for the effect size (delta). Default is 0.707.
#' @param paired Logical; whether to perform a paired t-test.
#' @param ID Character; name of the ID variable for paired t-tests (required for formula input with paired = TRUE).
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Required when using weakly informative priors.
#' @param prior An object of class `"rtmb_prior"`.
#' Use `prior_flat()` for no prior, `prior_normal()` for default normal/exponential priors,
#' `prior_jzs()` for JZS t-test priors, or `prior_weak()` for weakly
#' informative Bayesian inference.
#' Default is `prior_flat()`.
#' @param init List of initial values.
#' @param var.equal Logical; whether to assume equal variances. Default is TRUE.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param missing Missing value handling strategy: "listwise".
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Reserved; unused arguments are rejected.
#' @return An `RTMB_Model` object.
#'
#' @details
#' For classic inference, heteroscedastic two-sample t-tests use the same
#' RTMB Satterthwaite machinery as `optimize(marginal = ..., df_method = "satterthwaite")`.
#' The result is model-based and reproducible from the printed model code.
#' This corresponds to the Welch-type unequal-variance t-test, but the
#' degrees of freedom are computed by the package's internal Satterthwaite
#' procedure rather than by a separate closed-form formula.
#' When `prior_jzs()` is combined with `var.equal = FALSE`, BayesRTMB uses
#' a Welch-style JZS extension: the effect size `delta` is an explicit
#' parameter with a Cauchy prior, and the group mean difference is scaled by
#' the root-mean-square of the two group standard deviations.
#'
#' @example inst/examples/ex_ttest.R
#' @export
rtmb_ttest <- function(x, y = NULL, data = NULL, r = 0.707,
                       paired = FALSE, ID = NULL,
                       y_range = NULL,
                       prior = prior_flat(),
                       init = NULL, fixed = NULL,
                       var.equal = TRUE, missing = c("listwise", "fiml"),
                       WAIC = FALSE, ...) {

  .check_unused_dots(..., .fn = "rtmb_ttest()")
  missing <- match.arg(missing)
  x_expr <- substitute(x)
  y_expr <- substitute(y)

  # --- 1. Data Extraction ---
  is_formula <- is.call(x_expr) && (x_expr[[1]] == quote(`~`))
  formula_input <- is_formula
  raw_response <- NULL
  raw_group_idx <- NULL
  raw_id <- NULL
  if (is_formula) {
    x <- eval(x_expr, parent.frame())
    mf <- if (is.null(data)) model.frame(x, parent.frame()) else model.frame(x, data)
    response <- mf[[1]]; group <- as.factor(mf[[2]])
    levs <- levels(group)
    if (length(levs) != 2) {
      stop(
        sprintf("The grouping variable must have exactly 2 levels, but found %d: %s",
                length(levs), paste(levs, collapse = ", ")),
        call. = FALSE
      )
    }
    x_label <- levs[1]; y_label <- levs[2]
    
    if (!is.numeric(response) && !is.logical(response)) {
      stop("The response variable must be numeric. Character or factor variables are not supported.", call. = FALSE)
    }
    raw_response <- as.numeric(response)
    raw_group_idx <- as.integer(group)
    if (paired) {
      if (is.null(ID)) {
        stop("'ID' is required for paired t-tests with formula input (e.g., ID = 'subject_id').", call. = FALSE)
      }
      id_val <- if (!is.null(data)) data[[ID]] else eval(substitute(ID), parent.frame())
      if (is.null(id_val)) {
        stop(sprintf("ID variable '%s' not found in data.", as.character(substitute(ID))), call. = FALSE)
      }
      if (length(id_val) != nrow(mf)) {
        mf_rows <- suppressWarnings(as.integer(rownames(mf)))
        if (length(mf_rows) == nrow(mf) && !anyNA(mf_rows) && max(mf_rows) <= length(id_val)) {
          id_val <- id_val[mf_rows]
        } else {
          stop("Could not align 'ID' with the rows used in the paired t-test.", call. = FALSE)
        }
      }
      raw_id <- id_val
      d1 <- mf[group == levs[1], ]; d1$id <- id_val[group == levs[1]]
      d2 <- mf[group == levs[2], ]; d2$id <- id_val[group == levs[2]]
      common_ids <- intersect(d1$id, d2$id)
      Y1 <- d1[match(common_ids, d1$id), 1]; Y2 <- d2[match(common_ids, d2$id), 1]
    } else {
      Y1 <- as.numeric(response[group == levs[1]])
      Y2 <- as.numeric(response[group == levs[2]])
      if (missing == "listwise") {
         Y1 <- as.numeric(na.omit(Y1))
         Y2 <- as.numeric(na.omit(Y2))
      }
    }
  } else {
    x_label <- deparse(x_expr); y_label <- deparse(y_expr)
    if (identical(y_expr, quote(NULL))) {
      stop("'y' is required when 'x' is not a formula.", call. = FALSE)
    }
    if (!is.null(data)) {
      Y1 <- if (x_label %in% names(data)) data[[x_label]] else eval(x_expr, parent.frame())
      Y2 <- if (y_label %in% names(data)) data[[y_label]] else eval(y_expr, parent.frame())
    } else {
      Y1 <- eval(x_expr, parent.frame()); Y2 <- eval(y_expr, parent.frame())
    }
    
    if ((!is.numeric(Y1) && !is.logical(Y1)) || (!is.numeric(Y2) && !is.logical(Y2))) {
      stop("Both variables (x and y) must be numeric. Character or factor variables are not supported.", call. = FALSE)
    }
    if (paired) {
       if (length(Y1) != length(Y2)) stop("For paired t-tests without a formula, x and y must have the same length.", call. = FALSE)
       if (missing == "listwise") {
         valid_idx <- !is.na(Y1) & !is.na(Y2)
         Y1 <- Y1[valid_idx]
         Y2 <- Y2[valid_idx]
       }
    } else {
       if (missing == "listwise") {
         Y1 <- as.numeric(na.omit(Y1)); Y2 <- as.numeric(na.omit(Y2))
       }
    }
  }
  if (paired) {
    if (length(Y1) < 2L || length(Y2) < 2L) {
      stop("Paired t-tests require at least two complete pairs.", call. = FALSE)
    }
  } else {
    if (length(Y1) < 2L || length(Y2) < 2L) {
      stop("Two-sample t-tests require at least two observations in each group.", call. = FALSE)
    }
  }
  levs <- c(x_label, y_label)

  # --- 2. Prior Setup ---
  if (is.null(prior)) prior <- prior_flat()

  # Automatically switch to prior_weak() if y_range is provided and prior is default flat
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && identical(prior$type, "flat")) {
    prior <- prior_weak()
  }
  prior <- .validate_prior_type(
    prior,
    allowed = c("flat", "normal", "jzs", "weak"),
    context = "rtmb_ttest()"
  )
  prior_type <- prior$type
  is_jzs <- prior_type == "jzs"
  is_weak <- prior_type == "weak"
  is_normal <- prior_type == "normal"
  use_delta_param <- is_jzs
  if (is_jzs) r <- prior$r
  
  if (is_normal) {
    if (is.null(prior$mu_sd) && !is.null(prior$Intercept_sd)) {
      prior$mu_sd <- prior$Intercept_sd
    }
  }

  # --- 3. AST Construction ---
  setup_exprs <- list(as.name("{"))
  if (formula_input) {
    if (paired) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(d1 <- data.frame(Y = Y[G == 1], ID = ID[G == 1]))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(d2 <- data.frame(Y = Y[G == 2], ID = ID[G == 2]))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(common_ids <- intersect(d1$ID, d2$ID))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(Y1 <- d1$Y[match(common_ids, d1$ID)])
      setup_exprs[[length(setup_exprs) + 1]] <- quote(Y2 <- d2$Y[match(common_ids, d2$ID)])
      setup_exprs[[length(setup_exprs) + 1]] <- quote(diffs <- Y1 - Y2)
    } else {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(Y1 <- as.numeric(na.omit(Y[G == 1])))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(Y2 <- as.numeric(na.omit(Y[G == 2])))
    }
  }
  if (is_weak) {
    if (is.null(y_range)) {
      stop("Specifying 'y_range' (theoretical minimum and maximum values of the response) is required when using prior_weak().", call. = FALSE)
    }
    
    # Common setup values
    sd_ratio_val <- prior$sd_ratio
    max_beta_val <- prior$max_beta
    
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(half_d_y <- .(diff(y_range) / 2))
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(base_scale <- half_d_y * .(sd_ratio_val))
    setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- half_d_y)
    if (paired) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(x_sd <- 0.5)
    } else {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(n1 <- length(Y1))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(n2 <- length(Y2))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(p1 <- n1 / (n1 + n2))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(x_sd <- sqrt(p1 * (1 - p1)))
    }
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(diff_prior_sd <- .(max_beta_val) * base_scale / x_sd)
    setup_exprs[[length(setup_exprs) + 1]] <- bquote(mid_y <- .(if (paired) 0 else mean(y_range)))
    setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate_weak <- 1.0 / base_scale)
  }
  setup_ast <- as.call(setup_exprs)

  if (paired) {
    dat <- if (formula_input) list(Y = raw_response, G = raw_group_idx, ID = raw_id, r = r) else list(diffs = Y1 - Y2, r = r)
    param_ast <- if (use_delta_param) {
      quote({ sd_diff = Dim(1, lower = 0); delta = Dim(1) })
    } else {
      quote({ sd_diff = Dim(1, lower = 0); diff = Dim(1) })
    }
    tran_ast <- if (use_delta_param) {
      quote({ diff <- delta * sd_diff })
    } else {
      quote({ delta <- diff / sd_diff })
    }
    model_body <- list(quote(diffs ~ normal(diff, sd_diff)))
    if (is_jzs) model_body[[length(model_body)+1]] <- quote(delta ~ cauchy(0, r))
    if (is_weak) {
      model_body[[length(model_body)+1]] <- quote(diff ~ normal(0, diff_prior_sd))
      model_body[[length(model_body)+1]] <- quote(sd_diff ~ exponential(sigma_rate_weak))
    }
    if (is_normal) {
      if (!is.null(prior$b_sd)) {
        model_body[[length(model_body)+1]] <- bquote(diff ~ normal(0, .(prior$b_sd)))
      }
      if (!is.null(prior$sigma_rate)) {
        model_body[[length(model_body)+1]] <- bquote(sd_diff ~ exponential(.(prior$sigma_rate)))
      }
    }
    view_vars <- c("diff", "delta", "sd_diff")
    marginal <- if (use_delta_param) "delta" else "diff"
    waic_ast <- quote({
      log_lik <- normal_lpdf(diffs, diff, sd_diff, sum = FALSE)
    })
  } else {
    dat <- if (formula_input) list(Y = raw_response, G = raw_group_idx, r = r) else list(Y1 = Y1, Y2 = Y2, r = r)
    
    # Structure for independent t-tests
    if (!var.equal) {
      if (use_delta_param) {
        param_ast <- quote({ total_mean = Dim(1); sd = Dim(2, lower = 0); delta = Dim(1) })
        tran_ast <- quote({
          sd_pooled <- sqrt((sd[1]^2 + sd[2]^2)/2)
          diff <- delta * sd_pooled
          mean0 <- total_mean + diff/2
          mean1 <- total_mean - diff/2
        })
        view_vars <- c("diff", "delta", "total_mean", "mean0", "mean1", "sd")
        marginal <- c("total_mean", "delta")
      } else {
        # Heteroscedastic (Welch style): Use group means as parameters to trigger Satterthwaite on diff.
        param_ast <- quote({ mean0 = Dim(1); mean1 = Dim(1); sd = Dim(2, lower = 0) })
        tran_ast <- quote({
          diff <- mean0 - mean1
          total_mean <- (mean0 + mean1)/2
          sd_pooled <- sqrt((sd[1]^2 + sd[2]^2)/2)
          delta <- diff / sd_pooled
        })
        view_vars <- c("diff", "delta", "mean0", "mean1", "sd")
        marginal <- c("mean0", "mean1")
      }
      model_body <- list(quote(Y1 ~ normal(mean0, sd[1])), quote(Y2 ~ normal(mean1, sd[2])))
      waic_ast <- quote({
        log_lik <- c(
          normal_lpdf(Y1, mean0, sd[1], sum = FALSE),
          normal_lpdf(Y2, mean1, sd[2], sum = FALSE)
        )
      })
    } else {
      # Homoscedastic (Standard pooled): Use total_mean and diff as parameters to fix DF at N-K
      param_ast <- if (use_delta_param) {
        quote({ total_mean = Dim(1); sd = Dim(1, lower = 0); delta = Dim(1) })
      } else {
        quote({ total_mean = Dim(1); sd = Dim(1, lower = 0); diff = Dim(1) })
      }
      tran_ast <- if (use_delta_param) {
        quote({ diff <- delta * sd; mean0 <- total_mean + diff/2; mean1 <- total_mean - diff/2 })
      } else {
        quote({ delta <- diff / sd; mean0 <- total_mean + diff/2; mean1 <- total_mean - diff/2 })
      }
      model_body <- list(quote(Y1 ~ normal(mean0, sd)), quote(Y2 ~ normal(mean1, sd)))
      view_vars <- c("diff", "delta", "total_mean", "sd")
      marginal <- if (use_delta_param) c("total_mean", "delta") else c("total_mean", "diff")
      waic_ast <- quote({
        log_lik <- c(
          normal_lpdf(Y1, mean0, sd, sum = FALSE),
          normal_lpdf(Y2, mean1, sd, sum = FALSE)
        )
      })
    }
    
    if (is_jzs) {
       model_body[[length(model_body)+1]] <- quote(delta ~ cauchy(0, r))
    }
    if (is_weak) {
      if (!var.equal) {
        model_body[[length(model_body)+1]] <- quote(mean0 ~ normal(mid_y, alpha_prior_sd))
        model_body[[length(model_body)+1]] <- quote(mean1 ~ normal(mid_y, alpha_prior_sd))
      } else {
        model_body[[length(model_body)+1]] <- quote(total_mean ~ normal(mid_y, alpha_prior_sd))
        model_body[[length(model_body)+1]] <- quote(diff ~ normal(0, diff_prior_sd))
      }
      model_body[[length(model_body)+1]] <- quote(sd ~ exponential(sigma_rate_weak))
    }
    if (is_normal) {
      if (!var.equal) {
        mu_sd <- if (!is.null(prior$mu_sd)) prior$mu_sd else prior$Intercept_sd
        if (!is.null(mu_sd)) {
          model_body[[length(model_body)+1]] <- bquote(mean0 ~ normal(0, .(mu_sd)))
          model_body[[length(model_body)+1]] <- bquote(mean1 ~ normal(0, .(mu_sd)))
        }
      } else {
        mu_sd <- if (!is.null(prior$mu_sd)) prior$mu_sd else prior$Intercept_sd
        if (!is.null(mu_sd)) model_body[[length(model_body)+1]] <- bquote(total_mean ~ normal(0, .(mu_sd)))
        if (!is.null(prior$b_sd)) model_body[[length(model_body)+1]] <- bquote(diff ~ normal(0, .(prior$b_sd)))
      }
      if (!is.null(prior$sigma_rate)) model_body[[length(model_body)+1]] <- bquote(sd ~ exponential(.(prior$sigma_rate)))
    }
  }

  model_ast <- as.call(c(list(as.name("{")), model_body))
  code_obj <- if (isTRUE(WAIC)) {
    waic_ast <- .rtmb_waic_generate_ast(NULL, waic_ast)
    eval(substitute(
      rtmb_code(setup = S, parameters = P, transform = T, model = M, generate = G),
      list(S = setup_ast, P = param_ast, T = tran_ast, M = model_ast, G = waic_ast)
    ))
  } else {
    eval(substitute(
      rtmb_code(setup = S, parameters = P, transform = T, model = M),
      list(S = setup_ast, P = param_ast, T = tran_ast, M = model_ast)
    ))
  }

  tmp_env <- list2env(dat)
  ordered_data <- env_to_ordered_list(tmp_env, dat, setup_ast)
  obj <- rtmb_model(data = ordered_data, code = code_obj, fixed = fixed, view = view_vars)

  obj$type <- "ttest"
  obj$extra <- list(
    source = "wrapper",
    prior_type = prior$type,
    marginal = marginal,
    levs = levs,
    paired = paired,
    var_equal = var.equal
  )

  return(obj)
}
