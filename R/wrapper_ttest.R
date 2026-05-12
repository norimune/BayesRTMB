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
#' @param prior An object of class "rtmb_prior" (e.g., `prior_uniform()`, `prior_jzs()`, or `prior_weak()`).
#' @param init List of initial values.
#' @param null Character string specifying the target parameter for the null model (e.g., "delta").
#' @param var.equal Logical; whether to assume equal variances. Default is TRUE.
#' @param ... Additional arguments.
#' @return An `RTMB_Model` object.
#' @export
rtmb_ttest <- function(x, y = NULL, data = NULL, r = 0.707,
                       paired = FALSE, ID = NULL,
                       y_range = NULL,
                       prior = prior_uniform(),
                       init = NULL, fixed = NULL, null = NULL,
                       var.equal = TRUE, ...) {

  x_expr <- substitute(x)
  y_expr <- substitute(y)

  # --- 1. Data Extraction ---
  is_formula <- is.call(x_expr) && (x_expr[[1]] == quote(`~`))
  if (is_formula) {
    x <- eval(x_expr, parent.frame())
    mf <- if (is.null(data)) model.frame(x, parent.frame()) else model.frame(x, data)
    response <- mf[[1]]; group <- as.factor(mf[[2]])
    levs <- levels(group); x_label <- levs[1]; y_label <- levs[2]
    if (paired) {
      id_val <- if (!is.null(data)) data[[ID]] else eval(substitute(ID), parent.frame())
      d1 <- mf[group == levs[1], ]; d1$id <- id_val[group == levs[1]]
      d2 <- mf[group == levs[2], ]; d2$id <- id_val[group == levs[2]]
      common_ids <- intersect(d1$id, d2$id)
      Y1 <- d1[match(common_ids, d1$id), 1]; Y2 <- d2[match(common_ids, d2$id), 1]
    } else {
      Y1 <- as.numeric(na.omit(response[group == levs[1]]))
      Y2 <- as.numeric(na.omit(response[group == levs[2]]))
    }
  } else {
    x_label <- deparse(x_expr); y_label <- deparse(y_expr)
    if (!is.null(data)) {
      Y1 <- if (x_label %in% names(data)) data[[x_label]] else eval(x_expr, parent.frame())
      Y2 <- if (y_label %in% names(data)) data[[y_label]] else eval(y_expr, parent.frame())
    } else {
      Y1 <- eval(x_expr, parent.frame()); Y2 <- eval(y_expr, parent.frame())
    }
    Y1 <- as.numeric(na.omit(Y1)); Y2 <- as.numeric(na.omit(Y2))
  }
  levs <- c(x_label, y_label)

  # --- 2. Prior Setup ---
  if (is.null(prior)) prior <- prior_uniform()

  # Automatically switch to prior_weak() if y_range is provided and prior is default uniform
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && prior$type == "uniform") {
    prior <- prior_weak()
  }
  prior_type <- prior$type
  is_jzs <- prior_type == "jzs"
  is_weak <- prior_type == "weak"
  use_delta_param <- is_jzs
  if (is_jzs) r <- prior$r

  # --- 3. AST Construction ---
  setup_ast <- quote({})
  if (is_weak) {
    if (is.null(y_range)) {
      stop("Specifying 'y_range' (theoretical minimum and maximum values of the response) is required when using prior_weak().", call. = FALSE)
    }
    
    # Common setup values
    sd_ratio_val <- prior$sd_ratio
    max_beta_val <- prior$max_beta
    
    setup_ast <- bquote({
      half_d_y <- .(diff(y_range) / 2)
      base_scale <- half_d_y * .(sd_ratio_val)
      alpha_prior_sd <- half_d_y
      
      # Determine x_sd based on design
      .(if (paired) {
        quote(x_sd <- 0.5) # Balanced Condition effect (Condition 0/1)
      } else {
        quote({
          n1 <- length(Y1)
          n2 <- length(Y2)
          p1 <- n1 / (n1 + n2)
          x_sd <- sqrt(p1 * (1 - p1))
        })
      })
      
      diff_prior_sd <- .(max_beta_val) * base_scale / x_sd
      mid_y <- .(if(paired) 0 else mean(y_range))
      sigma_rate_weak <- 1.0 / base_scale
    })
  }

  if (paired) {
    dat <- list(diffs = Y1 - Y2, r = r)
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
    view_vars <- c("diff", "delta", "sd_diff")
    marginal <- if (use_delta_param) "delta" else "diff"
  } else {
    dat <- list(Y1 = Y1, Y2 = Y2, r = r)
    
    # Structure for independent t-tests
    if (!var.equal) {
      # Heteroscedastic (Welch style): Use group means as parameters to trigger Satterthwaite on diff
      param_ast <- quote({ mean0 = Dim(1); mean1 = Dim(1); sd = Dim(2, lower = 0) })
      tran_ast <- quote({ 
        diff <- mean0 - mean1
        mean <- (mean0 + mean1)/2
        sd_pooled <- sqrt((sd[1]^2 + sd[2]^2)/2)
        delta <- diff / sd_pooled
      })
      model_body <- list(quote(Y1 ~ normal(mean0, sd[1])), quote(Y2 ~ normal(mean1, sd[2])))
      view_vars <- c("diff", "delta", "mean0", "mean1", "sd")
      marginal <- c("mean0", "mean1")
    } else {
      # Homoscedastic (Standard pooled): Use mean and diff as parameters to fix DF at N-K
      param_ast <- if (use_delta_param) {
        quote({ mean = Dim(1); sd = Dim(1, lower = 0); delta = Dim(1) })
      } else {
        quote({ mean = Dim(1); sd = Dim(1, lower = 0); diff = Dim(1) })
      }
      tran_ast <- if (use_delta_param) {
        quote({ diff <- delta * sd; mean0 <- mean + diff/2; mean1 <- mean - diff/2 })
      } else {
        quote({ delta <- diff / sd; mean0 <- mean + diff/2; mean1 <- mean - diff/2 })
      }
      model_body <- list(quote(Y1 ~ normal(mean0, sd)), quote(Y2 ~ normal(mean1, sd)))
      view_vars <- c("diff", "delta", "mean", "sd")
      marginal <- if (use_delta_param) c("mean", "delta") else c("mean", "diff")
    }
    
    if (is_jzs) {
       if (!var.equal) {
          model_body[[length(model_body)+1]] <- quote((mean0 - mean1)/sqrt((sd[1]^2 + sd[2]^2)/2) ~ cauchy(0, r))
       } else if (var.equal) {
          model_body[[length(model_body)+1]] <- quote(delta ~ cauchy(0, r))
       }
    }
    if (is_weak) {
      if (!var.equal) {
        model_body[[length(model_body)+1]] <- quote(mean0 ~ normal(mid_y, alpha_prior_sd))
        model_body[[length(model_body)+1]] <- quote(mean1 ~ normal(mid_y, alpha_prior_sd))
      } else {
        model_body[[length(model_body)+1]] <- quote(mean ~ normal(mid_y, alpha_prior_sd))
        model_body[[length(model_body)+1]] <- quote(diff ~ normal(0, diff_prior_sd))
      }
      model_body[[length(model_body)+1]] <- quote(sd ~ exponential(sigma_rate_weak))
    }
  }

  model_ast <- as.call(c(list(as.name("{")), model_body))
  code_obj <- eval(substitute(
    rtmb_code(setup = S, parameters = P, transform = T, model = M),
    list(S = setup_ast, P = param_ast, T = tran_ast, M = model_ast)
  ))

  tmp_env <- list2env(dat)
  ordered_data <- env_to_ordered_list(tmp_env, dat, setup_ast)
  obj <- rtmb_model(data = ordered_data, code = code_obj, fixed = fixed, view = view_vars)

  obj$type <- "ttest"
  obj$extra$marginal <- marginal
  obj$extra$levs <- levs
  if (!is.null(null)) obj <- obj$null_model(target = null)

  return(obj)
}
