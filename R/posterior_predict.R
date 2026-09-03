#' Posterior predictive simulation
#'
#' Generate replicated outcomes from an estimated BayesRTMB model. Regression
#' wrappers provide an automatic simulator. For a custom model, supply an
#' `rtmb_code()` object containing a `generate` block that creates and reports
#' one replicated outcome vector. MCMC and variational fits use posterior draws. A MAP fit uses
#' sampling-based uncertainty when available and otherwise conditions on its
#' point estimate.
#'
#' @param object A BayesRTMB fit object.
#' @param ... Arguments passed to the fit object's `posterior_predict()` method.
#'
#' @return A numeric matrix with posterior predictive draws in rows and
#'   observations in columns.
#' @export
posterior_predict <- function(object, ...) {
  UseMethod("posterior_predict")
}

#' @rdname posterior_predict
#' @export
posterior_predict.RTMB_Fit_Base <- function(object, ...) {
  object$posterior_predict(...)
}

#' Posterior predictive checks
#'
#' Compare observed data with replicated outcomes from an estimated BayesRTMB
#' model. Density overlays are used for continuous outcomes and binned bar plots
#' for discrete outcomes. When `type = "auto"`, the display is selected from the
#' likelihood's `_lpdf` or `_lpmf` implementation.
#'
#' Supplying `x` to the fit object's `pp_check()` method switches to a
#' scatter-based check. Observed outcomes are compared with posterior predictive
#' means and 95% predictive intervals along the selected predictor. The reserved
#' value `x = ".fitted"` uses the posterior predictive mean on the horizontal
#' axis for a model-wide calibration check.
#'
#' @param object A BayesRTMB fit object.
#' @param ... Arguments passed to the fit object's `pp_check()` method or to
#'   the plotting method.
#'
#' @return An object of class `rtmb_pp_check`, returned invisibly after plotting.
#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}

#' @rdname pp_check
#' @export
pp_check.RTMB_Fit_Base <- function(object, ...) {
  object$pp_check(...)
}

.rtmb_with_seed <- function(seed, code) {
  if (is.null(seed)) return(code())

  seed <- as.integer(seed)[1L]
  if (is.na(seed)) stop("'seed' must be a finite integer.", call. = FALSE)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  code()
}

.rtmb_resolve_generate_ast <- function(code_expr, env) {
  if (is.null(code_expr)) return(NULL)

  code <- code_expr
  if (is.name(code) ||
      (is.call(code) && !identical(code[[1L]], as.name("{")) &&
       !identical(code[[1L]], as.name("rtmb_code")))) {
    code <- tryCatch(eval(code, envir = env), error = function(e) e)
    if (inherits(code, "error")) {
      stop("Could not evaluate 'code': ", conditionMessage(code), call. = FALSE)
    }
  } else if (is.call(code) && identical(code[[1L]], as.name("rtmb_code"))) {
    code <- eval(code, envir = env)
  }

  if (is.null(code)) return(NULL)

  if (is.list(code) && "generate" %in% names(code)) {
    code <- code$generate
  }

  if (!is.call(code) || !identical(code[[1L]], as.name("{"))) {
    stop(
      "'code' must be rtmb_code(generate = { ... }), a generate block, or an object containing a generate block.",
      call. = FALSE
    )
  }

  code
}

.rtmb_compile_generate <- function(model, gen_ast, code_env) {
  gen_fn <- eval(bquote(transform_code(.(gen_ast))), envir = environment())

  setup_values <- model$code$setup_env
  if (is.environment(setup_values)) setup_values <- as.list(setup_values, all.names = TRUE)
  if (!is.list(setup_values)) setup_values <- list()

  eval_env <- list2env(setup_values, parent = code_env)
  environment(gen_fn) <- eval_env
  gen_fn
}

.rtmb_parameter_draws <- function(fit, draws, chains = NULL, best_chains = NULL) {
  if (is.null(draws)) {
    draws_requested <- NULL
  } else {
    draws_requested <- as.integer(draws)[1L]
    if (is.na(draws_requested) || draws_requested < 1L) {
      stop("'draws' must be NULL or a positive integer.", call. = FALSE)
    }
  }

  if (!is.function(fit$draws)) {
    if (inherits(fit, "Classic_Fit")) {
      stop(
        "posterior_predict() is not available for Classic_Fit objects because ",
        "classic() does not store posterior draws. Use the fit returned by ",
        "sample(), optimize(), or variational() instead.",
        call. = FALSE
      )
    }
    stop(
      "This fit object does not provide a draws() method. If it was created ",
      "with an older BayesRTMB version, run upgrade_fit(fit) or refit the model.",
      call. = FALSE
    )
  }

  draw_array <- fit$draws(
    chains = chains,
    best_chains = best_chains,
    inc_random = TRUE,
    inc_transform = FALSE,
    inc_generate = FALSE
  )

  d <- dim(draw_array)
  if (length(d) != 3L || d[1L] < 1L || d[2L] < 1L) {
    stop("No posterior draws are available.", call. = FALSE)
  }

  draw_matrix <- matrix(draw_array, nrow = d[1L] * d[2L], ncol = d[3L])
  colnames(draw_matrix) <- dimnames(draw_array)[[3L]]
  if (is.null(colnames(draw_matrix))) {
    stop("Posterior draws do not have parameter names.", call. = FALSE)
  }

  par_columns <- vector("list", length(fit$model$par_list))
  names(par_columns) <- names(fit$model$par_list)

  for (name in names(fit$model$par_list)) {
    info <- fit$model$par_list[[name]]
    expected <- generate_flat_names(name, info$dim, fit$model$par_names[[name]])
    idx <- match(expected, colnames(draw_matrix))

    if (anyNA(idx)) {
      idx_fallback <- grep(paste0("^", name, "(\\[|$)"), colnames(draw_matrix))
      if (length(idx_fallback) == info$length) idx <- idx_fallback
    }

    if (anyNA(idx) || length(idx) != info$length) {
      stop(
        "Could not reconstruct posterior draws for parameter '", name,
        "'. Refit the model with the current BayesRTMB version.",
        call. = FALSE
      )
    }
    par_columns[[name]] <- idx
  }

  total <- nrow(draw_matrix)
  if (is.null(draws_requested)) {
    selected <- seq_len(total)
  } else if (draws_requested <= total) {
    selected <- sample.int(total, draws_requested, replace = FALSE)
  } else {
    selected <- sample.int(total, draws_requested, replace = TRUE)
  }

  list(
    matrix = draw_matrix,
    selected = selected,
    columns = par_columns
  )
}

.rtmb_parameter_state <- function(fit, draw_info, row) {
  state <- vector("list", length(fit$model$par_list))
  names(state) <- names(fit$model$par_list)

  for (name in names(state)) {
    info <- fit$model$par_list[[name]]
    value <- as.numeric(draw_info$matrix[row, draw_info$columns[[name]], drop = TRUE])
    if (length(info$dim) > 1L) dim(value) <- info$dim
    state[[name]] <- value
  }

  if (!is.null(fit$model$transform)) {
    transformed <- fit$model$transform(fit$model$data, state)
    if (!is.null(transformed)) {
      if (!is.list(transformed)) {
        stop("The model transform block did not return a list.", call. = FALSE)
      }
      state[names(transformed)] <- transformed
    }
  }

  state
}

.rtmb_glmer_eta <- function(data, state, spec, random) {
  n <- length(data[[spec$response]])
  k <- spec$K

  if (identical(spec$family, "sequential")) {
    if (k > 0L) {
      eta <- data$X %*% state$b
    } else {
      eta <- matrix(0, nrow = n, ncol = spec$num_categories - 1L)
    }
  } else if (isTRUE(spec$has_intercept)) {
    if (isTRUE(spec$use_centering)) {
      eta <- if (k > 0L) {
        as.numeric(state$Intercept_c) + data$X_c %*% state$b
      } else {
        rep(as.numeric(state$Intercept_c), n)
      }
    } else {
      eta <- if (k > 0L) {
        as.numeric(state$Intercept) + data$X %*% state$b
      } else {
        rep(as.numeric(state$Intercept), n)
      }
    }
  } else {
    eta <- if (k > 0L) data$X %*% state$b else rep(0, n)
  }

  if (!identical(spec$family, "sequential")) eta <- as.numeric(eta)

  if (length(spec$random_terms) > 0L && !identical(random, "population")) {
    for (term in spec$random_terms) {
      z <- data[[term$z_name]]
      group <- data[[term$group_name]]
      sd_re <- as.numeric(state[[term$sd_name]])

      if (identical(random, "conditional")) {
        r_re <- state[[term$effect_name]]
      } else {
        if (term$num_ranef == 1L) {
          r_re <- stats::rnorm(term$num_groups)
        } else {
          independent <- matrix(
            stats::rnorm(term$num_groups * term$num_ranef),
            nrow = term$num_groups,
            ncol = term$num_ranef
          )
          cf_corr <- state[[term$corr_name]]
          r_re <- independent %*% t(cf_corr)
        }
      }

      if (term$num_ranef == 1L) {
        contribution <- z[, 1L] * as.numeric(r_re)[group] * sd_re[1L]
      } else {
        r_rows <- r_re[group, , drop = FALSE]
        contribution <- rowSums(z * sweep(r_rows, 2L, sd_re, `*`))
      }

      if (identical(spec$family, "sequential")) {
        eta <- eta + matrix(contribution, nrow = n, ncol = ncol(eta))
      } else {
        eta <- eta + contribution
      }
    }
  }

  if (isTRUE(spec$has_offset)) {
    if (identical(spec$family, "sequential")) {
      eta <- eta + matrix(data$offset, nrow = n, ncol = ncol(eta))
    } else {
      eta <- eta + data$offset
    }
  }

  eta
}

.rtmb_ordered_rng <- function(eta, cutpoints) {
  eta <- as.numeric(eta)
  cutpoints <- as.numeric(cutpoints)
  n <- length(eta)
  k <- length(cutpoints) + 1L
  out <- integer(n)

  for (i in seq_len(n)) {
    cdf <- stats::plogis(cutpoints - eta[i])
    prob <- c(cdf[1L], diff(cdf), 1 - cdf[length(cdf)])
    prob <- pmax(prob, 0)
    prob <- prob / sum(prob)
    out[i] <- sample.int(k, size = 1L, prob = prob)
  }

  out
}

.rtmb_sequential_rng <- function(eta, cutpoints) {
  cutpoints <- as.numeric(cutpoints)
  k <- length(cutpoints) + 1L
  eta_matrix <- if (is.matrix(eta)) eta else matrix(eta, ncol = k - 1L)
  n <- nrow(eta_matrix)
  out <- rep.int(k, n)

  for (i in seq_len(n)) {
    for (stage in seq_len(k - 1L)) {
      advance <- stats::plogis(eta_matrix[i, stage] - cutpoints[stage])
      if (stats::runif(1L) > advance) {
        out[i] <- stage
        break
      }
    }
  }

  out
}

.rtmb_residual_covariance <- function(data, state, spec, idx) {
  times <- data[[spec$resid_time_name]][idx]
  sigma <- as.numeric(state$sigma)
  sigma_obs <- if (isTRUE(spec$has_sigma_idx)) {
    sigma[data$sigma_idx[idx]]
  } else {
    rep(sigma[1L], length(idx))
  }

  nr <- length(idx)
  covariance <- matrix(0, nr, nr)
  corr_un <- if (identical(spec$resid_corr, "un")) {
    state$L_resid %*% t(state$L_resid)
  } else {
    NULL
  }

  for (r in seq_len(nr)) {
    for (cc in seq_len(nr)) {
      lag <- abs(times[r] - times[cc])
      covariance[r, cc] <- switch(
        spec$resid_corr,
        ar1 = sigma_obs[r] * sigma_obs[cc] * state$rho_resid^lag,
        cs = if (r == cc) sigma_obs[r]^2 else sigma_obs[r] * sigma_obs[cc] * state$rho_resid,
        toep = if (lag == 0) {
          sigma_obs[r] * sigma_obs[cc]
        } else if (lag <= length(state$rho_resid)) {
          sigma_obs[r] * sigma_obs[cc] * state$rho_resid[lag]
        } else {
          0
        },
        un = {
          ir <- times[r] + 1L
          ic <- times[cc] + 1L
          sigma_obs[r] * sigma_obs[cc] * corr_un[ir, ic]
        }
      )
    }
  }
  diag(covariance) <- diag(covariance) + 1e-3
  covariance
}

.rtmb_glmer_rng <- function(data, state, spec, random) {
  eta <- .rtmb_glmer_eta(data, state, spec, random)
  n <- length(data[[spec$response]])

  if (!is.null(spec$resid_corr)) {
    group <- data[[spec$resid_group_name]]
    out <- numeric(n)
    for (g in unique(group)) {
      idx <- which(group == g)
      covariance <- .rtmb_residual_covariance(data, state, spec, idx)
      out[idx] <- as.numeric(MASS::mvrnorm(1L, mu = eta[idx], Sigma = covariance))
    }
    return(out)
  }

  sigma <- if (isTRUE(spec$has_sigma_idx)) {
    as.numeric(state$sigma)[data$sigma_idx]
  } else {
    as.numeric(state$sigma)
  }

  switch(
    spec$family,
    gaussian = stats::rnorm(n, mean = eta, sd = sigma),
    lognormal = stats::rlnorm(n, meanlog = eta, sdlog = sigma),
    student_t = eta + sigma * stats::rt(n, df = as.numeric(state$nu)),
    gamma = stats::rgamma(
      n,
      shape = as.numeric(state$shape),
      rate = as.numeric(state$shape) / exp(eta)
    ),
    bernoulli = stats::rbinom(n, size = 1L, prob = stats::plogis(eta)),
    binomial = stats::rbinom(n, size = data$trials, prob = stats::plogis(eta)),
    poisson = stats::rpois(n, lambda = exp(eta)),
    neg_binomial = stats::rnbinom(n, size = as.numeric(state$phi), mu = exp(eta)),
    ordered = .rtmb_ordered_rng(eta, state$cutpoints),
    sequential = .rtmb_sequential_rng(eta, state$cutpoints),
    stop("Unsupported regression family for posterior prediction: ", spec$family, call. = FALSE)
  )
}

.rtmb_extract_prediction <- function(result, variable) {
  if (!is.list(result)) {
    if (!is.null(variable)) {
      warning("'variable' was ignored because the generate block returned a single value.", call. = FALSE)
    }
    value <- result
    selected_name <- "y_rep"
  } else {
    result_names <- names(result)
    if (is.null(result_names)) result_names <- rep("", length(result))

    if (is.null(variable)) {
      if ("y_rep" %in% result_names) {
        variable <- "y_rep"
      } else if (length(result) == 1L) {
        variable <- result_names[1L]
      } else {
        stop(
          "The generate block returned multiple quantities. Specify the replicated outcome with 'variable'.",
          call. = FALSE
        )
      }
    }

    if (!nzchar(variable) || !variable %in% result_names) {
      stop("Generated quantity '", variable, "' was not returned by 'code'.", call. = FALSE)
    }
    value <- result[[variable]]
    selected_name <- variable
  }

  if (!is.numeric(value) && !is.integer(value) && !is.logical(value)) {
    stop("The replicated outcome must be numeric, integer, or logical.", call. = FALSE)
  }

  list(value = as.numeric(value), name = selected_name)
}

.rtmb_pp_call_name <- function(call) {
  if (!is.call(call)) return(NULL)
  head <- call[[1L]]
  if (is.name(head)) return(as.character(head))
  if (is.call(head) && as.character(head[[1L]]) %in% c("::", ":::")) {
    return(as.character(head[[3L]]))
  }
  NULL
}

.rtmb_expression_symbols <- function(expr) {
  if (is.name(expr)) return(as.character(expr))
  if (!is.call(expr)) return(character(0))
  unique(unlist(lapply(as.list(expr)[-1L], .rtmb_expression_symbols), use.names = FALSE))
}

.rtmb_setup_density_functions <- function(setup) {
  found <- character(0)
  walk <- function(expr) {
    if (!is.call(expr)) return()
    call_name <- .rtmb_pp_call_name(expr)
    if (!is.null(call_name) && call_name %in% c("<-", "=") &&
        length(expr) >= 3L && is.name(expr[[2L]])) {
      name <- as.character(expr[[2L]])
      if (grepl("_(lpdf|lpmf)$", name)) found <<- c(found, name)
    }
    for (part in as.list(expr)[-1L]) walk(part)
  }
  if (!is.null(setup)) walk(setup)
  unique(found)
}

.rtmb_distribution_density <- function(distribution, model) {
  if (is.null(distribution) || length(distribution) != 1L || !nzchar(distribution)) {
    return(NULL)
  }
  discrete <- c(
    "bernoulli", "bernoulli_logit", "binomial", "binomial_logit",
    "poisson", "neg_binomial", "neg_binomial_2", "categorical",
    "categorical_logit", "multinomial", "ordered_logistic",
    "sequential_logistic"
  )
  if (distribution %in% discrete) return("lpmf")

  setup_functions <- .rtmb_setup_density_functions(model$code$setup)
  lpdf <- paste0(distribution, "_lpdf")
  lpmf <- paste0(distribution, "_lpmf")
  has_lpdf <- lpdf %in% setup_functions ||
    exists(lpdf, envir = asNamespace("BayesRTMB"), inherits = FALSE) ||
    is.function(model$data[[lpdf]])
  has_lpmf <- lpmf %in% setup_functions ||
    exists(lpmf, envir = asNamespace("BayesRTMB"), inherits = FALSE) ||
    is.function(model$data[[lpmf]])

  if (has_lpdf && !has_lpmf) return("lpdf")
  if (has_lpmf && !has_lpdf) return("lpmf")
  NULL
}

.rtmb_infer_density <- function(model) {
  spec <- model$extra$posterior_predict
  if (!is.null(spec$density)) return(spec$density)

  model_ast <- model$code$model
  if (is.null(model_ast)) return(NULL)
  data_names <- names(model$data)
  found <- character(0)

  walk <- function(expr) {
    if (!is.call(expr)) return()
    call_name <- .rtmb_pp_call_name(expr)

    if (identical(call_name, "~") && length(expr) >= 3L) {
      lhs_names <- .rtmb_expression_symbols(expr[[2L]])
      if (any(lhs_names %in% data_names)) {
        distribution <- .rtmb_pp_call_name(expr[[3L]])
        density <- .rtmb_distribution_density(distribution, model)
        if (!is.null(density)) found <<- c(found, density)
      }
    } else if (!is.null(call_name) && grepl("_(lpdf|lpmf)$", call_name)) {
      arg_names <- unique(unlist(lapply(as.list(expr)[-1L], .rtmb_expression_symbols), use.names = FALSE))
      if (any(arg_names %in% data_names)) {
        found <<- c(found, sub("^.*_(lpdf|lpmf)$", "\\1", call_name))
      }
    }

    for (part in as.list(expr)[-1L]) walk(part)
  }
  walk(model_ast)

  found <- unique(found)
  if (length(found) == 1L) found else NULL
}

.rtmb_observed <- function(fit, observed, spec) {
  if (is.null(observed)) {
    response_name <- spec$response
    if (is.null(response_name) && "Y" %in% names(fit$model$data)) response_name <- "Y"
    if (is.null(response_name)) {
      stop("Specify observed data with 'observed' for this custom model.", call. = FALSE)
    }
    observed <- fit$model$data[[response_name]]
  } else if (is.character(observed) && length(observed) == 1L &&
             observed %in% names(fit$model$data)) {
    observed <- fit$model$data[[observed]]
  }

  if (!is.numeric(observed) && !is.integer(observed) && !is.logical(observed)) {
    stop("'observed' must be numeric or the name of numeric model data.", call. = FALSE)
  }
  observed_names <- names(observed)
  observed <- as.numeric(observed)
  names(observed) <- observed_names
  observed
}

.rtmb_posterior_predict <- function(fit, code_expr = NULL, code_env = parent.frame(),
                                    variable = NULL, draws = 100L, seed = NULL,
                                    random = c("conditional", "population", "simulate"),
                                    chains = NULL, best_chains = NULL,
                                    observed = NULL) {
  random <- match.arg(random)
  spec <- fit$model$extra$posterior_predict
  gen_ast <- .rtmb_resolve_generate_ast(code_expr, code_env)

  if (is.null(gen_ast) && (is.null(spec) || !identical(spec$kind, "glmer"))) {
    stop(
      "Automatic posterior prediction is available for regression wrappers. For this model, supply 'code = rtmb_code(generate = { ... })'.",
      call. = FALSE
    )
  }

  .rtmb_with_seed(seed, function() {
    draw_info <- .rtmb_parameter_draws(fit, draws, chains, best_chains)
    gen_fn <- if (!is.null(gen_ast)) .rtmb_compile_generate(fit$model, gen_ast, code_env) else NULL
    predictions <- vector("list", length(draw_info$selected))
    prediction_name <- if (is.null(variable)) "y_rep" else variable

    for (i in seq_along(draw_info$selected)) {
      state <- .rtmb_parameter_state(fit, draw_info, draw_info$selected[i])
      if (is.null(gen_fn)) {
        predictions[[i]] <- .rtmb_glmer_rng(fit$model$data, state, spec, random)
      } else {
        generated <- gen_fn(fit$model$data, state)
        extracted <- .rtmb_extract_prediction(generated, variable)
        predictions[[i]] <- extracted$value
        prediction_name <- extracted$name
      }
    }

    lengths <- lengths(predictions)
    if (length(unique(lengths)) != 1L) {
      stop("The replicated outcome must have the same length for every draw.", call. = FALSE)
    }

    out <- do.call(rbind, predictions)
    if (is.null(dim(out))) out <- matrix(out, nrow = length(predictions))
    rownames(out) <- paste0("draw", seq_len(nrow(out)))

    observed_value <- tryCatch(
      .rtmb_observed(fit, observed, if (is.null(spec)) list() else spec),
      error = function(e) NULL
    )
    if (!is.null(observed_value) && length(observed_value) == ncol(out)) {
      if (!is.null(names(observed_value))) colnames(out) <- names(observed_value)
      attr(out, "observed") <- observed_value
    }
    attr(out, "variable") <- prediction_name
    attr(out, "family") <- spec$family
    attr(out, "density") <- .rtmb_infer_density(fit$model)
    attr(out, "random") <- random
    class(out) <- c("rtmb_posterior_predict", "matrix", "array")
    out
  })
}

.rtmb_resolve_stat <- function(stat, env) {
  if (is.null(stat)) return(NULL)
  if (is.function(stat)) return(list(fun = stat, label = "statistic"))
  if (!is.character(stat) || length(stat) != 1L) {
    stop("'stat' must be NULL, a function, or the name of a function.", call. = FALSE)
  }
  fun <- tryCatch(
    get(stat, envir = env, mode = "function", inherits = TRUE),
    error = function(e) NULL
  )
  if (is.null(fun)) {
    fun <- tryCatch(match.fun(stat), error = function(e) NULL)
  }
  if (is.null(fun)) stop("Could not find statistic function '", stat, "'.", call. = FALSE)
  list(fun = fun, label = stat)
}

.rtmb_apply_stat <- function(x, stat) {
  x <- x[is.finite(x)]
  value <- stat$fun(x)
  if (!is.numeric(value) || length(value) != 1L) {
    stop("'stat' must return one numeric value.", call. = FALSE)
  }
  as.numeric(value)
}

.rtmb_resolve_pp_x <- function(fit, x, x_label, observed) {
  if (is.null(x)) return(NULL)

  if (is.character(x) && length(x) == 1L && identical(x, ".fitted")) {
    return(list(
      value = NULL,
      label = "Posterior predictive mean",
      name = ".fitted",
      fitted = TRUE
    ))
  }

  if (is.character(x) && length(x) == 1L && x %in% names(fit$model$data)) {
    x_name <- x
    x <- fit$model$data[[x_name]]
    x_label <- x_name
  } else {
    x_name <- NULL
  }

  if (is.data.frame(x) || is.matrix(x) || is.list(x) ||
      !(is.numeric(x) || is.integer(x) || is.factor(x) ||
        is.character(x) || is.logical(x))) {
    stop(
      "'x' must be a model-data column name or an atomic vector.",
      call. = FALSE
    )
  }
  if (length(x) != length(observed)) {
    if (is.character(x) && length(x) == 1L && is.null(x_name)) {
      stop("Predictor '", x, "' was not found in the model data.", call. = FALSE)
    }
    stop(
      "'x' has ", length(x), " values, but the observed outcome has ",
      length(observed), ".",
      call. = FALSE
    )
  }

  if (is.null(x_label) || !nzchar(x_label) || identical(x_label, "x")) {
    x_label <- if (is.null(x_name)) "Predictor" else x_name
  }

  list(value = x, label = x_label, name = x_name, fitted = FALSE)
}

.rtmb_scatter_positions <- function(x) {
  categorical <- is.factor(x) || is.character(x) || is.logical(x)
  if (categorical) {
    factor_x <- if (is.factor(x)) droplevels(x) else factor(x)
    base <- as.numeric(factor_x)
    ticks <- seq_along(levels(factor_x))
    labels <- levels(factor_x)
  } else {
    base <- as.numeric(x)
    ticks <- labels <- NULL
  }

  position <- base
  finite <- is.finite(base)
  unique_values <- sort(unique(base[finite]))
  spacing <- if (length(unique_values) > 1L) {
    min(diff(unique_values))
  } else {
    1
  }
  jitter_width <- spacing * if (categorical) 0.18 else 0.04
  duplicate_groups <- split(which(finite), base[finite], drop = TRUE)
  for (indices in duplicate_groups) {
    if (length(indices) > 1L) {
      position[indices] <- position[indices] +
        seq(-jitter_width, jitter_width, length.out = length(indices))
    }
  }

  list(
    position = position,
    categorical = categorical,
    ticks = ticks,
    labels = labels
  )
}

.rtmb_scatter_summary <- function(yrep, interval = 0.95) {
  alpha <- (1 - interval) / 2
  center <- colMeans(yrep, na.rm = TRUE)
  interval <- apply(
    yrep,
    2L,
    stats::quantile,
    probs = c(alpha, 1 - alpha),
    na.rm = TRUE,
    names = FALSE
  )
  if (is.null(dim(interval))) interval <- matrix(interval, nrow = 2L)
  list(center = center, lower = interval[1L, ], upper = interval[2L, ])
}

.rtmb_pp_check <- function(fit, type = c("auto", "dens", "bars"), stat = NULL,
                           x = NULL, x_label = NULL,
                           code_expr = NULL, code_env = parent.frame(),
                           observed = NULL, variable = NULL, draws = 100L,
                           seed = NULL,
                           random = c("conditional", "population", "simulate"),
                           chains = NULL, best_chains = NULL, plot = TRUE, ...) {
  type <- match.arg(type)
  random <- match.arg(random)
  spec <- fit$model$extra$posterior_predict
  observed_value <- .rtmb_observed(fit, observed, if (is.null(spec)) list() else spec)
  stat_info <- .rtmb_resolve_stat(stat, code_env)
  x_info <- .rtmb_resolve_pp_x(fit, x, x_label, observed_value)
  if (!is.null(x_info) && !is.null(stat_info)) {
    stop("'x' and 'stat' cannot be used together.", call. = FALSE)
  }

  yrep <- .rtmb_posterior_predict(
    fit = fit,
    code_expr = code_expr,
    code_env = code_env,
    variable = variable,
    draws = draws,
    seed = seed,
    random = random,
    chains = chains,
    best_chains = best_chains,
    observed = observed_value
  )

  if (ncol(yrep) != length(observed_value)) {
    stop(
      "The replicated outcome has ", ncol(yrep), " values, but 'observed' has ",
      length(observed_value), ".",
      call. = FALSE
    )
  }
  if (!is.null(x_info) && isTRUE(x_info$fitted)) {
    x_info$value <- colMeans(yrep, na.rm = TRUE)
  }

  density_type <- attr(yrep, "density")
  if (!is.null(x_info)) {
    type <- "scatter"
  } else if (identical(type, "auto")) {
    if (identical(density_type, "lpdf")) {
      type <- "dens"
    } else if (identical(density_type, "lpmf")) {
      type <- "bars"
    } else if (!is.null(stat_info)) {
      # The distribution display is not used when a statistic is requested.
      type <- "auto"
    } else {
      stop(
        "Could not determine whether the observation model uses an _lpdf or _lpmf function. Set 'type' to 'dens' or 'bars'.",
        call. = FALSE
      )
    }
  }

  out <- list(
    observed = observed_value,
    yrep = unclass(yrep),
    type = type,
    density = density_type,
    family = attr(yrep, "family"),
    random = random,
    stat = stat_info,
    predictor = if (is.null(x_info)) NULL else x_info$value,
    predictor_label = if (is.null(x_info)) NULL else x_info$label,
    predictor_name = if (is.null(x_info)) NULL else x_info$name,
    predictor_is_fitted = !is.null(x_info) && isTRUE(x_info$fitted)
  )

  if (!is.null(stat_info)) {
    out$observed_stat <- .rtmb_apply_stat(observed_value, stat_info)
    out$replicated_stat <- apply(yrep, 1L, .rtmb_apply_stat, stat = stat_info)
  }

  class(out) <- "rtmb_pp_check"
  if (isTRUE(plot)) graphics::plot(out, ...)
  invisible(out)
}

.rtmb_bar_summary <- function(observed, yrep, max_bins = 30L, interval = 0.95) {
  values <- c(observed, as.numeric(yrep))
  values <- values[is.finite(values)]
  if (length(values) == 0L) stop("No finite values are available for plotting.", call. = FALSE)

  value_range <- range(values)
  is_integer <- all(abs(values - round(values)) < 1e-8)
  if (is_integer && diff(value_range) <= max_bins - 1L) {
    centers <- seq.int(floor(value_range[1L]), ceiling(value_range[2L]))
    breaks <- c(centers - 0.5, centers[length(centers)] + 0.5)
    labels <- as.character(centers)
  } else {
    breaks <- pretty(value_range, n = max_bins)
    breaks <- sort(unique(c(min(values) - 1e-8, breaks, max(values) + 1e-8)))
    centers <- breaks[-length(breaks)] + diff(breaks) / 2
    labels <- format(centers, trim = TRUE, digits = 4L)
  }

  probabilities <- function(x) {
    graphics::hist(x[is.finite(x)], breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts /
      sum(is.finite(x))
  }

  obs_prob <- probabilities(observed)
  rep_prob <- t(apply(yrep, 1L, probabilities))
  alpha <- (1 - interval) / 2
  list(
    observed = obs_prob,
    median = apply(rep_prob, 2L, stats::median, na.rm = TRUE),
    lower = apply(rep_prob, 2L, stats::quantile, probs = alpha, na.rm = TRUE),
    upper = apply(rep_prob, 2L, stats::quantile, probs = 1 - alpha, na.rm = TRUE),
    labels = labels
  )
}

#' @param x An `rtmb_pp_check` object.
#' @param main Optional plot title.
#' @param xlab Optional x-axis label.
#' @param ylab Optional y-axis label.
#' @param observed_col Color for observed data.
#' @param predictive_col Color for predictive data and intervals.
#' @param show_legend Logical; display the plot legend.
#' @param legend_position Legend position. `"auto"` moves scatter-plot legends
#'   away from the fitted trend; the other values are passed to `legend()`.
#' @param legend_cex Relative text and symbol size for the legend.
#' @param interval Probability covered by posterior predictive intervals.
#' @rdname pp_check
#' @export
plot.rtmb_pp_check <- function(x, main = NULL, xlab = NULL, ylab = NULL,
                               observed_col = "#1B1B1B", predictive_col = "#2C7FB8",
                               show_legend = TRUE,
                               legend_position = c("auto", "topleft", "topright",
                                                   "bottomleft", "bottomright"),
                               legend_cex = 0.78,
                               interval = 0.95,
                               ...) {
  legend_position <- match.arg(legend_position)
  if (!is.numeric(legend_cex) || length(legend_cex) != 1L ||
      !is.finite(legend_cex) || legend_cex <= 0) {
    stop("'legend_cex' must be one positive number.", call. = FALSE)
  }
  if (!is.numeric(interval) || length(interval) != 1L ||
      !is.finite(interval) || interval <= 0 || interval >= 1) {
    stop("'interval' must be one number between 0 and 1.", call. = FALSE)
  }
  interval_label <- paste0(
    format(100 * interval, trim = TRUE, scientific = FALSE, digits = 4L),
    "%"
  )

  if (!is.null(x$predictor)) {
    positions <- .rtmb_scatter_positions(x$predictor)
    prediction <- .rtmb_scatter_summary(x$yrep, interval = interval)
    valid <- is.finite(positions$position) & is.finite(x$observed) &
      is.finite(prediction$center) & is.finite(prediction$lower) &
      is.finite(prediction$upper)
    if (!any(valid)) {
      stop("No complete predictor and outcome values are available for plotting.", call. = FALSE)
    }

    plot_x <- positions$position[valid]
    observed <- x$observed[valid]
    center <- prediction$center[valid]
    lower <- prediction$lower[valid]
    upper <- prediction$upper[valid]
    if (is.null(main)) main <- "Posterior predictive scatter"
    if (is.null(xlab)) xlab <- x$predictor_label
    if (is.null(ylab)) ylab <- "Outcome"
    y_range <- range(c(observed, lower, upper), finite = TRUE)
    if (diff(y_range) == 0) y_range <- y_range + c(-0.5, 0.5)
    x_range <- range(plot_x, finite = TRUE)
    if (diff(x_range) == 0) x_range <- x_range + c(-0.5, 0.5)

    graphics::plot(
      plot_x,
      observed,
      type = "n",
      xlim = x_range,
      ylim = y_range,
      xaxt = if (positions$categorical) "n" else "s",
      main = main,
      xlab = xlab,
      ylab = ylab,
      ...
    )
    if (positions$categorical) {
      graphics::axis(1L, at = positions$ticks, labels = positions$labels)
    }
    interval_col <- grDevices::adjustcolor(predictive_col, alpha.f = 0.28)
    identity_col <- grDevices::adjustcolor(observed_col, alpha.f = 0.5)
    if (isTRUE(x$predictor_is_fitted)) {
      graphics::abline(a = 0, b = 1, col = identity_col, lty = 2L, lwd = 1.2)
    }
    graphics::segments(plot_x, lower, plot_x, upper, col = interval_col, lwd = 1.2)
    graphics::points(plot_x, center, pch = 4L, col = predictive_col, lwd = 1.2)
    graphics::points(
      plot_x,
      observed,
      pch = 16L,
      col = grDevices::adjustcolor(observed_col, alpha.f = 0.65),
      cex = 0.75
    )
    legend_text <- c("Observed", "Predictive mean", paste(interval_label, "PI"))
    legend_pch <- c(16L, 4L, NA_integer_)
    legend_lty <- c(NA_integer_, NA_integer_, 1L)
    legend_col <- c(observed_col, predictive_col, interval_col)
    if (isTRUE(x$predictor_is_fitted)) {
      legend_text <- c(legend_text, "Identity")
      legend_pch <- c(legend_pch, NA_integer_)
      legend_lty <- c(legend_lty, 2L)
      legend_col <- c(legend_col, identity_col)
    }
    if (isTRUE(show_legend)) {
      legend_at <- legend_position
      if (identical(legend_at, "auto")) {
        trend <- suppressWarnings(stats::cor(plot_x, center, use = "complete.obs"))
        legend_at <- if (is.finite(trend) && trend < 0) "topright" else "topleft"
      }
      graphics::legend(
        legend_at,
        legend = legend_text,
        pch = legend_pch,
        lty = legend_lty,
        col = legend_col,
        lwd = rep(1.2, length(legend_text)),
        cex = legend_cex,
        x.intersp = 0.7,
        y.intersp = 0.85,
        bty = "n"
      )
    }
    return(invisible(x))
  }

  if (!is.null(x$stat)) {
    values <- x$replicated_stat[is.finite(x$replicated_stat)]
    if (length(values) == 0L) stop("The replicated statistic has no finite values.", call. = FALSE)
    if (is.null(main)) main <- "Posterior predictive statistic"
    if (is.null(xlab)) xlab <- x$stat$label
    if (is.null(ylab)) ylab <- "Frequency"
    graphics::hist(
      values,
      col = grDevices::adjustcolor(predictive_col, alpha.f = 0.55),
      border = "white",
      main = main,
      xlab = xlab,
      ylab = ylab,
      ...
    )
    graphics::abline(v = x$observed_stat, col = observed_col, lwd = 2L)
    if (isTRUE(show_legend)) {
      legend_at <- if (identical(legend_position, "auto")) "topright" else legend_position
      graphics::legend(
        legend_at, legend = "Observed", col = observed_col, lwd = 2L,
        cex = legend_cex, bty = "n"
      )
    }
    return(invisible(x))
  }

  if (identical(x$type, "dens")) {
    observed <- x$observed[is.finite(x$observed)]
    if (length(unique(observed)) < 2L) {
      stop("Density checks require at least two distinct observed values; use type = 'bars'.", call. = FALSE)
    }
    replicated_density <- lapply(seq_len(nrow(x$yrep)), function(i) {
      values <- x$yrep[i, ]
      values <- values[is.finite(values)]
      if (length(unique(values)) < 2L) NULL else stats::density(values)
    })
    replicated_density <- Filter(Negate(is.null), replicated_density)
    if (length(replicated_density) == 0L) {
      stop("Replicated outcomes do not contain enough distinct values for a density plot.", call. = FALSE)
    }

    observed_density <- stats::density(observed)
    x_range <- range(c(observed_density$x, unlist(lapply(replicated_density, `[[`, "x"))))
    y_max <- max(c(observed_density$y, unlist(lapply(replicated_density, `[[`, "y"))))
    if (is.null(main)) main <- "Posterior predictive density"
    if (is.null(xlab)) xlab <- "Outcome"
    if (is.null(ylab)) ylab <- "Density"
    graphics::plot(
      observed_density,
      type = "n",
      xlim = x_range,
      ylim = c(0, y_max * 1.05),
      main = main,
      xlab = xlab,
      ylab = ylab,
      ...
    )
    line_col <- grDevices::adjustcolor(predictive_col, alpha.f = 0.18)
    for (density_i in replicated_density) {
      graphics::lines(density_i, col = line_col, lwd = 1L)
    }
    graphics::lines(observed_density, col = observed_col, lwd = 2.5)
    if (isTRUE(show_legend)) {
      legend_at <- if (identical(legend_position, "auto")) "topright" else legend_position
      graphics::legend(
        legend_at, legend = c("Observed", "Replicated"),
        col = c(observed_col, predictive_col), lwd = c(2.5, 1.5),
        cex = legend_cex, bty = "n"
      )
    }
    return(invisible(x))
  }

  summary <- .rtmb_bar_summary(x$observed, x$yrep, interval = interval)
  labels <- summary$labels
  if (length(labels) > 15L) {
    keep <- unique(round(seq(1L, length(labels), length.out = 10L)))
    labels[-keep] <- ""
  }
  if (is.null(main)) main <- "Posterior predictive distribution"
  if (is.null(xlab)) xlab <- "Outcome"
  if (is.null(ylab)) ylab <- "Proportion"
  y_max <- max(c(summary$observed, summary$upper), na.rm = TRUE)
  mids <- graphics::barplot(
    summary$observed,
    names.arg = labels,
    col = grDevices::adjustcolor(observed_col, alpha.f = 0.65),
    border = NA,
    ylim = c(0, y_max * 1.12),
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
  )
  graphics::segments(mids, summary$lower, mids, summary$upper, col = predictive_col, lwd = 2L)
  graphics::points(mids, summary$median, pch = 19L, col = predictive_col)
  if (isTRUE(show_legend)) {
    legend_at <- if (identical(legend_position, "auto")) "topright" else legend_position
    graphics::legend(
      legend_at,
      legend = c("Observed", paste0("Replicated median (", interval_label, " interval)")),
      fill = c(grDevices::adjustcolor(observed_col, alpha.f = 0.65), NA),
      border = c(NA, NA), pch = c(NA, 19L), col = c(NA, predictive_col),
      cex = legend_cex, bty = "n"
    )
  }
  invisible(x)
}
