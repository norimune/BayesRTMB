#' Create a null model by fixing specified parameters
#'
#' @param model An \code{RTMB_Model} object.
#' @param target Character string specifying the target parameter and its prior.
#' @param value Numeric value to fix parameters to. Default is 0.
#' @return A new \code{RTMB_Model} object with the specified parameters fixed.
#' @keywords internal
.null_model_impl <- function(model, target, value = 0) {
  # (RTMB_Model.R の L2782-3101 のロジックを抽出)
  # self を model に置き換えて実装
  
  if (!is.character(target) || length(target) != 1) {
    stop("target must be specified as a single string (e.g., 'delta ~ cauchy(0, r)' or 'delta').", call. = FALSE)
  }

  # Auto-completion processing...
  if (!grepl("~", target)) {
    target_ast <- str2lang(target)
    has_index <- is.call(target_ast) && identical(target_ast[[1]], as.name("["))
    if (has_index) {
      base_name <- as.character(target_ast[[2]])
      idx_expr <- target_ast[[3]]
    } else {
      base_name <- as.character(target_ast)
      idx_expr <- NULL
    }

    find_prior <- function(expr) {
      if (is.call(expr)) {
        if (identical(expr[[1]], as.name("~"))) {
          lhs <- expr[[2]]
          if (identical(lhs, as.name(base_name)) ||
            (is.call(lhs) && identical(lhs[[1]], as.name("[")) && identical(lhs[[2]], as.name(base_name)))) {
            return(expr[[3]])
          }
        }
        for (i in seq_along(expr)) {
          res <- find_prior(expr[[i]])
          if (!is.null(res)) return(res)
        }
      }
      return(NULL)
    }

    prior_expr <- find_prior(model$code$model)
    if (is.null(prior_expr)) {
      stop(sprintf("The prior distribution definition for '%s' was not found in the model code.", base_name), call. = FALSE)
    }

    if (has_index) {
      eval_env <- list2env(model$data, parent = if (!is.null(model$code$env)) model$code$env else parent.frame())
      for (n in names(model$init)) assign(n, model$init[[n]], envir = eval_env)

      idx_val <- tryCatch({
        res <- eval(idx_expr, envir = eval_env)
        if (length(res) == 1 && (is.numeric(res) || is.character(res))) res else as.character(idx_expr)
      }, error = function(e) as.character(idx_expr))

      for (i in 2:length(prior_expr)) {
        arg_expr <- prior_expr[[i]]
        if (!is.numeric(arg_expr)) {
          arg_val <- tryCatch(eval(arg_expr, envir = eval_env), warning = function(w) NULL, error = function(e) NULL)
          if (!is.null(arg_val) && length(arg_val) > 1) {
            ext_val <- arg_val[idx_val]
            if (length(ext_val) == 1 && !is.na(ext_val)) prior_expr[[i]] <- ext_val else prior_expr[[i]] <- call("[", arg_expr, idx_expr)
          }
        }
      }
    }
    new_formula <- call("~", target_ast, prior_expr)
    target <- paste(deparse(new_formula), collapse = " ")
    cat(sprintf("Auto-completed target: %s\n", target))
  } else {
    target_f <- str2lang(target)
    target_ast <- target_f[[2]]
    lhs <- target_ast
    has_index <- is.call(lhs) && identical(lhs[[1]], as.name("["))
    if (has_index) {
      base_name <- as.character(lhs[[2]])
      idx_expr <- lhs[[3]]
      eval_env <- list2env(model$data, parent = if (!is.null(model$code$env)) model$code$env else parent.frame())
      for (n in names(model$init)) assign(n, model$init[[n]], envir = eval_env)
      idx_val <- tryCatch({
        res <- eval(idx_expr, envir = eval_env)
        if (length(res) == 1 && (is.numeric(res) || is.character(res))) res else as.character(idx_expr)
      }, error = function(e) as.character(idx_expr))
    } else {
      base_name <- as.character(lhs)
      idx_expr <- idx_val <- NULL
    }
  }

  f <- tryCatch(as.formula(target), error = function(e) stop("Failed to parse 'target'.", call. = FALSE))
  spec <- deparse(f[[2]])
  prior_expr <- f[[3]]
  prior_str <- deparse(prior_expr)

  structural_bounds <- c("ordered", "positive_ordered", "simplex", "corr_matrix", "CF_corr", "cov_matrix", "CF_cov", "centered", "centered_matrix", "centered_tri", "positive_centered_tri", "lower_tri_stz")
  no_fix_bounds <- c("simplex", "corr_matrix", "cov_matrix", "CF_cov", "centered", "centered_matrix", "centered_tri", "positive_centered_tri")
  elementwise_bounds <- c("none", "lower", "upper", "interval")

  all_flat_names <- character(0)
  flat_to_info <- list()
  for (name in names(model$par_list)) {
    p <- model$par_list[[name]]
    fnames <- generate_flat_names(name, p$dim, model$par_names[[name]])
    for (i in seq_along(fnames)) flat_to_info[[fnames[i]]] <- list(par = name, idx = i)
    all_flat_names <- c(all_flat_names, fnames)
  }

  targets <- list()
  if (spec %in% names(model$par_list)) {
    p <- model$par_list[[spec]]
    if (p$bounds %in% no_fix_bounds) stop(sprintf("Parameter '%s' cannot be fixed.", spec), call. = FALSE)
    targets[[spec]] <- 1:p$length
  } else if (spec %in% all_flat_names) {
    info <- flat_to_info[[spec]]
    p <- model$par_list[[info$par]]
    if (p$bounds %in% structural_bounds) stop(sprintf("Elements of '%s' cannot be individually fixed.", spec), call. = FALSE)
    targets[[info$par]] <- unique(c(targets[[info$par]], info$idx))
  } else {
    stop(sprintf("'%s' is not a valid parameter.", spec), call. = FALSE)
  }

  par_name <- names(targets)[1]
  fix_indices <- targets[[par_name]]
  k_fixed <- length(fix_indices)

  if (is.call(prior_expr)) {
    fn_name <- as.character(prior_expr[[1]])
    prior_expr[[1]] <- as.name(paste0(fn_name, "_lpdf"))
    prior_expr <- as.call(append(as.list(prior_expr), call("rep", value, k_fixed), after = 1))
    eval_env <- list2env(model$data, parent = parent.frame())
    init_list <- constrained_vector_to_list(model$init, model$par_list)
    for (n in names(init_list)) assign(n, init_list[[n]], envir = eval_env)
    if (has_index) {
      assign("k", idx_val, envir = eval_env)
      if (is.symbol(idx_expr)) assign(as.character(idx_expr), idx_val, envir = eval_env)
    }
    correction_val <- tryCatch(sum(eval(prior_expr, envir = eval_env)), error = function(e) stop(sprintf("Failed to evaluate prior: %s", e$message)))
  } else {
    stop("Priors must be calls.")
  }

  p <- model$par_list[[par_name]]
  if (!is.null(p$lower) && any(value <= p$lower)) stop("Value must be > lower bound.", call. = FALSE)
  if (!is.null(p$upper) && any(value >= p$upper)) stop("Value must be < upper bound.", call. = FALSE)

  b_type <- p$bounds
  lj_val <- 0
  if (b_type %in% c("lower", "upper", "interval")) {
    val_vec <- rep(value, length.out = k_fixed)
    if (b_type == "lower") lj_val <- sum(log(val_vec - p$lower))
    else if (b_type == "upper") lj_val <- sum(log(p$upper - val_vec))
    else if (b_type == "interval") {
      prob <- (val_vec - p$lower) / (p$upper - p$lower)
      u_val <- log(prob / (1 - prob))
      lj_val <- sum(log(p$upper - p$lower) - u_val - 2 * log(1 + exp(-u_val)))
    }
    correction_val <- correction_val + lj_val
  }

  new_model <- model$clone()
  prior_removed <- FALSE
  if (!is.null(model$code)) {
    new_model$code <- as.list(model$code)
    old_ast <- paste(deparse(new_model$code$model), collapse = "\n")
    new_model$code$model <- remove_prior_from_ast(model$code$model, spec, model$par_names[[base_name]])
    if (paste(deparse(new_model$code$model), collapse = "\n") != old_ast) {
      prior_removed <- TRUE
    } else {
      cancel_call <- prior_expr
      cancel_call[[2]] <- target_ast
      sub_expr <- bquote(lp <- lp - .(cancel_call))
      ast_list <- as.list(new_model$code$model)
      ret_idx <- which(sapply(ast_list, function(x) is.call(x) && identical(x[[1]], as.name("return"))))
      if (length(ret_idx) > 0) {
        new_model$code$model <- as.call(append(ast_list, list(sub_expr), after = ret_idx[1] - 1))
        prior_removed <- TRUE
      }
    }
    user_env <- if (!is.null(model$code$env)) model$code$env else parent.frame()
    new_model$log_prob <- eval(substitute(BayesRTMB::model_code(E, env = EN), list(E = new_model$code$model, EN = user_env)))
  }

  is_full <- length(fix_indices) == p$length
  map_list <- if (is.null(model$map)) list() else model$map
  if (is_full) {
    map_list[[par_name]] <- factor(rep(NA, p$unc_length))
  } else {
    if (!(par_name %in% names(map_list))) map_list[[par_name]] <- factor(1:p$unc_length)
    current <- as.integer(map_list[[par_name]])
    current[fix_indices] <- NA
    non_na <- !is.na(current)
    if (any(non_na)) current[non_na] <- seq_len(sum(non_na))
    map_list[[par_name]] <- factor(ifelse(is.na(current), NA, paste0(par_name, "_", current)))
  }

  adjusted_init <- if (is.list(model$init)) model$init else constrained_vector_to_list(model$init, model$par_list)
  val <- as.numeric(adjusted_init[[par_name]])
  val[fix_indices] <- value
  adjusted_init[[par_name]] <- val
  if (length(p$dim) > 1) dim(adjusted_init[[par_name]]) <- p$dim

  new_model$map <- map_list
  new_model$init <- adjusted_init
  new_model$prior_correction <- if (prior_removed) lj_val else correction_val + lj_val

  cat("Null model created.\n")
  return(new_model)
}

#' Evaluate the univariate prior density contribution
#'
#' @param model An \code{RTMB_Model} object.
#' @param p_name Flattened parameter name (e.g., "beta[1]").
#' @param u_val Value in unconstrained space.
#' @param par_list Optional parameter list (constrained).
#' @return Numeric log-prior density.
#' @keywords internal
.evaluate_univariate_prior_impl <- function(model, p_name, u_val, par_list = NULL) {
  ast <- model$code$model
  if (is.null(ast)) return(0)

  p_base <- gsub("\\[.*\\]$", "", p_name)
  p_info <- model$par_list[[p_base]]
  curr_par_list <- if (!is.null(par_list)) par_list else to_constrained(unconstrained_vector_to_list(model$init, model$par_list), model$par_list)

  p_idx <- 1
  if (grepl("\\[", p_name)) {
    idx_str <- gsub(".*\\[(.*)\\].*", "\\1", p_name)
    p_idx <- suppressWarnings(as.integer(idx_str))
    if (is.na(p_idx)) {
      p_names <- model$par_names[[p_base]]
      if (!is.null(p_names)) p_idx <- which(p_names == idx_str)
      if (length(p_idx) == 0) p_idx <- 1
    }
  }

  lower <- if (is.numeric(p_info$lower)) p_info$lower else -Inf
  upper <- if (is.numeric(p_info$upper)) p_info$upper else Inf
  c_val <- if (is.infinite(lower) && is.infinite(upper)) u_val else if (is.finite(lower) && is.infinite(upper)) exp(u_val) + lower else lower + (upper - lower) / (1 + exp(-u_val))
  curr_par_list[[p_base]][p_idx] <- c_val

  total_lp <- 0
  resolve_name <- function(expr, env) {
    if (is.name(expr)) return(as.character(expr))
    if (is.call(expr) && identical(expr[[1]], as.name("["))) {
      base <- as.character(expr[[2]])
      idx_vals <- sapply(as.list(expr)[-(1:2)], function(idx) tryCatch(eval(idx, envir = env), error = function(e) as.character(idx)))
      p_names <- model$par_names[[base]]
      if (!is.null(p_names) && length(idx_vals) == 1 && is.numeric(idx_vals[1])) {
        ii <- as.integer(idx_vals[1])
        if (ii >= 1 && ii <= length(p_names)) return(paste0(base, "[", p_names[ii], "]"))
      }
      return(paste0(base, "[", paste(idx_vals, collapse = ","), "]"))
    }
    return(deparse(expr))
  }

  eval_env <- new.env(parent = if (!is.null(model$code$env)) model$code$env else parent.frame())
  eval_env$dat <- model$data
  eval_env$par <- curr_par_list
  eval_env$lp <- 0

  eval_env$`~` <- function(lhs, rhs) {
    lhs_expr <- substitute(lhs)
    rhs_expr <- substitute(rhs)
    res_name <- resolve_name(lhs_expr, env = parent.frame())
    dist_func <- if (is.call(rhs_expr)) deparse(rhs_expr[[1]]) else "none"

    if (dist_func %in% c("multi_normal", "multi_student_t", "multi_cauchy")) {
      if (res_name == p_base && p_name != p_base) {
        args <- as.list(rhs_expr)[-1]
        eval_args <- lapply(args, function(a) eval(a, envir = parent.frame()))
        if (dist_func == "multi_normal") {
          total_lp <<- total_lp + dnorm(u_val, eval_args[[1]][p_idx], sqrt(eval_args[[2]][p_idx, p_idx]), log = TRUE)
        } else {
          idx_df <- if (dist_func == "multi_student_t") 1 else NULL
          idx_mu <- if (dist_func == "multi_student_t") 2 else 1
          idx_Sigma <- if (dist_func == "multi_student_t") 3 else 2
          df_val <- if (!is.null(idx_df)) eval_args[[idx_df]] else 1
          mu <- eval_args[[idx_mu]][p_idx]
          Sigma_kk <- eval_args[[idx_Sigma]][p_idx, p_idx]
          scale <- sqrt(as.numeric(Sigma_kk))
          total_lp <<- total_lp + dt((u_val - mu) / scale, df = as.numeric(df_val), log = TRUE) - log(scale)
        }
      } else if (res_name == p_name) {
        full_call <- as.call(c(rhs_expr[[1]], lhs_expr, as.list(rhs_expr)[-1]))
        total_lp <<- total_lp + sum(eval(full_call, envir = parent.frame()))
      }
    } else {
      full_call <- if (is.call(rhs_expr)) as.call(c(rhs_expr[[1]], lhs_expr, as.list(rhs_expr)[-1])) else rhs_expr
      log_dens <- tryCatch(eval(full_call, envir = parent.frame()), error = function(e) NULL)
      if (!is.null(log_dens)) {
        if (res_name == p_name) total_lp <<- total_lp + sum(log_dens)
        else if (res_name == p_base && p_name != p_base && length(log_dens) >= p_idx) total_lp <<- total_lp + log_dens[p_idx]
        else if (p_name == p_base && grepl(paste0("^", p_base, "\\["), res_name)) total_lp <<- total_lp + sum(log_dens)
      }
    }
  }

  eval_env$normal <- function(x, mean, sd) dnorm(x, mean, sd, log = TRUE)
  eval_env$cauchy <- function(x, location, scale) BayesRTMB::cauchy_lpdf(x, location, scale, sum = FALSE)
  eval_env$exponential <- function(x, rate) dexp(x, rate, log = TRUE)
  eval_env$gamma <- function(x, shape, rate) dgamma(x, shape, rate, log = TRUE)
  eval_env$inverse_gamma <- function(x, shape, scale) BayesRTMB::dinverse_gamma(x, shape, scale, log = TRUE)
  eval_env$student_t <- function(x, df, mean, sd) BayesRTMB::dstudent_t(x, df, mean, sd, log = TRUE)
  eval_env$beta <- function(x, shape1, shape2) dbeta(x, shape1, shape2, log = TRUE)
  eval_env$bernoulli <- function(x, prob) dbinom(x, 1, prob, log = TRUE)
  eval_env$binomial <- function(x, size, prob) dbinom(x, size, prob, log = TRUE)
  eval_env$poisson <- function(x, lambda) dpois(x, lambda, log = TRUE)
  eval_env$multi_normal <- function(x, mean, Sigma) BayesRTMB::multi_normal_lpdf(x, mean, Sigma, sum = FALSE)
  eval_env$multi_student_t <- function(x, df, mean, Sigma) BayesRTMB::multi_student_t_lpdf(x, df, mean, Sigma, sum = FALSE)
  eval_env$multi_cauchy <- function(x, mean, Sigma) BayesRTMB::multi_cauchy_lpdf(x, mean, Sigma, sum = FALSE)
  eval_env$getAll <- function(d, p) { for(n in names(d)) assign(n, d[[n]], envir=parent.frame()); for(n in names(p)) assign(n, p[[n]], envir=parent.frame()) }
  eval_env$Dim <- function(...) NULL

  tryCatch(eval(as.call(c(list(as.name("{")), quote(getAll(dat, par)), ast)), envir = eval_env), error = function(e) NULL)

  lj <- 0
  if (is.finite(lower) && is.infinite(upper)) lj <- u_val
  else if (is.finite(lower) && is.finite(upper)) lj <- log(upper - lower) - u_val - 2 * log(1 + exp(-u_val))
  return(total_lp + lj)
}

#' Internal function to remove a prior from AST
#' @keywords internal
remove_prior_from_ast <- function(expr, target_spec, p_names = NULL) {
  target_base <- gsub("\\[.*\\]$", "", target_spec)
  has_index <- grepl("\\[", target_spec)

  if (is.call(expr)) {
    if (identical(expr[[1]], as.name("{"))) {
      new_args <- lapply(as.list(expr)[-1], function(e) remove_prior_from_ast(e, target_spec, p_names))
      return(as.call(c(as.name("{"), Filter(Negate(is.null), new_args))))
    }
    if (identical(expr[[1]], as.name("~"))) {
      lhs <- expr[[2]]
      lhs_str <- deparse(lhs)
      if (!has_index) {
        if (identical(lhs, as.name(target_base)) || (is.call(lhs) && identical(lhs[[1]], as.name("[")) && identical(lhs[[2]], as.name(target_base)))) return(NULL)
      } else {
        if (lhs_str == target_spec) return(NULL)
      }
    }
    if (identical(expr[[1]], as.name("for"))) {
      var_name <- as.character(expr[[2]])
      body <- expr[[4]]
      if (has_index && is.call(body) && identical(body[[1]], as.name("~"))) {
        if (is.call(body[[2]]) && identical(body[[2]][[2]], as.name(target_base))) {
          target_idx_str <- gsub(".*\\[(.*)\\].*", "\\1", target_spec)
          if (!is.null(p_names)) {
            idx_match <- which(p_names == target_idx_str)
            if (length(idx_match) > 0) target_idx_str <- idx_match
          }
          expr[[4]] <- bquote(if (.(as.name(var_name)) != .(target_idx_str)) .(body))
          return(expr)
        }
      }
      new_body <- remove_prior_from_ast(body, target_spec, p_names)
      if (is.null(new_body)) return(NULL)
      expr[[4]] <- new_body
      return(expr)
    }
    return(as.call(lapply(as.list(expr), function(e) remove_prior_from_ast(e, target_spec, p_names))))
  }
  return(expr)
}
