#' @noRd
.print_code_impl <- function(self, private) {
  if (is.null(self$code)) {
    cat("Model code is not available (Not built with rtmb_code).\n")
    return(invisible(self))
  }

  cat("=== RTMB Model Code ===\n\n")
  cat("rtmb_code(\n")

  blocks <- names(self$code)
  n_blocks <- length(blocks)

  for (i in seq_along(blocks)) {
    block <- blocks[i]
    expr <- self$code[[block]]

    lines <- deparse(expr, width.cutoff = 500L, control = "useSource")
    lines <- sapply(lines, function(x) {
      m <- regexpr("^ +", x)
      if (m > 0) {
        n_spaces <- attr(m, "match.length")
        new_spaces <- strrep(" ", n_spaces / 2)
        sub("^ +", new_spaces, x)
      } else {
        x
      }
    }, USE.NAMES = FALSE)

    lines <- gsub('^(\\s*)"#(.*)"$', "\\1#\\2", lines)
    lines <- paste0("  ", lines)

    if (trimws(lines[1]) == "{") {
      lines[1] <- paste0("  ", block, " = {")
    }

    if (i < n_blocks) {
      lines[length(lines)] <- paste0(lines[length(lines)], ", ")
    }

    cat(paste(lines, collapse = "\n"))
    cat("\n")
  }
  cat(")\n")
  return(invisible(self))
}

#' @noRd
.get_n_obs_impl <- function(self, private) {
  if (is.null(self$code$model)) return(NULL)

  get_base_name <- function(x) {
    if (is.name(x)) return(as.character(x))
    if (is.call(x) && identical(x[[1]], as.name("["))) return(get_base_name(x[[2]]))
    return(NULL)
  }

  data_names <- names(self$data)

  rewrite_to_count <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("~"))) {
        target <- x[[2]]
        base_name <- get_base_name(target)
        if (!is.null(base_name) && base_name %in% data_names) {
          return(call("<-", as.name("n_obs"),
                      call("+", as.name("n_obs"),
                           call("if", call("is.matrix", target),
                                call("nrow", target), call("length", target)))))
        } else {
          return(quote(NULL))
        }
      }
      x[] <- lapply(as.list(x), rewrite_to_count)
    }
    return(x)
  }

  raw_expr <- self$code$model
  count_expr <- rewrite_to_count(raw_expr)

  body_list <- c(
    list(as.name("{")),
    quote(RTMB::getAll(dat, par)),
    quote(n_obs <- 0),
    if (is.call(count_expr) && identical(count_expr[[1]], as.name("{"))) as.list(count_expr)[-1] else list(count_expr),
    quote(return(n_obs))
  )

  fn_expr <- call("function", as.pairlist(alist(dat = , par = )), as.call(body_list))
  count_fn <- eval(fn_expr, envir = asNamespace("BayesRTMB"))
  test_para <- self$get_par_list()

  n_total <- tryCatch({
    count_fn(self$data, test_para)
  }, error = function(e) return(NA_integer_))

  return(n_total)
}

#' @noRd
.evaluate_univariate_prior_impl <- function(self, private, p_name, u_val, par_list, caller_env) {
  ast <- self$code$model
  if (is.null(ast)) return(0)

  p_base <- gsub("\\[.*\\]$", "", p_name)
  p_info <- self$par_list[[p_base]]
  curr_par_list <- if (!is.null(par_list)) par_list else self$to_constrained(self$unconstrained_vector_to_list(self$init))

  p_idx <- 1
  if (grepl("\\[", p_name)) {
    idx_str <- gsub(".*\\[(.*)\\].*", "\\1", p_name)
    p_idx <- suppressWarnings(as.integer(idx_str))
    if (is.na(p_idx)) {
      p_names <- self$par_names[[p_base]]
      if (!is.null(p_names)) p_idx <- which(p_names == idx_str)
      if (length(p_idx) == 0) p_idx <- 1
    }
  }

  c_val <- u_val
  lower <- if (is.numeric(p_info$lower)) p_info$lower else -Inf
  upper <- if (is.numeric(p_info$upper)) p_info$upper else Inf

  if (is.infinite(lower) && is.infinite(upper)) {
    c_val <- u_val
  } else if (is.finite(lower) && is.infinite(upper)) {
    c_val <- exp(u_val) + lower
  } else if (is.finite(lower) && is.finite(upper)) {
    c_val <- lower + (upper - lower) / (1 + exp(-u_val))
  }
  curr_par_list[[p_base]][p_idx] <- c_val

  total_lp <- 0
  resolve_name <- function(expr, env) {
    if (is.name(expr)) return(as.character(expr))
    if (is.call(expr) && identical(expr[[1]], as.name("["))) {
      base <- as.character(expr[[2]])
      indices <- as.list(expr)[-(1:2)]
      idx_vals <- sapply(indices, function(idx) {
        val <- try(eval(idx, envir = env), silent = TRUE)
        if (inherits(val, "try-error")) return(as.character(idx))
        val
      })
      p_names <- self$par_names[[base]]
      if (!is.null(p_names) && length(idx_vals) == 1 && is.numeric(idx_vals[1])) {
        ii <- as.integer(idx_vals[1])
        if (ii >= 1 && ii <= length(p_names)) return(paste0(base, "[", p_names[ii], "]"))
      }
      return(paste0(base, "[", paste(idx_vals, collapse = ","), "]"))
    }
    return(deparse(expr))
  }

  eval_env <- new.env(parent = if (!is.null(self$code$env)) self$code$env else caller_env)
  eval_env$dat <- self$data
  eval_env$par <- curr_par_list
  eval_env$lp <- 0

  eval_env$`~` <- function(lhs, rhs) {
    lhs_expr <- substitute(lhs)
    rhs_expr <- substitute(rhs)
    res_name <- resolve_name(lhs_expr, env = parent.frame())
    dist_func <- if (is.call(rhs_expr)) deparse(rhs_expr[[1]]) else "none"

    is_multivariate <- dist_func %in% c("multi_normal", "multi_student_t", "multi_cauchy")

    if (is_multivariate) {
      base_match <- (res_name == p_base)
      if (base_match && p_name != p_base) {
        args <- as.list(rhs_expr)[-1]
        eval_args <- lapply(args, function(a) tryCatch(eval(a, envir = eval_env), error = function(e) e$message))
        if (dist_func == "multi_normal") {
          mu <- eval_args[[1]][p_idx]
          Sigma_kk <- eval_args[[2]][p_idx, p_idx]
          log_dens <- dnorm(u_val, mu, sqrt(Sigma_kk), log = TRUE)
        } else {
          idx_df <- if (dist_func == "multi_student_t") 1 else NULL
          idx_mu <- if (dist_func == "multi_student_t") 2 else 1
          idx_Sigma <- if (dist_func == "multi_student_t") 3 else 2
          df_val <- if (!is.null(idx_df)) eval_args[[idx_df]] else 1
          mu <- eval_args[[idx_mu]][p_idx]
          Sigma_kk <- eval_args[[idx_Sigma]][p_idx, p_idx]
          scale <- sqrt(as.numeric(Sigma_kk))
          val <- (u_val - mu) / scale
          log_dens <- dt(as.numeric(val), df = as.numeric(df_val), log = TRUE) - log(scale)
        }
        total_lp <<- total_lp + log_dens
      } else if (res_name == p_name) {
        full_call <- as.call(c(rhs_expr[[1]], lhs_expr, as.list(rhs_expr)[-1]))
        log_dens <- tryCatch(eval(full_call, envir = parent.frame()), warning = function(w) NULL, error = function(e) NULL)
        if (!is.null(log_dens)) total_lp <<- total_lp + sum(log_dens)
      }
    } else {
      if (is.call(rhs_expr)) {
        full_call <- as.call(c(rhs_expr[[1]], lhs_expr, as.list(rhs_expr)[-1]))
        log_dens <- tryCatch(eval(full_call, envir = parent.frame()), warning = function(w) NULL, error = function(e) NULL)
      } else {
        log_dens <- tryCatch(eval(rhs_expr, envir = parent.frame()), warning = function(w) NULL, error = function(e) NULL)
      }
      if (!is.null(log_dens)) {
        if (res_name == p_name) total_lp <<- total_lp + sum(log_dens)
        else if (res_name == p_base && p_name != p_base && length(log_dens) >= p_idx) total_lp <<- total_lp + log_dens[p_idx]
        else if (p_name == p_base && grepl(paste0("^", p_base, "\\["), res_name)) total_lp <<- total_lp + sum(log_dens)
      }
    }
    invisible(NULL)
  }

  eval_env$normal <- function(x, mean, sd) dnorm(x, mean, sd, log = TRUE)
  eval_env$cauchy <- function(x, location, scale) cauchy_lpdf(x, location, scale, sum = FALSE)
  eval_env$exponential <- function(x, rate) dexp(x, rate, log = TRUE)
  eval_env$gamma <- function(x, shape, rate) dgamma(x, shape, rate, log = TRUE)
  eval_env$inverse_gamma <- function(x, shape, scale) inverse_gamma_lpdf(x, shape, scale, sum = FALSE)
  eval_env$student_t <- function(x, df, mean, sd) student_t_lpdf(x, df, mu = mean, sigma = sd, sum = FALSE)
  eval_env$beta <- function(x, shape1, shape2) dbeta(x, shape1, shape2, log = TRUE)
  eval_env$bernoulli <- function(x, prob) dbinom(x, 1, prob, log = TRUE)
  eval_env$binomial <- function(x, size, prob) dbinom(x, size, prob, log = TRUE)
  eval_env$poisson <- function(x, lambda) dpois(x, lambda, log = TRUE)
  eval_env$multi_normal <- function(x, mean, Sigma) multi_normal_lpdf(x, mean, Sigma, sum = FALSE)
  eval_env$multi_student_t <- function(x, df, mean, Sigma) multi_student_t_lpdf(x, df, mean, Sigma, sum = FALSE)
  eval_env$multi_cauchy <- function(x, mean, Sigma) multi_cauchy_lpdf(x, mean, Sigma, sum = FALSE)

  eval_env$getAll <- function(d, p) {
    for (n in names(d)) assign(n, d[[n]], envir = parent.frame())
    for (n in names(p)) assign(n, p[[n]], envir = parent.frame())
  }
  eval_env$Dim <- function(...) NULL

  tryCatch({
    eval(as.call(c(list(as.name("{")), quote(getAll(dat, par)), ast)), envir = eval_env)
  }, error = function(e) NULL)

  lj <- 0
  if (is.finite(lower) && is.infinite(upper)) lj <- u_val
  else if (is.finite(lower) && is.finite(upper)) lj <- log(upper - lower) - u_val - 2 * log(1 + exp(-u_val))

  return(total_lp + lj)
}

#' @noRd
remove_prior_from_ast <- function(expr, target_spec, p_names = NULL) {
  # target_spec is like "b" or "b[talk]"
  target_base <- gsub("\\[.*\\]$", "", target_spec)
  has_index <- grepl("\\[", target_spec)
  target_idx_str <- if (has_index) gsub(".*\\[(.*)\\].*", "\\1", target_spec) else NULL

  if (is.call(expr)) {
    # 1. Handle { ... } blocks
    if (identical(expr[[1]], as.name("{"))) {
      new_args <- lapply(as.list(expr)[-1], function(e) remove_prior_from_ast(e, target_spec, p_names))
      new_args <- Filter(Negate(is.null), new_args)
      return(as.call(c(as.name("{"), new_args)))
    }

    # 2. Handle sampling statements b[k] ~ ...
    if (identical(expr[[1]], as.name("~"))) {
      lhs <- expr[[2]]
      lhs_str <- deparse(lhs)

      # If target is "b" (full vector), remove any prior for "b" or "b[...]"
      if (!has_index) {
        if (identical(lhs, as.name(target_base)) ||
          (is.call(lhs) && identical(lhs[[1]], as.name("[")) && identical(lhs[[2]], as.name(target_base)))) {
          return(NULL)
        }
      } else {
        # Target is specific element "b[talk]"
        # If LHS is exactly "b[talk]", remove it.
        if (lhs_str == target_spec) {
          return(NULL)
        }

        # --- Advanced Slicing for Vectorized Priors ---
        if (identical(lhs, as.name(target_base))) {
          rhs <- expr[[3]]
          idx_val <- target_idx_str
          if (!is.null(p_names)) {
              idx_match <- which(p_names == target_idx_str)
              if (length(idx_match) > 0) idx_val <- idx_match
          }
          neg_idx <- bquote(-.(as.numeric(idx_val)))
          new_lhs <- bquote(.(as.name(target_base))[.(neg_idx)])
          new_rhs <- rhs
          for (i in 2:length(rhs)) {
              arg <- rhs[[i]]
              if (is.call(arg) || is.symbol(arg) || (is.numeric(arg) && length(arg) > 1)) {
                  new_rhs[[i]] <- bquote(.(arg)[.(neg_idx)])
              }
          }
          return(as.call(list(as.name("~"), new_lhs, new_rhs)))
        }
      }
    }

    # 3. Handle for loops
    if (identical(expr[[1]], as.name("for"))) {
      var_name <- as.character(expr[[2]])
      body <- expr[[4]]

      if (has_index && is.call(body) && identical(body[[1]], as.name("~"))) {
        if (is.call(body[[2]]) && identical(body[[2]][[2]], as.name(target_base))) {
          idx_val <- target_idx_str
          if (!is.null(p_names)) {
            idx_match <- which(p_names == target_idx_str)
            if (length(idx_match) > 0) idx_val <- idx_match
          }
          new_body <- bquote(if (.(as.name(var_name)) != .(as.numeric(idx_val))) .(body))
          expr[[4]] <- new_body
          return(expr)
        }
      }

      new_body <- remove_prior_from_ast(body, target_spec, p_names)
      if (is.null(new_body)) return(NULL)
      expr[[4]] <- new_body
      return(expr)
    }

    new_args <- lapply(as.list(expr), function(e) remove_prior_from_ast(e, target_spec, p_names))
    return(as.call(new_args))
  }
  return(expr)
}
