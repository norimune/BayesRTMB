#' @noRd
.print_code_impl <- function(self, private) {
  if (is.null(self$code)) {
    cat("Model code is not available (Not built with rtmb_code).\n")
    return(invisible(self))
  }

  cat("=== RTMB Model Code ===\n\n")
  cat("rtmb_code(\n")

  blocks <- setdiff(names(self$code), "env")
  is_empty_block <- function(expr) {
    is.call(expr) &&
      identical(expr[[1]], as.name("{")) &&
      length(as.list(expr)) == 1L
  }
  blocks <- blocks[!(blocks == "setup" & vapply(self$code[blocks], is_empty_block, logical(1)))]
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


