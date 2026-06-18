.rtmb_resolve_tape <- function(tape) {
  if (missing(tape) || is.null(tape)) return("auto")
  if (is.logical(tape)) return(if (isTRUE(tape)) "auto" else "none")
  match.arg(tape, c("auto", "none", "force"))
}

.rtmb_report_named_list <- function(x, block = "derived") {
  if (is.null(x)) return(invisible(NULL))
  if (!is.list(x)) {
    stop(sprintf("The '%s' block must return a named list.", block), call. = FALSE)
  }
  nms <- names(x)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop(sprintf("The '%s' block must return a named list.", block), call. = FALSE)
  }
  for (nm in nms) {
    assign(nm, x[[nm]], inherits = FALSE)
    eval(substitute(RTMB::REPORT(v), list(v = as.name(nm))))
  }
  invisible(NULL)
}

.rtmb_try_make_taped_derived <- function(data, par_template, derived_fn,
                                         transform_fn = NULL,
                                         extra_template = NULL,
                                         block = "derived",
                                         tape = c("auto", "none", "force"),
                                         progress = "none") {
  tape <- .rtmb_resolve_tape(tape)
  if (identical(tape, "none") || is.null(derived_fn)) return(NULL)

  if (is.null(extra_template)) extra_template <- list()
  base_names <- names(par_template)
  extra_names <- names(extra_template)
  parameters <- c(par_template, extra_template)

  if (length(parameters) == 0L || is.null(names(parameters)) || any(!nzchar(names(parameters)))) {
    if (identical(tape, "force")) {
      stop("Cannot tape derived quantities because the parameter template is empty or unnamed.", call. = FALSE)
    }
    return(NULL)
  }

  make_error <- NULL
  out <- tryCatch({
    f_ad <- function(theta) {
      para <- theta[base_names]
      if (!is.null(transform_fn)) {
        tran <- transform_fn(data, para)
        if (is.null(tran)) tran <- list()
        if (!is.list(tran)) {
          stop("'transform' must return a named list.", call. = FALSE)
        }
        para <- c(para, tran)
      }
      if (length(extra_names) > 0L) {
        para <- c(para, theta[extra_names])
      }
      res <- derived_fn(data, para)
      if (is.null(res)) res <- list()
      .rtmb_report_named_list(res, block = block)
      0
    }

    ad_obj <- RTMB::MakeADFun(
      func = f_ad,
      parameters = parameters,
      silent = TRUE
    )
    x0 <- as.numeric(unlist(parameters, use.names = FALSE))
    test_report <- ad_obj$report(x0)

    list(
      report = function(par, extra = NULL) {
        if (is.null(extra)) extra <- list()
        x <- as.numeric(unlist(c(par, extra), use.names = FALSE))
        ad_obj$report(x)
      },
      test = test_report
    )
  }, error = function(e) {
    make_error <<- e
    NULL
  })

  if (is.null(out) && identical(tape, "force")) {
    stop(
      sprintf("Failed to tape the '%s' block. Original error: %s", block, conditionMessage(make_error)),
      call. = FALSE
    )
  }
  out
}

.rtmb_tape_fallback_message <- function(block, err, progress = "none") {
  if (is.null(err)) return(invisible(NULL))
  .rtmb_progress_line(
    sprintf("Falling back to R evaluation for %s; tape evaluation failed: %s", block, conditionMessage(err)),
    progress
  )
  invisible(NULL)
}

.rtmb_skip_generate_error <- function(err, progress = "none") {
  warning(
    "Generated quantities were skipped because the 'generate' block failed: ",
    conditionMessage(err),
    call. = FALSE
  )
  .rtmb_progress_line("Generated quantities skipped.", progress)
  invisible(NULL)
}
