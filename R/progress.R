#' @noRd
.rtmb_resolve_progress <- function(progress = c("auto", "none", "bar", "message")) {
  progress <- match.arg(progress)
  if (isTRUE(getOption("BayesRTMB.silent", FALSE))) return("none")
  if (identical(progress, "none")) return("none")
  "message"
}

#' @noRd
.rtmb_progress_line <- function(msg, progress) {
  if (!identical(progress, "none") && nzchar(msg)) {
    cat(msg, "\n", sep = "")
    .rtmb_flush_console()
  }
}

#' @noRd
.rtmb_progress_start_line <- function(msg) {
  if (!isTRUE(getOption("BayesRTMB.silent", FALSE)) && nzchar(msg)) {
    cat(msg, "\n", sep = "")
    .rtmb_flush_console()
  }
}

#' @noRd
.rtmb_flush_console <- function() {
  if (interactive()) {
    try(utils::flush.console(), silent = TRUE)
  }
  invisible(NULL)
}

#' @noRd
.rtmb_progress_meter <- function(total, progress, label = "progress", message_step = 20L) {
  total <- max(1L, as.integer(total))
  current <- 0L
  mode <- match.arg(progress, c("none", "message"))
  next_percent <- as.integer(message_step)

  advance <- function(amount = 1L, msg = NULL, force = FALSE) {
    amount <- max(0L, as.integer(amount))
    current <<- min(total, current + amount)

    if (identical(mode, "none")) return(invisible(NULL))

    if (identical(mode, "message")) {
      if (!is.null(msg) && length(msg) > 0L && nzchar(msg[1])) {
        cat(msg[1], "\n", sep = "")
        .rtmb_flush_console()
      } else {
        pct <- floor(100 * current / total)
        while (pct >= next_percent && next_percent <= 100L) {
          cat(sprintf("%s: %d%%\n", label, next_percent))
          .rtmb_flush_console()
          next_percent <<- next_percent + as.integer(message_step)
        }
      }
      return(invisible(NULL))
    }
    invisible(NULL)
  }

  finish <- function() {
    if (identical(mode, "message") && next_percent <= 100L) {
      cat(sprintf("%s: 100%%\n", label))
      .rtmb_flush_console()
      next_percent <<- 110L
    }
    invisible(NULL)
  }

  list(advance = advance, finish = finish)
}

#' @noRd
.rtmb_progress_file_dir <- function() {
  stamp <- format(Sys.time(), "%Y%m%d%H%M%OS3")
  stamp <- gsub("[^0-9]", "", stamp)
  path <- file.path(tempdir(), paste0("BayesRTMB-progress-", Sys.getpid(), "-", stamp))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

#' @noRd
.rtmb_write_progress_file <- function(path, msg) {
  if (is.numeric(msg)) msg <- ""
  msg <- as.character(msg[1])
  if (!nzchar(msg)) return(invisible(FALSE))

  ok <- tryCatch({
    cat(msg, "\n", file = path, append = TRUE, sep = "")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  invisible(isTRUE(ok))
}

#' @noRd
.rtmb_read_progress_file <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  x <- x[nzchar(x)]
  if (length(x) == 0L) return(character(0))
  x
}

#' @noRd
.rtmb_report_progress_files <- function(files, line_counts) {
  for (i in seq_along(files)) {
    messages <- .rtmb_read_progress_file(files[i])
    n_messages <- length(messages)
    if (n_messages > line_counts[i]) {
      new_messages <- messages[(line_counts[i] + 1L):n_messages]
      cat(paste0(new_messages, "\n"), sep = "")
      .rtmb_flush_console()
      line_counts[i] <- n_messages
    }
  }
  line_counts
}

#' @noRd
.rtmb_collect_progress_futures <- function(futures, files, line_counts = NULL) {
  poll_interval <- getOption("BayesRTMB.progress_poll_interval", 0.2)
  if (!is.numeric(poll_interval) || length(poll_interval) != 1L ||
      !is.finite(poll_interval) || poll_interval <= 0) {
    poll_interval <- 0.2
  }

  if (is.null(line_counts) || length(line_counts) != length(files)) {
    line_counts <- integer(length(files))
  }
  repeat {
    line_counts <- .rtmb_report_progress_files(files, line_counts)
    if (all(vapply(futures, future::resolved, logical(1)))) break
    Sys.sleep(poll_interval)
  }
  line_counts <- .rtmb_report_progress_files(files, line_counts)
  lapply(futures, future::value)
}
