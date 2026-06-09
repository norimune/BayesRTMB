#' @noRd
.rtmb_resolve_progress <- function(progress = c("auto", "none", "bar", "message")) {
  progress <- match.arg(progress)
  if (isTRUE(getOption("BayesRTMB.silent", FALSE))) return("none")
  if (identical(progress, "auto")) {
    return(if (interactive()) "bar" else "message")
  }
  progress
}

#' @noRd
.rtmb_progress_line <- function(msg, progress) {
  if (!identical(progress, "none") && nzchar(msg)) cat(msg, "\n", sep = "")
}

#' @noRd
.rtmb_progress_meter <- function(total, progress, label = "progress", message_step = 10L) {
  total <- max(1L, as.integer(total))
  current <- 0L
  mode <- match.arg(progress, c("none", "bar", "message"))
  update_every <- max(1L, floor(total / 100L))
  next_bar <- update_every
  next_percent <- as.integer(message_step)
  pb <- NULL

  if (identical(mode, "bar")) {
    pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  }

  advance <- function(amount = 1L, msg = NULL, force = FALSE) {
    amount <- max(0L, as.integer(amount))
    current <<- min(total, current + amount)

    if (identical(mode, "none")) return(invisible(NULL))

    if (identical(mode, "message")) {
      if (!is.null(msg) && length(msg) > 0L && nzchar(msg[1])) {
        cat(msg[1], "\n", sep = "")
      } else {
        pct <- floor(100 * current / total)
        while (pct >= next_percent && next_percent <= 100L) {
          cat(sprintf("%s: %d%%\n", label, next_percent))
          next_percent <<- next_percent + as.integer(message_step)
        }
      }
      return(invisible(NULL))
    }

    if ((current >= next_bar || isTRUE(force)) && !is.null(pb)) {
      utils::setTxtProgressBar(pb, current)
      while (next_bar <= current) next_bar <<- next_bar + update_every
    }
    invisible(NULL)
  }

  finish <- function() {
    if (identical(mode, "bar") && !is.null(pb)) {
      utils::setTxtProgressBar(pb, total)
      close(pb)
      pb <<- NULL
    } else if (identical(mode, "message") && next_percent <= 100L) {
      cat(sprintf("%s: 100%%\n", label))
      next_percent <<- 110L
    }
    invisible(NULL)
  }

  list(advance = advance, finish = finish)
}

#' @noRd
.rtmb_progressr_message_handler <- function() {
  progressr::make_progression_handler(
    name = "bayesrtmb_message",
    reporter = list(
      update = function(config, state, progression, ...) {
        msg <- state$message
        if (!is.null(msg) && length(msg) > 0L && nzchar(msg[1])) {
          cat(msg[1], "\n", sep = "")
        }
      }
    ),
    clear = FALSE,
    intrusiveness = 0,
    target = "terminal"
  )
}

#' @noRd
.rtmb_progressr_handler <- function(progress) {
  progress <- match.arg(progress, c("none", "bar", "message"))
  if (identical(progress, "none")) return(NULL)
  if (identical(progress, "message")) return(.rtmb_progressr_message_handler())
  progressr::handler_txtprogressbar(style = 3L, clear = FALSE)
}
