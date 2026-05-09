#' Print model code and related information
#'
#' @param model An \code{RTMB_Model} object.
#' @keywords internal
.print_code_impl <- function(model) {
  if (is.null(model$code)) {
    cat("Model code is not available (likely built from a simple objective function).\n")
    return(invisible(NULL))
  }
  
  cat("RTMB Model Code:\n")
  cat("--------------------------------------------------\n")
  # Data & Parameters headers
  cat("data {\n")
  for (n in names(model$data)) {
    d <- model$data[[n]]
    cat(sprintf("  %s (%s)\n", n, paste(dim(d) %||% length(d), collapse = "x")))
  }
  cat("}\n\n")
  
  cat("parameters {\n")
  for (n in names(model$par_list)) {
    p <- model$par_list[[n]]
    cat(sprintf("  %s (%s) [type: %s, random: %s]\n", 
                n, paste(p$dim %||% p$length, collapse = "x"), p$type, p$random))
  }
  cat("}\n\n")
  
  # Model block
  cat("model {\n")
  code_str <- deparse(model$code$model)
  # Remove the outer { } if present
  if (code_str[1] == "{") code_str <- code_str[-c(1, length(code_str))]
  cat(paste0("  ", code_str, collapse = "\n"))
  cat("\n}\n")
  cat("--------------------------------------------------\n")
}

#' @keywords internal
.print_log_prob_impl <- function(model) {
  if (is.null(model$log_prob)) {
    cat("Log-probability function is not defined.\n")
    return(invisible(NULL))
  }
  print(model$log_prob)
}

#' @keywords internal
.print_transform_impl <- function(model) {
  if (is.null(model$transform)) {
    cat("No transformation block defined.\n")
    return(invisible(NULL))
  }
  print(model$transform)
}

#' @keywords internal
.print_generate_impl <- function(model) {
  if (is.null(model$generate)) {
    cat("No generate block defined.\n")
    return(invisible(NULL))
  }
  print(model$generate)
}

# Helper for NULL dimensions
`%||%` <- function(a, b) if (!is.null(a)) a else b
