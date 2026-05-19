#' Convert Wide Data to Long Format
#'
#' @description
#' A highly intuitive wrapper around \code{stats::reshape} designed for 
#' psychological research. It converts data from wide format to long format
#' by identifying within-subjects factors.
#'
#' @param data A data frame in wide format.
#' @param within The columns to gather into long format. Can be:
#'   \itemize{
#'     \item A character vector of column names.
#'     \item A range string like \code{"time1:time4"}.
#'     \item A prefix string like \code{"time"} (matches all columns starting with "time").
#'   }
#' @param label The name for the new column that will contain the level names (e.g., "Time"). 
#'   If NULL, it defaults to the prefix used in \code{within} or "Condition".
#' @param value The name for the new column that will contain the measurement values. 
#'   Default is "Value".
#' @param id The identifier columns that should be repeated for each row. 
#'   If NULL (default), all columns NOT specified in \code{within} are treated as IDs.
#' @return A data frame in long format.
#' @export
to_long <- function(data, within = NULL, label = NULL, value = "Value", id = NULL) {
  if (is.null(within)) {
    stop("Argument 'within' is required. Specify columns to gather (e.g., 'time1:time4' or 'time').")
  }
  
  all_names <- names(data)
  target_cols <- NULL
  inferred_label <- "Condition"

  # 1. Handle 'within' specification
  if (is.character(within) && length(within) == 1) {
    if (grepl(":", within)) {
      # Range notation: "time1:time4"
      parts <- strsplit(within, ":")[[1]]
      start_idx <- which(all_names == parts[1])
      end_idx <- which(all_names == parts[2])
      if (length(start_idx) == 1 && length(end_idx) == 1) {
        target_cols <- all_names[start_idx:end_idx]
      } else {
        stop(sprintf("Range columns '%s' or '%s' not found.", parts[1], parts[2]))
      }
    } else if (within %in% all_names) {
      # Single column name
      target_cols <- within
      inferred_label <- within
    } else {
      # Prefix matching: "time" -> time1, time2, etc.
      target_cols <- all_names[grepl(paste0("^", within), all_names)]
      if (length(target_cols) > 0) {
        inferred_label <- within
      } else {
        stop(sprintf("No columns found starting with prefix '%s'.", within))
      }
    }
  } else {
    # Vector of names or indices
    if (is.numeric(within)) {
      target_cols <- all_names[within]
    } else {
      target_cols <- within
    }
  }
  
  if (is.null(target_cols) || length(target_cols) == 0) {
    stop("Could not identify columns to gather based on 'within' specification.")
  }

  # 2. Determine label column name
  if (is.null(label)) {
    label <- inferred_label
  }

  # 3. Determine ID columns (everything not in target_cols)
  if (is.null(id)) {
    id <- setdiff(all_names, target_cols)
  }

  # 4. Perform reshape
  data_res <- as.data.frame(data)
  
  # stats::reshape requires a list of varying column names
  res <- stats::reshape(data_res, 
                        direction = "long", 
                        varying = list(target_cols), 
                        v.names = value, 
                        timevar = label, 
                        times = target_cols, 
                        idvar = id)
  
  # Sort by ID and label for readability
  res <- res[do.call(order, res[c(id, label)]), ]
  rownames(res) <- NULL
  
  # Ensure the label column is a factor with levels in original column order
  res[[label]] <- factor(res[[label]], levels = target_cols)
  
  return(res)
}

#' Convert Long Data to Wide Format
#'
#' @description
#' A user-friendly wrapper around \code{stats::reshape} to convert data from
#' long format back to wide format.
#'
#' @param data A data frame in long format.
#' @param within The name of the column containing within-subjects factor levels (e.g., "Condition").
#' @param value The name of the column containing measurement values (e.g., "Value").
#' @param id The identifier columns (e.g., "Subject", "Group"). If NULL, 
#'   inferred as all columns except \code{within} and \code{value}.
#' @return A data frame in wide format.
#' @export
to_wide <- function(data, within = "Condition", value = "Value", id = NULL) {
  df <- as.data.frame(data)
  if (is.null(id)) {
    id <- setdiff(names(df), c(within, value))
  }

  res <- stats::reshape(df, 
                        direction = "wide", 
                        idvar = id, 
                        timevar = within, 
                        v.names = value)
  
  # Clean up column names (remove the 'value.' prefix added by reshape)
  new_names <- names(res)
  prefix <- paste0(value, ".")
  idx_to_change <- grepl(paste0("^", prefix), new_names)
  new_names[idx_to_change] <- gsub(paste0("^", prefix), "", new_names[idx_to_change])
  names(res) <- new_names
  
  rownames(res) <- NULL
  return(res)
}
