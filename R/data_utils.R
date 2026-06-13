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
#'     \item A list of column-name vectors or range/prefix specifications to
#'       concatenate into one long value column.
#'     \item Multiple range/prefix/name specifications when \code{value} has the same length.
#'   }
#' @param label The name for the new column that will contain the level names (e.g., "Time"). 
#'   If NULL, it defaults to the prefix used in \code{within} or "Condition".
#' @param value The name for the new column that will contain the measurement values.
#'   Can be a character vector for multiple measurement groups. Default is "Value".
#' @param id The identifier columns that should be repeated for each row. 
#'   If NULL (default), all columns NOT specified in \code{within} are treated as IDs.
#' @param sort Logical; if TRUE, sort the output by \code{id} and the within
#'   label. If FALSE (default), preserve the original row order of \code{data}
#'   and the within-column order inside each row.
#' @return A data frame in long format.
#' @export
to_long <- function(data, within = NULL, label = NULL, value = "Value", id = NULL, sort = FALSE) {
  if (is.null(within)) {
    stop("Argument 'within' is required. Specify columns to gather (e.g., 'time1:time4' or 'time').")
  }
  if (!is.logical(sort) || length(sort) != 1L || is.na(sort)) {
    stop("'sort' must be TRUE or FALSE.", call. = FALSE)
  }
  
  all_names <- names(data)
  target_cols <- NULL
  inferred_label <- "Condition"

  resolve_within <- function(spec) {
    inferred <- "Condition"

    if (is.character(spec) && length(spec) == 1) {
      if (grepl(":", spec)) {
        # Range notation: "time1:time4"
        parts <- strsplit(spec, ":")[[1]]
        start_idx <- which(all_names == parts[1])
        end_idx <- which(all_names == parts[2])
        if (length(start_idx) == 1 && length(end_idx) == 1) {
          cols <- all_names[start_idx:end_idx]
        } else {
          stop(sprintf("Range columns '%s' or '%s' not found.", parts[1], parts[2]))
        }
      } else if (spec %in% all_names) {
        # Single column name
        cols <- spec
        inferred <- spec
      } else {
        # Prefix matching: "time" -> time1, time2, etc.
        cols <- all_names[grepl(paste0("^", spec), all_names)]
        if (length(cols) > 0) {
          inferred <- spec
        } else {
          stop(sprintf("No columns found starting with prefix '%s'.", spec))
        }
      }
    } else {
      # Vector of names or indices
      if (is.numeric(spec)) {
        cols <- all_names[spec]
      } else {
        cols <- spec
      }
    }

    list(cols = cols, label = inferred)
  }

  # 1. Handle 'within' specification
  if (length(value) > 1) {
    within_specs <- if (is.list(within)) within else as.list(within)
    if (length(within_specs) != length(value)) {
      stop(sprintf(
        "When using multiple value columns, 'within' and 'value' must have the same length (%d vs %d).",
        length(within_specs), length(value)
      ))
    }

    resolved <- lapply(within_specs, resolve_within)
    target_cols <- lapply(resolved, `[[`, "cols")
    inferred_labels <- vapply(resolved, `[[`, character(1), "label")
    inferred_label <- inferred_labels[1]

    group_lengths <- lengths(target_cols)
    if (length(unique(group_lengths)) != 1) {
      warning(
        sprintf(
          "All 'within' groups should have the same number of columns; got lengths: %s.",
          paste(group_lengths, collapse = ", ")
        ),
        call. = FALSE
      )
      stop("Cannot reshape multiple 'within' groups unless they have the same number of columns.", call. = FALSE)
    }
  } else if (is.list(within)) {
    resolved <- lapply(within, resolve_within)
    target_cols <- unlist(lapply(resolved, `[[`, "cols"), use.names = FALSE)
    inferred_labels <- vapply(resolved, `[[`, character(1), "label")
    non_default_labels <- inferred_labels[inferred_labels != "Condition"]
    inferred_label <- if (length(non_default_labels) > 0L) non_default_labels[1] else "Condition"
  } else {
    resolved <- resolve_within(within)
    target_cols <- resolved$cols
    inferred_label <- resolved$label
  }
  
  if (is.null(target_cols) || length(target_cols) == 0) {
    stop("Could not identify columns to gather based on 'within' specification.")
  }
  if (is.list(target_cols) && any(lengths(target_cols) == 0)) {
    stop("Could not identify columns to gather based on one or more 'within' specifications.")
  }

  # 2. Determine label column name
  if (is.null(label)) {
    label <- inferred_label
  }

  # 3. Determine ID columns (everything not in target_cols)
  if (is.null(id)) {
    id <- setdiff(all_names, unlist(target_cols, use.names = FALSE))
  }

  # 4. Perform reshape
  data_res <- as.data.frame(data)
  order_col <- "..rtmb_row_order"
  while (order_col %in% names(data_res)) {
    order_col <- paste0(".", order_col)
  }
  data_res[[order_col]] <- seq_len(nrow(data_res))
  label_levels <- if (is.list(target_cols)) target_cols[[1]] else target_cols
  
  # stats::reshape requires a list of varying column names
  res <- stats::reshape(data_res, 
                        direction = "long", 
                        varying = if (is.list(target_cols)) target_cols else list(target_cols), 
                        v.names = value, 
                        timevar = label, 
                        times = label_levels, 
                        idvar = c(id, order_col))
  
  # Ensure the label column is a factor with levels in original column order
  res[[label]] <- factor(res[[label]], levels = label_levels)

  if (sort) {
    # Sort by ID and label for readability.
    res <- res[do.call(order, res[c(id, label)]), ]
  } else {
    # Preserve the source row order, then the within-column order inside each row.
    res <- res[do.call(order, res[c(order_col, label)]), ]
  }
  res[[order_col]] <- NULL
  rownames(res) <- NULL
  
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
#'   Can be a character vector for multiple measurement columns.
#' @param id The identifier columns (e.g., "Subject", "Group"). If NULL, 
#'   inferred as all columns except \code{within} and \code{value}.
#' @return A data frame in wide format.
#' @export
to_wide <- function(data, within = "Condition", value = "Value", id = NULL) {
  df <- as.data.frame(data)
  if (!within %in% names(df)) {
    stop(sprintf("Column '%s' specified in 'within' was not found.", within), call. = FALSE)
  }
  missing_value <- setdiff(value, names(df))
  if (length(missing_value) > 0) {
    stop(
      sprintf("Column(s) specified in 'value' not found: %s.", paste(missing_value, collapse = ", ")),
      call. = FALSE
    )
  }

  if (is.null(id)) {
    id <- setdiff(names(df), c(within, value))
  }

  res <- stats::reshape(df, 
                        direction = "wide", 
                        idvar = id, 
                        timevar = within, 
                        v.names = value)
  
  if (length(value) == 1) {
    # Clean up column names (remove the 'value.' prefix added by reshape)
    new_names <- names(res)
    prefix <- paste0(value, ".")
    idx_to_change <- grepl(paste0("^", prefix), new_names)
    new_names[idx_to_change] <- gsub(paste0("^", prefix), "", new_names[idx_to_change])
    names(res) <- new_names
  }
  
  rownames(res) <- NULL
  return(res)
}
