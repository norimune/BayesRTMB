# Helper to handle wide-to-long conversion for repeated measures
.handle_wide_to_long <- function(formula, data, within, factors) {
  lhs <- formula[[2]]
  if (!(is.call(lhs) && identical(lhs[[1]], as.name("cbind")))) {
    return(list(formula = formula, data = data, factors = factors))
  }

  # 1. Expand columns in cbind()
  # Support cbind(v1:v4) or cbind(v1, v2, ...)
  data_names <- names(data)
  resp_cols <- character(0)
  
  for (i in 2:length(lhs)) {
    arg <- lhs[[i]]
    if (is.call(arg) && identical(arg[[1]], as.name(":"))) {
      start_name <- as.character(arg[[2]])
      end_name <- as.character(arg[[3]])
      idx_start <- which(data_names == start_name)
      idx_end <- which(data_names == end_name)
      if (length(idx_start) == 0 || length(idx_end) == 0) {
        stop(sprintf("Range columns '%s' or '%s' not found in data.", start_name, end_name))
      }
      resp_cols <- c(resp_cols, data_names[idx_start:idx_end])
    } else {
      resp_cols <- c(resp_cols, as.character(arg))
    }
  }
  
  # Validate columns exist
  missing_cols <- setdiff(resp_cols, data_names)
  if (length(missing_cols) > 0) {
    stop(sprintf("Response columns not found in data: %s", paste(missing_cols, collapse = ", ")))
  }

  n_cols <- length(resp_cols)

  # 2. Handle 'within' inference or validation
  if (is.null(within)) {
    # Infer from RHS
    rhs_vars <- all.vars(formula[[3]])
    missing_vars <- setdiff(rhs_vars, data_names)
    
    if (length(missing_vars) == 1) {
      within <- setNames(list(n_cols), missing_vars)
    } else if (length(missing_vars) > 1) {
      stop(sprintf("Multiple variables in RHS are not in data (%s). Please specify 'within' list.", paste(missing_vars, collapse = ", ")))
    } else {
      stop("LHS is cbind() but no within-factor found in RHS that is not in data. Please specify 'within'.")
    }
  }

  n_levels <- prod(as.numeric(within))
  if (n_cols != n_levels) {
    stop(sprintf("Number of response columns (%d) does not match total levels in 'within' (%d).", n_cols, n_levels))
  }

  # 3. Create within-subjects grid (last factor changes fastest)
  within_list <- lapply(within, function(x) seq_len(x))
  within_grid <- expand.grid(rev(within_list))
  within_grid <- within_grid[, rev(seq_along(within)), drop = FALSE]
  colnames(within_grid) <- names(within)

  # 4. Reshape data to long format
  N_rows <- nrow(data)
  
  # Repeat rows of original data
  long_data <- data[rep(seq_len(N_rows), each = n_cols), , drop = FALSE]
  # Remove wide columns
  long_data <- long_data[, setdiff(names(long_data), resp_cols), drop = FALSE]
  
  # Add within-factors
  for (fn in names(within)) {
    long_data[[fn]] <- rep(within_grid[[fn]], N_rows)
  }
  
  # Extract response values row by row
  Y_wide <- as.matrix(data[, resp_cols, drop = FALSE])
  long_data$RTMB_Response <- as.vector(t(Y_wide))
  
  # 5. Remove NA trials
  long_data <- long_data[!is.na(long_data$RTMB_Response), , drop = FALSE]
  if (nrow(long_data) == 0) stop("All observations removed due to NA in response columns.")

  # 6. Update formula and factors
  new_formula <- formula
  new_formula[[2]] <- as.name("RTMB_Response")
  
  new_factors <- unique(c(factors, names(within)))

  return(list(formula = new_formula, data = long_data, factors = new_factors))
}
