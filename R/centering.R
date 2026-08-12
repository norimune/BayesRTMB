#' Center a Numeric Variable Around Its Grand Mean
#'
#' @description
#' Subtracts the overall mean from a numeric vector. This helper is used by the
#' regression wrappers in their generated `setup` code and can also be used in
#' hand-written `rtmb_code()` models.
#'
#' @param x Numeric vector to center.
#' @param na.rm Logical; whether missing values are removed when calculating the
#'   mean.
#'
#' @return A numeric vector with the same length as `x`.
#' @export
#'
#' @examples
#' center_grand_mean(c(1, 2, 3))
center_grand_mean <- function(x, na.rm = TRUE) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("'x' must be a numeric vector.", call. = FALSE)
  }
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
  }

  x - mean(x, na.rm = na.rm)
}

#' Center a Numeric Variable Within Clusters
#'
#' @description
#' Subtracts each cluster mean from a numeric vector. Missing cluster values
#' produce missing centered values. This helper is used by the regression
#' wrappers in their generated `setup` code and can also be used in hand-written
#' `rtmb_code()` models.
#'
#' @param x Numeric vector to center.
#' @param cluster Vector identifying the cluster for each element of `x`.
#' @param na.rm Logical; whether missing values are removed when calculating
#'   cluster means.
#'
#' @return A numeric vector with the same length as `x`.
#' @export
#'
#' @examples
#' center_within_cluster(c(1, 3, 2, 6), c("a", "a", "b", "b"))
center_within_cluster <- function(x, cluster, na.rm = TRUE) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("'x' must be a numeric vector.", call. = FALSE)
  }
  if (!is.atomic(cluster) || !is.null(dim(cluster))) {
    stop("'cluster' must be a vector.", call. = FALSE)
  }
  if (length(x) != length(cluster)) {
    stop("'x' and 'cluster' must have the same length.", call. = FALSE)
  }
  if (!is.logical(na.rm) || length(na.rm) != 1L || is.na(na.rm)) {
    stop("'na.rm' must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(x) == 0L) return(as.numeric(x))

  group_mean <- rep(NA_real_, length(x))
  observed_cluster <- !is.na(cluster)
  cluster_rows <- split(
    which(observed_cluster),
    as.character(cluster[observed_cluster]),
    drop = TRUE
  )
  for (rows in cluster_rows) {
    group_mean[rows] <- mean(x[rows], na.rm = na.rm)
  }

  x - group_mean
}
