#' Find random effects terms (bars) in a formula
#'
#' @param term A formula or expression.
#' @return A list of bar terms.
#' @importFrom stats as.formula model.frame model.matrix model.response model.offset
#' @keywords internal
#' @noRd
findbars <- function(term) {
  if (is.name(term) || !is.call(term)) return(NULL)
  if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
  if (term[[1]] == as.name("|") || term[[1]] == as.name("||")) return(list(term))
  res <- lapply(term, findbars)
  return(unlist(res, recursive = FALSE))
}

#' Remove random effects terms from a formula
#'
#' @param term A formula or expression.
#' @return A formula or expression with bars removed.
#' @keywords internal
#' @noRd
nobars <- function(term) {
  if (!is.language(term)) return(term)
  if (is.call(term) && (term[[1]] == as.name("|") || term[[1]] == as.name("||"))) return(NULL)
  if (length(term) == 1) return(term)
  if (length(term) == 2) {
    nb <- nobars(term[[2]])
    if (is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars(term[[2]])
  nb3 <- nobars(term[[3]])
  if (is.null(nb3)) {
    if (is.call(term) && term[[1]] == as.name("~")) {
      term[[2]] <- nb2
      term[[3]] <- as.name("1")
      return(term)
    }
    return(nb2)
  }
  if (is.null(nb2)) return(nb3)
  term[[2]] <- nb2
  term[[3]] <- nb3
  return(term)
}

#' Replace bars with pluses to get all variables
#'
#' @param term A formula or expression.
#' @return A formula or expression with bars replaced by pluses.
#' @keywords internal
#' @noRd
subbars <- function(term) {
  if (is.name(term) || !is.call(term)) return(term)
  if (term[[1]] == as.name("(")) return(subbars(term[[2]]))
  if (term[[1]] == as.name("|") || term[[1]] == as.name("||")) {
    return(as.call(list(as.name("+"), subbars(term[[2]]), subbars(term[[3]]))))
  }
  if (length(term) == 1) return(term)
  res <- lapply(term, subbars)
  return(as.call(res))
}
