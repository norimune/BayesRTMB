#' @importFrom stats .getXlevels ave chisq.test fisher.test kmeans lm na.pass optim prcomp qt rchisq setNames update
#' @importFrom utils write.csv
NULL

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# Add global variables to avoid R CMD check notes
utils::globalVariables(c("N", "Y", "iteration", "chain", "value", "variable"))
