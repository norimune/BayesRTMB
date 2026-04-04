#' Mathematical helper functions
#'
#' `math` is an environment containing mathematical helper functions
#' used to build `log_prob()` functions.
#'
#' @export
math <- new.env(parent = baseenv())

math$inv_logit <- function(x) plogis(x)

math$logit <- function(x) qlogis(x)

math$distance <- function(x, y, eps = 1e-8) {
  sqrt(sum((x - y)^2) + eps)
}

math$squared_distance <- function(x, y, eps = 1e-8) {
  sum((x - y)^2) + eps
}

math$log_sum_exp <- function(x, y = NULL) {
  if (is.null(y)) {
    max_val <- max(as.numeric(x))
    return(max_val + log(sum(exp(x - max_val))))
  } else {
    max_val <- (x + y + abs(x - y)) / 2
    return(max_val + log(exp(x - max_val) + exp(y - max_val)))
  }
}

math$log1m <- function(x) {
  log1p(-x)
}

math$log1p_exp <- function(x) {
  max_val <- (x + sqrt(x^2 + 1e-10)) / 2
  return(max_val + log(exp(x - max_val) + exp(-max_val)))
}

math$log1m_exp <- function(x) {
  return(log(1 - exp(x) + 1e-10))
}

math$log_mix <- function(theta, lp1, lp2) {
  # log_sum_expを使って安定的に計算
  l1 <- log(theta) + lp1
  l2 <- log(1 - theta) + lp2
  return(math$log_sum_exp(l1, l2))
}

math$softmax <- function(x) {
  max_x <- max(x)
  exp_x <- exp(x - max_x)
  return(exp_x / sum(exp_x))
}

math$log_softmax <- function(x) {
  return(x - math$log_sum_exp(x))
}

math$fabs <- function(x, eps = 1e-8) {
  sqrt(x^2 + eps)
}

math$log_det_chol <- function(L) {
  return(2 * sum(log(diag(L))))
}

math$quad_form_chol <- function(x, L) {
  z <- solve(L, x)
  return(sum(z^2))
}

math$quad_form_diag <- function(m, v) {
  v_col <- matrix(as.vector(v), ncol = 1)
  return(m * (v_col %*% t(v_col)))
}
