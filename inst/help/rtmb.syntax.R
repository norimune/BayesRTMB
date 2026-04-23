#' Guidelines for Writing RTMB-Compatible Code
#'
#' @description
#' Models defined in \code{rtmb_code} rely on Automatic Differentiation (AD) via the
#' \code{RTMB} package. To ensure the model is differentiable and numerically stable,
#' specific coding practices must be followed.
#'
#' @details
#' \strong{1. Automatic Differentiation and \code{advector}:}
#' Parameters and intermediate calculations involving them are treated as \code{advector}
#' objects. RTMB "records" all mathematical operations to compute derivatives
#' automatically. If a function strips these attributes or is not differentiable, the
#' gradient calculation will fail.
#'
#' \strong{2. Avoid Discrete Branching:}
#' Standard R conditional statements like \code{if (x > 0)} or \code{ifelse()} based
#' on parameter values do not provide derivatives.
#' \itemize{
#'   \item \strong{Problem:} They create "jumps" in the likelihood surface.
#'   \item \strong{Solution:} Use smooth approximations. For example, use \code{\link{fabs}}
#'   instead of \code{abs}, or \code{\link{log_mix}} for mixture logic.
#' }
#'
#' \strong{3. Numerical Stability:}
#' Computations in log-space are preferred to prevent overflow or underflow.
#' Use the following stable utilities instead of raw algebraic expressions:
#' \itemize{
#'   \item \code{\link{log_sum_exp}}: For summing probabilities in log-space.
#'   \item \code{\link{log1p_exp}}: For \code{log(1 + exp(x))}.
#'   \item \code{\link{inv_logit}}: For mapping real numbers to probabilities.
#' }
#'
#' \strong{4. Vectorization vs. Loops:}
#' \itemize{
#'   \item \strong{Vectorization:} Highly recommended for performance. Standard R
#'   vectorized arithmetic (\code{+}, \code{-}, \code{*}, \code{/}, \code{log}, \code{exp})
#'   works seamlessly with AD.
#'   \item \strong{Loops:} Standard \code{for} loops are safe as long as the operations
#'   inside are differentiable.
#'   \item \strong{Avoid \code{apply}:} Functions like \code{apply}, \code{sapply},
#'   or \code{lapply} may sometimes strip AD attributes and should be replaced with
#'   vectorized operations or explicit loops.
#' }
#'
#' \strong{5. Matrix Operations:}
#' Use specialized functions for matrix algebra to maintain efficiency:
#' \itemize{
#'   \item \code{\link{quad_form_chol}}: For efficient quadratic forms.
#'   \item \code{\link{log_det_chol}}: For log-determinants via Cholesky factors.
#' }
#'
#' @name rtmb_syntax
#' @seealso \code{\link{math_functions}}, \code{\link{distributions}}, \code{\link{rtmb_code}}
NULL
