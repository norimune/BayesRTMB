#' RTMB-based Contingency Table Analysis (Chi-squared Test)
#'
#' @description
#' `rtmb_table` performs a chi-squared test of independence between two categorical variables.
#' It provides both classic (frequentist) Pearson chi-squared tests and Bayesian multinomial-style models.
#'
#' @param x Variable name or formula.
#' @param y Variable name (optional if x is a formula).
#' @param data A data frame.
#' @param classic Logical; if TRUE, perform frequentist chi-squared and Fisher's exact tests.
#' @param correct Logical; if TRUE, apply Yates' continuity correction (for 2x2 classic only).
#' @param prior Prior specification (Bayesian mode). Default is `prior_uniform()`.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param ... Additional arguments.
#'
#' @return A `Classic_Fit` or `MCMC_Fit` object.
#'
#' @examples
#' \donttest{
#' # Classic chi-squared test
#' rtmb_table(skill, cond, data = debate, classic = TRUE)
#' }
#' @export
rtmb_table <- function(x, y = NULL, data = NULL, classic = FALSE, correct = TRUE, prior = prior_uniform(), fixed = NULL, ...) {

  x_expr <- substitute(x)
  y_expr <- substitute(y)

  # 1. Extract Variables
  if (is.null(data)) {
    v1 <- eval(x_expr, parent.frame())
    v2 <- eval(y_expr, parent.frame())
    v1_name <- deparse(x_expr)
    v2_name <- deparse(y_expr)
  } else {
    v1 <- if (is.name(x_expr)) data[[as.character(x_expr)]] else eval(x_expr, data)
    v2 <- if (is.name(y_expr)) data[[as.character(y_expr)]] else eval(y_expr, data)
    v1_name <- as.character(x_expr)
    v2_name <- as.character(y_expr)
  }

  # Ensure factors
  v1 <- as.factor(v1)
  v2 <- as.factor(v2)
  tab <- table(v1, v2)

  # Prepare row/column names for labels
  R_names <- rownames(tab)
  C_names <- colnames(tab)
  grid <- expand.grid(Row = R_names, Col = C_names)
  cell_labels <- paste0(v1_name, ":", grid$Row, ", ", v2_name, ":", grid$Col)

  # 2. Classic Mode (Frequentist)
  if (isTRUE(classic)) {
    chisq_res <- stats::chisq.test(tab, correct = correct)
    fisher_res <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    
    res <- Classic_Fit$new(
      model = list(type = "table", data = list(tab = tab), par_names = list(p = cell_labels)),
      df_fixed = data.frame(
        Statistic = c("Chi-squared", if (!is.null(fisher_res)) "Fisher's Exact" else NULL),
        Estimate = c(chisq_res$statistic, if (!is.null(fisher_res)) NA else NULL),
        df = c(chisq_res$parameter, if (!is.null(fisher_res)) NA else NULL),
        Pr = c(chisq_res$p.value, if (!is.null(fisher_res)) fisher_res$p.value else NULL),
        check.names = FALSE
      ),
      test_results = list(chisq = chisq_res, fisher = fisher_res)
    )
    class(res) <- c("rtmb_table", class(res))
    res$type <- "table"
    res$extra <- list(tab = tab, correct = correct)
    return(res)
  }

  # 3. Bayes Mode (Multinomial Model)
  # Data for RTMB
  Y_vec <- as.vector(tab)
  R <- nrow(tab)
  C <- ncol(tab)
  N_total <- sum(Y_vec)

  setup <- list(
    Y = Y_vec,
    R = R,
    C = C,
    N = N_total
  )

  rtmb_model_code <- rtmb_code(
    setup = {
      # No special setup needed for now, Y and N are provided
    },
    parameters = {
      # Simplex for probabilities (automatically constrained to sum to 1)
      # We give it meaningful names during rtmb_model() call
      p <- Dim(R * C, type = "simplex")
    },
    transform = {
      # Derived quantities for reporting
      mu <- p * N
      # Pearson Chi-squared statistic for the estimated probabilities
      # Sum (O - E)^2 / E
      chisq_val <- sum((Y - mu)^2 / mu)
    },
    model = {
      # Likelihood
      Y ~ multinomial(N, p)

      # Prior (Flat Dirichlet by default)
      p ~ dirichlet(rep(1, R * C))
    },
    generate = {
      # mu and chisq_val are already reported in transform
    }
  )

  # Create the model object with explicit parameter names for 'p' and 'mu'
  res <- rtmb_model(
    code = rtmb_model_code,
    data = setup,
    par_names = list(p = cell_labels, mu = cell_labels),
    fixed = fixed
  )
  class(res) <- c("rtmb_table", class(res))

  res$type <- "table"
  res$extra <- list(tab = tab, correct = correct)
  return(res)
}
