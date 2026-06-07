#' RTMB-based Contingency Table Analysis (Chi-squared Test)
#'
#' @description
#' `rtmb_table` performs a chi-squared test of independence between two categorical variables.
#' It provides both classic (frequentist) Pearson chi-squared tests and Bayesian multinomial-style models.
#'
#' @param x Variable name, formula, table, or matrix.
#' @param y Variable name (optional if x is a formula).
#' @param data A data frame.
#' @param correct Logical; if TRUE, apply Yates' continuity correction for 2x2 classic analyses.
#' @param prior Prior specification (Bayesian mode). Use `prior_flat()` for a
#'   uniform Dirichlet prior or `prior_normal(dirichlet_alpha = ...)` to set the
#'   Dirichlet concentration. Default is `prior_flat()`.
#' @param fixed Optional named list of fixed values for specific parameters.
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @param ... Reserved; unused arguments are rejected.
#'
#' @return An `RTMB_Model` object.
#'
#' @examples
#' \donttest{
#' # Classic chi-squared test
#' rtmb_table(skill, cond, data = debate)$classic()
#' rtmb_table(table(debate$skill, debate$cond))$classic()
#' }
#' @export
rtmb_table <- function(x, y = NULL, data = NULL, correct = TRUE, prior = prior_flat(),
                       fixed = NULL, WAIC = FALSE, ...) {

  .check_unused_dots(..., .fn = "rtmb_table()")
  x_expr <- substitute(x)
  y_expr <- substitute(y)
  prior <- .validate_prior_type(
    prior,
    allowed = c("flat", "normal"),
    context = "rtmb_table()"
  )

  # 1. Extract table or variables
  if (!is.null(data)) {
    x_val <- if (is.name(x_expr)) data[[as.character(x_expr)]] else eval(x_expr, data)
  } else {
    x_val <- eval(x_expr, parent.frame())
  }

  if (inherits(x_val, "table") || (is.matrix(x_val) && is.numeric(x_val))) {
    tab <- as.table(x_val)
    if (length(dim(tab)) != 2L) {
      stop("`x` must be a two-way table or matrix when supplied as a table-like object.", call. = FALSE)
    }
    dimn <- dimnames(tab)
    if (is.null(dimn)) dimn <- vector("list", 2L)
    if (is.null(dimn[[1]])) dimn[[1]] <- paste0("Row", seq_len(nrow(tab)))
    if (is.null(dimn[[2]])) dimn[[2]] <- paste0("Col", seq_len(ncol(tab)))
    dimnames(tab) <- dimn
    table_name <- deparse(x_expr)
    dim_names <- names(dimn)
    v1_name <- if (!is.null(dim_names) && nzchar(dim_names[1])) dim_names[1] else paste0(table_name, "_row")
    v2_name <- if (!is.null(dim_names) && length(dim_names) >= 2L && nzchar(dim_names[2])) dim_names[2] else paste0(table_name, "_col")
  } else {
    if (missing(y) || identical(y_expr, quote(NULL))) {
      stop("`y` is required unless `x` is a two-way table or matrix.", call. = FALSE)
    }
    if (is.null(data)) {
      v1 <- x_val
      v2 <- eval(y_expr, parent.frame())
      v1_name <- deparse(x_expr)
      v2_name <- deparse(y_expr)
    } else {
      v1 <- x_val
      v2 <- if (is.name(y_expr)) data[[as.character(y_expr)]] else eval(y_expr, data)
      v1_name <- as.character(x_expr)
      v2_name <- as.character(y_expr)
    }

    v1 <- as.factor(v1)
    v2 <- as.factor(v2)
    tab <- table(v1, v2)
  }

  # Prepare row/column names for labels
  R_names <- rownames(tab)
  C_names <- colnames(tab)
  grid <- expand.grid(Row = R_names, Col = C_names)
  cell_labels <- paste0(v1_name, ":", grid$Row, ", ", v2_name, ":", grid$Col)

  # 2. Model (Multinomial Model)
  # Data for RTMB
  Y_vec <- as.vector(tab)
  R <- nrow(tab)
  C <- ncol(tab)
  N_total <- sum(Y_vec)
  if (R < 2L || C < 2L) {
    stop("rtmb_table() requires a table with at least two rows and two columns.", call. = FALSE)
  }
  if (any(is.na(Y_vec)) || any(!is.finite(Y_vec))) {
    stop("Table counts must be finite and non-missing.", call. = FALSE)
  }
  if (any(Y_vec < 0)) {
    stop("Table counts must be non-negative.", call. = FALSE)
  }
  if (any(abs(Y_vec - round(Y_vec)) > .Machine$double.eps^0.5)) {
    stop("Table counts must be whole numbers.", call. = FALSE)
  }
  if (N_total <= 0) {
    stop("Table counts must sum to a positive total.", call. = FALSE)
  }

  setup <- list(
    Y = Y_vec,
    R = R,
    C = C,
    N = N_total
  )

  # Pre-resolve prior settings (avoid runtime NSE dependency on prior object)
  dirichlet_alpha_val <- if (inherits(prior, "rtmb_prior") && prior$type == "normal" &&
                              !is.null(prior$dirichlet_alpha)) {
    prior$dirichlet_alpha
  } else {
    1
  }
  setup$dirichlet_alpha <- dirichlet_alpha_val

  gen_ast <- if (isTRUE(WAIC)) {
    .rtmb_waic_generate_ast(NULL, quote({
      log_lik <- multinomial_lpmf(Y, N, p)
    }))
  } else {
    quote({
      # mu and chisq_val are already reported in transform
    })
  }

  rtmb_model_code <- eval(substitute(
    rtmb_code(
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

        # Prior (alpha value is embedded as data at construction time)
        p ~ dirichlet(rep(dirichlet_alpha, R * C))
      },
      generate = G
    ),
    list(G = gen_ast)
  ))

  # Create the model object with explicit parameter names for 'p' and 'mu'
  res <- rtmb_model(
    code = rtmb_model_code,
    data = setup,
    par_names = list(p = cell_labels, mu = cell_labels),
    fixed = fixed
  )
  # Set metadata for special print/summary dispatch
  res$type <- "table"
  res$extra <- list(
    source = "wrapper",
    prior_type = if (inherits(prior, "rtmb_prior")) prior$type else "flat",
    marginal = "p",
    tab = tab, 
    correct = correct
  )
  
  class(res) <- c("rtmb_table", class(res))
  return(res)
}
