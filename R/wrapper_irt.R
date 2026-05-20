#' RTMB-based IRT (Item Response Theory) Wrapper
#'
#' @description
#' Performs Item Response Theory modeling (1PL, 2PL, 3PL, and Graded Response Model) using RTMB.
#' Missing values (NA) in the data are automatically removed internally for efficient computation.
#'
#' @param data A data frame or matrix of item responses (N persons x J items).
#' @param model Character string for the model type: "1PL", "2PL", or "3PL".
#' @param type Character string for the data type: "binary" or "ordered".
#' @param prior Prior configuration: `prior_flat()`, `prior_normal()`, or `prior_weak()`.
#'   Hyperparameters can be specified within these functions (e.g., `prior_normal(b_sd = 5)`).
#'   Available parameters for IRT: `a_rate` (discrimination), `b_sd` (difficulty), `c_alpha`/`c_beta` (guessing).
#' @param init List of initial values.
#' @param fixed A named list of parameter values to fix (optional).
#' @param view Character vector of parameter names to prioritize in summary.
#' @param missing Missing value handling strategy: "listwise" (default) or "fiml" (Full Information Maximum Likelihood).
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @example inst/examples/ex_irt.R
#' @export
rtmb_irt <- function(data, model = c("2PL", "1PL", "3PL"), type = c("binary", "ordered"),
                     prior = prior_flat(), 
                     init = NULL, fixed = NULL, view = NULL,
                     missing = c("listwise", "fiml"), WAIC = FALSE) {

  missing <- match.arg(missing)

  model <- match.arg(model)
  type <- match.arg(type)

  if (model == "3PL" && type == "ordered") {
    stop("The 3PL model only supports 'binary' data.")
  }

  if (is.data.frame(data)) {
    if (!all(sapply(data, is.numeric))) {
      stop("All variables in the data must be numeric. Character or factor variables are not supported.", call. = FALSE)
    }
  } else if (!is.numeric(data) && !is.logical(data)) {
    stop("The data matrix must be numeric.", call. = FALSE)
  }

  Y <- as.matrix(data)
  
  if (missing == "listwise") {
    Y <- na.omit(Y)
  }

  if (type == "binary") {
    unique_vals <- na.omit(unique(as.vector(Y)))
    if (length(unique_vals) > 2 || any(!unique_vals %in% c(0, 1))) {
      stop("Binary IRT models (type = 'binary') require responses to be exactly 0, 1, or NA. Found other values.", call. = FALSE)
    }
  }

  # Get row and column names
  item_names <- colnames(Y)
  if (is.null(item_names)) item_names <- paste0("Item", 1:ncol(Y))
  person_names <- rownames(Y)
  if (is.null(person_names)) person_names <- paste0("Person", 1:nrow(Y))

  # Prior Handling
  if (is.null(prior)) prior <- prior_flat()

  if (!inherits(prior, "rtmb_prior")) {
    stop(
      "prior must be an object of class 'rtmb_prior'. ",
      "Use prior_flat(), prior_normal(), or prior_weak().",
      call. = FALSE
    )
  }

  prior_type <- prior$type
  
  # Extract hyperparameters from prior object if they exist
  a_rate  <- prior$a_rate
  b_sd    <- if (!is.null(prior$b_sd)) prior$b_sd else prior$mu_sd
  c_alpha <- prior$c_alpha
  c_beta  <- prior$c_beta

  # Default weak values
  if (prior_type == "weak") {
    if (is.null(a_rate)) a_rate <- 1
    if (is.null(b_sd)) b_sd <- 3
    if (is.null(c_alpha)) c_alpha <- 1
    if (is.null(c_beta)) c_beta <- 9
  }
  theta_sd <- 1 # Fixed to 1 for identification by default

  # --- Construction of Setup Block ---
  setup_exprs <- list()
  setup_exprs[[1]] <- quote(obs_data <- which(!is.na(Y), arr.ind = TRUE))
  setup_exprs[[2]] <- quote(person_idx <- as.integer(obs_data[, "row"]))
  setup_exprs[[3]] <- quote(item_idx <- as.integer(obs_data[, "col"]))
  setup_exprs[[4]] <- quote(Y_obs <- Y[obs_data])
  setup_exprs[[5]] <- quote(N_persons <- nrow(Y))
  setup_exprs[[6]] <- quote(N_items <- ncol(Y))
  setup_exprs[[7]] <- quote(N_obs <- length(Y_obs))

  if (type == "ordered") {
    setup_exprs[[8]] <- quote(if (min(Y_obs) == 0) Y_obs <- Y_obs + 1)
    setup_exprs[[9]] <- quote(K_cat <- max(Y_obs))
  }
  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(list(Y = Y))
  eval(setup_ast, tmp_env)

  # --- Construction of Parameters Block ---
  param_exprs <- list()
  param_exprs[[length(param_exprs) + 1]] <- quote(theta <- Dim(N_persons, random = TRUE))

  if (type == "binary") {
    param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(N_items))
  } else {
    param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(c(N_items, K_cat - 1), type = "ordered"))
  }

  if (model %in% c("2PL", "3PL")) {
    param_exprs[[length(param_exprs) + 1]] <- quote(a <- Dim(N_items, lower = 0))
  }
  if (model == "3PL") {
    param_exprs[[length(param_exprs) + 1]] <- quote(c <- Dim(N_items, lower = 0, upper = 1))
  }
  param_ast <- as.call(c(list(as.name("{")), param_exprs))

  # --- Construction of Model Block ---
  loop_body <- list()
  loop_body[[1]] <- quote(p <- person_idx[i])
  loop_body[[2]] <- quote(j <- item_idx[i])
  loop_body[[3]] <- quote(y <- Y_obs[i])

  if (type == "binary") {
    if (model == "1PL") {
      loop_body[[4]] <- quote(eta <- theta[p] - b[j])
    } else {
      loop_body[[4]] <- quote(eta <- a[j] * (theta[p] - b[j]))
    }
    if (model == "3PL") {
      loop_body[[5]] <- quote(prob <- c[j] + (1 - c[j]) * plogis(eta))
      loop_body[[6]] <- quote(y ~ bernoulli(prob))
    } else {
      loop_body[[5]] <- quote(y ~ bernoulli_logit(eta))
    }
  } else if (type == "ordered") {
    if (model == "1PL") {
      loop_body[[4]] <- quote(eta <- theta[p])
    } else {
      loop_body[[4]] <- quote(eta <- a[j] * theta[p])
    }
    loop_body[[5]] <- quote(y ~ ordered_logistic(eta, b[j, ]))
  }

  loop_ast <- quote(for (i in 1:N_obs) {})
  loop_ast[[4]] <- as.call(c(list(as.name("{")), loop_body))

  model_exprs <- list()
  model_exprs[[length(model_exprs) + 1]] <- "# Likelihood"
  model_exprs[[length(model_exprs) + 1]] <- loop_ast

  model_exprs[[length(model_exprs) + 1]] <- "# Priors"
  if (prior_type %in% c("normal", "weak")) {
    if (model %in% c("2PL", "3PL") && !is.null(a_rate)) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(a ~ exponential(.(a_rate)))
    }
    if (!is.null(b_sd)) {
      if (type == "binary") {
        model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(b_sd)))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(for (j in 1:N_items) b[j, ] ~ normal(0, .(b_sd)))
      }
    }
    if (model == "3PL" && !is.null(c_alpha) && !is.null(c_beta)) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(c ~ beta(.(c_alpha), .(c_beta)))
    }
  }
  model_exprs[[length(model_exprs) + 1]] <- bquote(theta ~ normal(0, .(theta_sd)))

  model_ast <- as.call(c(list(as.name("{")), model_exprs))

  # --- Mapping of Parameter Names ---
  code_obj <- list(setup = setup_ast, parameters = param_ast, model = model_ast)

  if (isTRUE(WAIC)) {
    gq_body <- list(
      quote(log_lik <- numeric(N_obs)),
      quote(for (i in 1:N_obs) {})
    )
    gq_loop_body <- list(
      quote(p <- person_idx[i]),
      quote(j <- item_idx[i]),
      quote(y <- Y_obs[i])
    )

    if (type == "binary") {
      if (model == "1PL") {
        gq_loop_body[[length(gq_loop_body) + 1]] <- quote(eta <- theta[p] - b[j])
      } else {
        gq_loop_body[[length(gq_loop_body) + 1]] <- quote(eta <- a[j] * (theta[p] - b[j]))
      }
      if (model == "3PL") {
        gq_loop_body[[length(gq_loop_body) + 1]] <- quote(prob <- c[j] + (1 - c[j]) * plogis(eta))
        gq_loop_body[[length(gq_loop_body) + 1]] <- quote(log_lik[i] <- bernoulli_lpmf(y, prob))
      } else {
        gq_loop_body[[length(gq_loop_body) + 1]] <- quote(log_lik[i] <- bernoulli_logit_lpmf(y, eta))
      }
    } else {
      if (model == "1PL") {
        gq_loop_body[[length(gq_loop_body) + 1]] <- quote(eta <- theta[p])
      } else {
        gq_loop_body[[length(gq_loop_body) + 1]] <- quote(eta <- a[j] * theta[p])
      }
      gq_loop_body[[length(gq_loop_body) + 1]] <- quote(log_lik[i] <- ordered_logistic_lpmf(y, eta, b[j, ]))
    }

    gq_body[[2]][[4]] <- as.call(c(list(as.name("{")), gq_loop_body))
    code_obj$generate <- .rtmb_waic_generate_ast(NULL, as.call(c(list(as.name("{")), gq_body)))
  }

  par_names_list <- list()
  par_names_list$theta <- person_names
  if (type == "binary") {
    par_names_list$b <- item_names
  } else {
    par_names_list$b <- list(item_names, paste0("Threshold", 1:(tmp_env$K_cat - 1)))
  }
  if (model %in% c("2PL", "3PL")) par_names_list$a <- item_names
  if (model == "3PL") par_names_list$c <- item_names

  view_vars <- c()
  if (model %in% c("2PL", "3PL")) view_vars <- c(view_vars, "a")
  view_vars <- c(view_vars, "b")
  if (model == "3PL") view_vars <- c(view_vars, "c")

  # --- Automatic Setting of Initial Values ---
  if (is.null(init)) {
    init <- list()
    if (type == "binary") {
      init$b <- rep(0, length(item_names))
    }
    if (model %in% c("2PL", "3PL")) init$a <- rep(1.0, length(item_names))
    if (model == "3PL") init$c <- rep(0.1, length(item_names))
  }

  obj <- rtmb_model(data = list(Y = Y), code = code_obj, par_names = par_names_list,
                    init = init, fixed = fixed, view = view_vars)
  obj$type <- "irt"
  obj$extra <- list(
    source = "wrapper",
    prior_type = prior$type,
    marginal = "theta"
  )

  return(obj)
}
