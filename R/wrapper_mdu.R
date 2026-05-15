#' Create Initial Values for Multidimensional Unfolding
#'
#' @param Y Numeric matrix or data frame (N rows x M items).
#' @param D Number of unfolding dimensions.
#' @param distance Character; `"squared"` or `"euclidean"`.
#' @param distance_eps Small positive constant added to the distance.
#' @param min_sigma Minimum initial residual standard deviation.
#'
#' @return A named list of initial values.
#' @export
make_init_mdu <- function(Y, D, distance = c("squared", "euclidean"),
                          distance_eps = 1e-4, min_sigma = 0.1) {
  distance <- match.arg(distance)
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"

  N <- nrow(Y)
  M <- ncol(Y)
  if (!is.numeric(D) || length(D) != 1L || D < 1) {
    stop("'D' must be a positive integer.", call. = FALSE)
  }
  D <- as.integer(D)
  if (D >= M) {
    stop("'D' must be smaller than the number of items/columns.", call. = FALSE)
  }

  row_means <- rowMeans(Y, na.rm = TRUE)
  col_means <- colMeans(Y, na.rm = TRUE)
  grand_mean <- mean(Y, na.rm = TRUE)
  row_means[!is.finite(row_means)] <- grand_mean
  col_means[!is.finite(col_means)] <- grand_mean

  Y_centered <- Y -
    matrix(row_means, N, M) -
    matrix(col_means, N, M, byrow = TRUE) +
    grand_mean
  Y_centered[!is.finite(Y_centered)] <- 0

  svd_res <- svd(Y_centered, nu = D, nv = D)
  sv <- pmax(svd_res$d[seq_len(D)], 0)
  D_mat <- diag(sqrt(sv), nrow = D, ncol = D)

  init_theta <- svd_res$u[, seq_len(D), drop = FALSE] %*% D_mat
  init_delta <- svd_res$v[, seq_len(D), drop = FALSE] %*% D_mat

  all_coords <- rbind(init_theta, init_delta)
  scale_factor <- stats::sd(as.vector(all_coords), na.rm = TRUE)
  if (!is.finite(scale_factor) || scale_factor <= 1e-8) scale_factor <- 1

  init_theta <- init_theta / scale_factor
  init_delta <- init_delta / scale_factor

  init_theta <- sweep(init_theta, 2, colMeans(init_theta), "-")
  init_delta <- sweep(init_delta, 2, colMeans(init_delta), "-")

  qr_res <- qr(t(init_delta[seq_len(D), , drop = FALSE]))
  Q <- qr.Q(qr_res)
  R <- qr.R(qr_res)
  sign_diag <- sign(diag(R))
  sign_diag[sign_diag == 0] <- 1
  Q <- sweep(Q, 2, sign_diag, "*")
  init_delta <- init_delta %*% Q
  init_theta <- init_theta %*% Q

  # Match Dim(type = "centered_tri"): zero above the diagonal and center
  # the free lower part of each dimension.
  init_delta_tri <- matrix(0, nrow = M, ncol = D)
  for (d in seq_len(D)) {
    vals <- init_delta[d:M, d]
    vals <- vals - mean(vals)
    init_delta_tri[d:M, d] <- vals
  }
  init_delta <- init_delta_tri
  init_theta <- sweep(init_theta, 2, colMeans(init_theta), "-")

  init_alpha <- max(Y, na.rm = TRUE)

  D_sq <- matrix(NA_real_, N, M)
  for (i in seq_len(N)) {
    for (j in seq_len(M)) {
      D_sq[i, j] <- sum((init_theta[i, ] - init_delta[j, ])^2)
    }
  }
  D_used <- if (distance == "squared") D_sq + distance_eps else sqrt(D_sq + distance_eps)

  diff_Y <- init_alpha - Y
  init_beta <- sum(diff_Y * D_used, na.rm = TRUE) / sum(D_used^2, na.rm = TRUE)
  if (!is.finite(init_beta) || init_beta <= 0) init_beta <- 1
  init_beta <- max(init_beta, 1e-4)

  mu <- init_alpha - init_beta * D_used
  init_sigma <- apply(Y - mu, 2, stats::sd, na.rm = TRUE)
  y_sd <- stats::sd(as.numeric(Y), na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
  init_sigma[!is.finite(init_sigma)] <- y_sd
  init_sigma <- pmax(init_sigma, min_sigma)

  list(
    alpha = init_alpha,
    beta = init_beta,
    sigma = init_sigma,
    theta = init_theta,
    delta = init_delta
  )
}

#' RTMB-based Multidimensional Unfolding Wrapper
#'
#' @description
#' Fits a multidimensional unfolding model for preference/rating data. Rows are
#' persons or observations and columns are stimuli/items. The model represents
#' both row scores (`theta`) and item locations (`delta`) in a shared D-dimensional
#' space.
#'
#' @param data Numeric matrix or data frame (N rows x M items).
#' @param ndim Number of unfolding dimensions.
#' @param distance Character; `"squared"` uses squared Euclidean distance
#'   (default, often easier for optimization), while `"euclidean"` uses Euclidean
#'   distance.
#' @param prior Prior configuration: `prior_flat()`, `prior_normal()`, or
#'   `prior_weak()`. `prior_flat()` creates a maximum-likelihood model suitable
#'   for `classic()`. The latent coordinates `delta` and `theta` are always
#'   treated as random effects with normal scale priors, similarly to IRT ability
#'   parameters.
#' @param y_range Optional response range. If supplied with the default flat
#'   prior, `prior_weak()` is used.
#' @param init Optional named list of initial values.
#' @param fixed Optional named list of parameter values to fix.
#' @param view Character vector of parameter names to prioritize in summaries.
#' @param distance_eps Small positive constant added to the distance.
#'
#' @return An `RTMB_Model` object.
#' @export
rtmb_mdu <- function(data, ndim = 2, distance = c("squared", "euclidean"),
                     prior = prior_flat(), y_range = NULL,
                     init = NULL, fixed = NULL, view = NULL,
                     distance_eps = 1e-4) {

  distance <- match.arg(distance)
  Y <- as.matrix(data)
  storage.mode(Y) <- "double"

  if (anyNA(Y)) {
    stop("rtmb_mdu() does not currently support missing values in 'data'.", call. = FALSE)
  }
  if (!is.numeric(ndim) || length(ndim) != 1L || ndim < 1) {
    stop("'ndim' must be a positive integer.", call. = FALSE)
  }
  D <- as.integer(ndim)
  N <- nrow(Y)
  M <- ncol(Y)
  if (D >= M) {
    stop("'ndim' must be smaller than the number of items/columns.", call. = FALSE)
  }

  item_names <- colnames(Y)
  if (is.null(item_names)) item_names <- paste0("Item", seq_len(M))
  person_names <- rownames(Y)
  if (is.null(person_names)) person_names <- paste0("Person", seq_len(N))
  dim_names <- paste0("Dim", seq_len(D))

  if (is.null(prior)) prior <- prior_flat()
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && identical(prior$type, "flat")) {
    prior <- prior_weak()
  }
  if (!inherits(prior, "rtmb_prior")) {
    stop(
      "prior must be an object of class 'rtmb_prior'. ",
      "Use prior_flat(), prior_normal(), or prior_weak().",
      call. = FALSE
    )
  }

  prior_type <- prior$type
  if (prior_type == "weak" && is.null(y_range)) {
    stop("When using prior_weak(), please specify 'y_range' (e.g., y_range = c(1, 7)).", call. = FALSE)
  }

  y_center <- mean(Y)
  y_sd <- stats::sd(as.numeric(Y))
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1

  alpha_center <- 0
  alpha_sd <- prior$alpha_sd %||% prior$mu_sd %||% 10
  beta_rate <- prior$beta_rate %||% 0.1
  sigma_rate <- prior$sigma_rate %||% 0.1
  delta_sd <- prior$delta_sd %||% prior$tau_sd %||% 1
  theta_sd <- prior$theta_sd %||% 1

  if (prior_type == "weak") {
    half_range <- diff(y_range) / 2
    if (!is.finite(half_range) || half_range <= 0) {
      stop("'y_range' must be an increasing numeric vector of length 2.", call. = FALSE)
    }
    sd_ratio <- prior$sd_ratio %||% 0.5
    alpha_center <- mean(y_range)
    alpha_sd <- half_range
    sigma_rate <- 1 / (half_range * sd_ratio)
    beta_rate <- 1 / (prior$max_beta %||% 1)
    delta_sd <- prior$delta_sd %||% 1
    theta_sd <- prior$theta_sd %||% 1
  }

  dat_mdu <- list(
    Y = Y,
    D = D,
    distance_eps = distance_eps,
    alpha_center = alpha_center,
    alpha_sd = alpha_sd,
    beta_rate = beta_rate,
    sigma_rate = sigma_rate,
    delta_sd = delta_sd,
    theta_sd = theta_sd
  )

  use_random_coords <- TRUE

  setup_ast <- quote({
    N <- nrow(Y)
    M <- ncol(Y)
  })

  param_ast <- bquote({
    alpha <- Dim()
    beta <- Dim(lower = 0)
    sigma <- Dim(M, lower = 0)
    delta <- Dim(c(M, D), type = "centered_tri", random = TRUE)
    theta <- Dim(c(N, D), random = TRUE)
  })

  loop_ast <- if (distance == "squared") {
    quote(for (i in 1:N) {
      for (j in 1:M) {
        d_ij <- squared_distance(theta[i, ], delta[j, ]) + distance_eps
        mu_ij <- alpha - beta * d_ij
        Y[i, j] ~ normal(mu_ij, sigma[j])
      }
    })
  } else {
    quote(for (i in 1:N) {
      for (j in 1:M) {
        d_ij <- distance(theta[i, ], delta[j, ]) + distance_eps
        mu_ij <- alpha - beta * d_ij
        Y[i, j] ~ normal(mu_ij, sigma[j])
      }
    })
  }

  model_exprs <- list(loop_ast)
  model_exprs[[length(model_exprs) + 1]] <- quote(delta ~ centered_tri_multi_normal(delta_sd))
  model_exprs[[length(model_exprs) + 1]] <- quote(for (d in 1:D) theta[, d] ~ normal(0, theta_sd))
  if (prior_type %in% c("normal", "weak")) {
    model_exprs[[length(model_exprs) + 1]] <- quote(alpha ~ normal(alpha_center, alpha_sd))
    model_exprs[[length(model_exprs) + 1]] <- quote(beta ~ exponential(beta_rate))
    model_exprs[[length(model_exprs) + 1]] <- quote(sigma ~ exponential(sigma_rate))
  }
  model_ast <- as.call(c(list(as.name("{")), model_exprs))

  code_obj <- list(setup = setup_ast, parameters = param_ast, model = model_ast, env = parent.frame())

  if (is.null(init)) {
    init <- make_init_mdu(Y, D, distance = distance, distance_eps = distance_eps)
  }

  par_names_list <- list(
    sigma = item_names,
    delta = list(item_names, dim_names),
    theta = list(person_names, dim_names)
  )

  if (is.null(view)) view <- c("alpha", "beta", "sigma", "delta")

  obj <- rtmb_model(
    data = dat_mdu,
    code = code_obj,
    par_names = par_names_list,
    init = init,
    fixed = fixed,
    view = view
  )
  obj$type <- "mdu"
  obj$extra <- list(
    source = "wrapper",
    prior_type = prior$type,
    distance = distance,
    marginal = c("delta", "theta")
  )
  obj
}
