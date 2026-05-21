#' Create Initial Values for Multidimensional Unfolding
#'
#' @param Y Numeric matrix or data frame (N rows x M items).
#' @param D Number of unfolding dimensions.
#' @param distance Character; `"squared"` or `"euclidean"`.
#' @param alpha Character; `"random"` for item-specific alpha initial values or
#'   `"fix"` for a common alpha initial value.
#' @param distance_eps Small positive constant added to the distance.
#' @param min_sigma Minimum initial residual standard deviation.
#'
#' @return A named list of initial values.
#' @export
make_init_mdu <- function(Y, D, distance = c("squared", "euclidean"),
                          alpha = c("random", "fix"),
                          distance_eps = 1e-4, min_sigma = 0.1) {
  distance <- match.arg(distance)
  alpha <- match.arg(alpha)
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

  init_alpha <- if (alpha == "random") {
    a <- colMeans(Y, na.rm = TRUE)
    a[!is.finite(a)] <- mean(Y, na.rm = TRUE)
    as.numeric(a)
  } else {
    mean(Y, na.rm = TRUE)
  }

  D_sq <- matrix(NA_real_, N, M)
  for (i in seq_len(N)) {
    for (j in seq_len(M)) {
      D_sq[i, j] <- sum((init_theta[i, ] - init_delta[j, ])^2)
    }
  }
  D_used <- if (distance == "squared") D_sq + distance_eps else sqrt(D_sq + distance_eps)

  alpha_for_mu <- if (length(init_alpha) == M) {
    matrix(init_alpha, N, M, byrow = TRUE)
  } else {
    matrix(init_alpha, N, M)
  }

  utility_init <- alpha_for_mu - D_used
  init_beta <- sum(Y * utility_init, na.rm = TRUE) / sum(utility_init^2, na.rm = TRUE)
  if (!is.finite(init_beta) || init_beta <= 0) init_beta <- 1
  init_beta <- max(init_beta, 1e-4)

  mu <- init_beta * utility_init
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

.validate_mds_distance_matrix <- function(data) {
  Y <- as.matrix(data)
  storage.mode(Y) <- "double"

  if (!is.matrix(Y) || nrow(Y) != ncol(Y)) {
    stop("For method = 'MDS', 'data' must be a square distance matrix.", call. = FALSE)
  }
  if (anyNA(Y)) {
    stop("For method = 'MDS', 'data' must not contain missing values.", call. = FALSE)
  }
  if (any(!is.finite(Y))) {
    stop("For method = 'MDS', 'data' must contain only finite values.", call. = FALSE)
  }
  if (any(Y < 0)) {
    stop("For method = 'MDS', 'data' must be a non-negative distance matrix.", call. = FALSE)
  }
  if (!isTRUE(all.equal(Y, t(Y), tolerance = 1e-8, check.attributes = FALSE))) {
    stop("For method = 'MDS', 'data' must be a symmetric distance matrix.", call. = FALSE)
  }
  if (any(abs(diag(Y)) > 1e-8)) {
    stop("For method = 'MDS', the diagonal of the distance matrix must be zero.", call. = FALSE)
  }
  Y
}

.to_centered_tri_init <- function(coords, D) {
  coords <- as.matrix(coords)
  N <- nrow(coords)
  if (ncol(coords) < D) {
    coords <- cbind(coords, matrix(0, nrow = N, ncol = D - ncol(coords)))
  }
  coords <- coords[, seq_len(D), drop = FALSE]
  coords[!is.finite(coords)] <- 0
  coords <- sweep(coords, 2, colMeans(coords), "-")

  # Rotate the leading D x D block into a Cholesky-like lower-triangular
  # orientation so the constrained initial values match type = "centered_tri".
  if (N >= D) {
    qr_res <- qr(t(coords[seq_len(D), , drop = FALSE]))
    Q <- qr.Q(qr_res)
    R <- qr.R(qr_res)
    sign_diag <- sign(diag(R))
    sign_diag[sign_diag == 0] <- 1
    Q <- sweep(Q, 2, sign_diag, "*")
    coords <- coords %*% Q
  }

  coords_tri <- matrix(0, nrow = N, ncol = D)
  for (d in seq_len(D)) {
    vals <- coords[d:N, d]
    vals <- vals - mean(vals)
    coords_tri[d:N, d] <- vals
  }
  coords_tri
}

.make_init_mds <- function(Y, D, distance = c("euclidean", "squared"), min_sigma = 0.1) {
  distance <- match.arg(distance)
  coords <- tryCatch(
    stats::cmdscale(stats::as.dist(Y), k = D, eig = FALSE),
    error = function(e) NULL
  )
  if (is.null(coords)) {
    coords <- matrix(stats::rnorm(nrow(Y) * D, sd = 0.1), nrow = nrow(Y), ncol = D)
  }
  delta <- .to_centered_tri_init(coords, D)

  pred <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  for (i in seq_len(nrow(Y) - 1L)) {
    for (j in (i + 1L):nrow(Y)) {
      d_sq <- sum((delta[i, ] - delta[j, ])^2)
      pred[i, j] <- if (distance == "squared") d_sq else sqrt(d_sq)
    }
  }
  obs <- Y[upper.tri(Y)]
  fitted <- pred[upper.tri(pred)]
  sigma <- stats::sd(obs - fitted, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) sigma <- stats::sd(obs, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  tau <- stats::sd(as.numeric(delta), na.rm = TRUE)
  if (!is.finite(tau) || tau <= 0) tau <- 1

  list(delta = delta, tau = tau, sigma = max(sigma, min_sigma))
}

#' Make Best and Worst Responses from Best-Worst Pair Indices
#'
#' @description
#' Converts a matrix of Best-Worst pair indices (`Y_dif`) into separate `Best`
#' and `Worst` data frames. The returned values are positions within each row of
#' `sets`, not the item labels stored in `sets`.
#'
#' @param Y_dif Matrix or data frame of pair indices (N persons x P tasks).
#'   Each value must be an integer from 1 to `C * (C - 1)`, where `C` is the
#'   number of items per task.
#' @param sets Matrix or data frame of presented item sets (P tasks x C items).
#'
#' @return A list with two data frames: `Best` and `Worst`.
#' @export
make_bw_from_ydif <- function(Y_dif, sets) {
  Y_dif <- .as_mdu_integer_matrix(Y_dif, "Y_dif")
  sets <- .as_mdu_integer_matrix(sets, "sets")

  N <- nrow(Y_dif)
  P <- ncol(Y_dif)
  C <- ncol(sets)

  if (nrow(sets) != P) {
    stop("nrow(sets) must be equal to ncol(Y_dif).", call. = FALSE)
  }
  if (C < 2) {
    stop("'sets' must have at least two columns.", call. = FALSE)
  }
  if (anyNA(Y_dif)) {
    stop("'Y_dif' cannot contain missing values.", call. = FALSE)
  }

  idx_i <- rep(seq_len(C), each = C)
  idx_j <- rep(seq_len(C), times = C)
  mask <- idx_i != idx_j
  best_rel <- idx_i[mask]
  worst_rel <- idx_j[mask]

  Y_dif_int <- matrix(as.integer(Y_dif), nrow = N, ncol = P)
  max_pair <- C * (C - 1)
  if (any(Y_dif_int < 1 | Y_dif_int > max_pair)) {
    stop("'Y_dif' must contain values from 1 to C * (C - 1).", call. = FALSE)
  }

  Best <- matrix(NA, nrow = N, ncol = P)
  Worst <- matrix(NA, nrow = N, ncol = P)
  for (n in seq_len(N)) {
    for (p in seq_len(P)) {
      ans <- Y_dif_int[n, p]
      Best[n, p] <- best_rel[ans]
      Worst[n, p] <- worst_rel[ans]
    }
  }

  colnames(Best) <- colnames(Y_dif)
  colnames(Worst) <- colnames(Y_dif)
  rownames(Best) <- rownames(Y_dif)
  rownames(Worst) <- rownames(Y_dif)

  list(
    Best = as.data.frame(Best),
    Worst = as.data.frame(Worst)
  )
}

#' @rdname make_bw_from_ydif
#' @export
restore_bw_from_ydif <- make_bw_from_ydif

#' Make Best-Worst Pair Indices from Best and Worst Responses
#'
#' @description
#' Converts separate `Best` and `Worst` response matrices into `Y_dif` pair
#' indices for Best-Worst MDU models. The returned index is the position of the
#' ordered `(best, worst)` pair among all `C * (C - 1)` position pairs generated
#' from each row of `sets`.
#'
#' @param Best Matrix or data frame of best responses (N persons x P tasks).
#'   Values must be positions within the corresponding row of `sets`, from
#'   `1` to `ncol(sets)`.
#' @param Worst Matrix or data frame of worst responses (N persons x P tasks).
#'   Values must be positions within the corresponding row of `sets`, from
#'   `1` to `ncol(sets)`.
#' @param sets Matrix or data frame of presented item sets (P tasks x C items).
#'
#' @return An integer matrix of pair indices (`Y_dif`).
#' @export
make_ydif_from_bw <- function(Best, Worst, sets) {
  Best <- .as_mdu_integer_matrix(Best, "Best")
  Worst <- .as_mdu_integer_matrix(Worst, "Worst")
  sets <- .as_mdu_integer_matrix(sets, "sets")

  if (!identical(dim(Best), dim(Worst))) {
    stop("'Best' and 'Worst' must have the same dimensions.", call. = FALSE)
  }
  if (anyNA(Best) || anyNA(Worst)) {
    stop("'Best' and 'Worst' cannot contain missing values.", call. = FALSE)
  }

  N <- nrow(Best)
  P <- ncol(Best)
  C <- ncol(sets)
  if (nrow(sets) != P) {
    stop("nrow(sets) must be equal to ncol(Best).", call. = FALSE)
  }
  if (C < 2) {
    stop("'sets' must have at least two columns.", call. = FALSE)
  }

  best_dimnames <- dimnames(Best)
  Best <- matrix(as.integer(Best), nrow = N, ncol = P, dimnames = best_dimnames)
  Worst <- matrix(as.integer(Worst), nrow = N, ncol = P, dimnames = dimnames(Worst))
  if (any(Best < 1 | Best > C)) {
    stop("Each Best response must be a position from 1 to ncol(sets).", call. = FALSE)
  }
  if (any(Worst < 1 | Worst > C)) {
    stop("Each Worst response must be a position from 1 to ncol(sets).", call. = FALSE)
  }

  idx_i <- rep(seq_len(C), each = C)
  idx_j <- rep(seq_len(C), times = C)
  mask <- idx_i != idx_j
  best_rel <- idx_i[mask]
  worst_rel <- idx_j[mask]

  Y_dif <- matrix(NA_integer_, nrow = N, ncol = P)
  for (n in seq_len(N)) {
    for (p in seq_len(P)) {
      b_rel <- Best[n, p]
      w_rel <- Worst[n, p]
      if (b_rel == w_rel) {
        stop("Best and Worst responses cannot be the same item within a task.", call. = FALSE)
      }
      Y_dif[n, p] <- which(best_rel == b_rel & worst_rel == w_rel)
    }
  }

  colnames(Y_dif) <- colnames(Best)
  rownames(Y_dif) <- rownames(Best)
  Y_dif
}

.as_mdu_integer_matrix <- function(x, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) {
    stop("'", name, "' must be a matrix or data frame.", call. = FALSE)
  }
  x
}

.normalize_mdu_choice_data <- function(data, sets, method) {
  if (is.null(sets)) {
    stop("'sets' must be supplied when method is 'Best' or 'Best-Worst'.", call. = FALSE)
  }
  y_dif_input <- method == "Best-Worst" && !is.list(data)
  if (!is.list(data) && !y_dif_input) {
    stop("For choice MDU, 'data' must be a list containing 'Best' and optionally 'Worst'.", call. = FALSE)
  }
  if (!y_dif_input && is.null(data$Best)) {
    stop("'data$Best' is required when method is 'Best' or 'Best-Worst'.", call. = FALSE)
  }
  if (!y_dif_input && method == "Best-Worst" && is.null(data$Worst)) {
    stop("'data$Worst' is required when method is 'Best-Worst'.", call. = FALSE)
  }

  S_raw <- .as_mdu_integer_matrix(sets, "sets")
  if (anyNA(S_raw)) stop("'sets' cannot contain missing values.", call. = FALSE)
  P <- nrow(S_raw)
  C <- ncol(S_raw)
  if (P < 1 || C < 2) {
    stop("'sets' must have at least one row and at least two columns.", call. = FALSE)
  }

  set_values <- as.vector(S_raw)
  use_names <- is.character(S_raw) || is.factor(S_raw)
  if (use_names) {
    item_names <- unique(as.character(set_values))
    S <- matrix(match(as.character(S_raw), item_names), nrow = P, ncol = C)
  } else {
    S <- matrix(as.integer(S_raw), nrow = P, ncol = C)
    if (any(!is.finite(S)) || any(S < 1)) {
      stop("'sets' must contain positive item indices or item names.", call. = FALSE)
    }
    item_names <- paste0("Item", seq_len(max(S)))
  }
  if (any(apply(S, 1, function(z) anyDuplicated(z) > 0L))) {
    stop("Each row of 'sets' must contain distinct items.", call. = FALSE)
  }
  M <- length(item_names)

  idx_i <- rep(seq_len(C), each = C)
  idx_j <- rep(seq_len(C), times = C)
  mask <- idx_i != idx_j
  best_rel_all <- idx_i[mask]
  worst_rel_all <- idx_j[mask]

  if (y_dif_input) {
    Y_dif <- .as_mdu_integer_matrix(data, "data")
    Y_dif <- matrix(as.integer(Y_dif), nrow = nrow(Y_dif), ncol = ncol(Y_dif))
    if (ncol(Y_dif) != P) {
      stop("'data' must have the same number of columns as nrow(sets).", call. = FALSE)
    }
    if (anyNA(Y_dif)) {
      stop("'data' cannot contain missing values when supplied as Y_dif.", call. = FALSE)
    }
    max_pair <- C * (C - 1)
    if (any(Y_dif < 1 | Y_dif > max_pair)) {
      stop(
        "'data' as Y_dif must contain pair indices from 1 to C * (C - 1).",
        call. = FALSE
      )
    }

    N <- nrow(Y_dif)
    Y_best <- matrix(NA_integer_, N, P)
    score_mat <- matrix(0, N, M)
    for (n in seq_len(N)) {
      for (p in seq_len(P)) {
        ans <- Y_dif[n, p]
        b <- S[p, best_rel_all[ans]]
        w <- S[p, worst_rel_all[ans]]
        Y_best[n, p] <- best_rel_all[ans]
        score_mat[n, b] <- score_mat[n, b] + 1
        score_mat[n, w] <- score_mat[n, w] - 1
      }
    }

    return(list(
      Y_best = Y_best,
      Y_dif = Y_dif,
      S = S,
      Best = NULL,
      Worst = NULL,
      score_mat = score_mat,
      N = N,
      M = M,
      P = P,
      C = C,
      item_names = item_names,
      input = "Y_dif"
    ))
  }

  convert_response <- function(x, name, as_position = FALSE) {
    x <- .as_mdu_integer_matrix(x, paste0("data$", name))
    if (ncol(x) != P) {
      stop("'", name, "' must have the same number of columns as nrow(sets).", call. = FALSE)
    }
    if (as_position) {
      if (is.character(x) || is.factor(x)) {
        stop("'", name, "' must contain set positions from 1 to ncol(sets).", call. = FALSE)
      }
      res <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
      if (any(!is.na(res) & (res < 1 | res > C))) {
        stop("'", name, "' must contain set positions from 1 to ncol(sets).", call. = FALSE)
      }
      return(res)
    }
    if (use_names || is.character(x) || is.factor(x)) {
      res <- matrix(match(as.character(x), item_names), nrow = nrow(x), ncol = ncol(x))
    } else {
      res <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
    }
    if (any(!is.na(res) & (res < 1 | res > M))) {
      stop("'", name, "' contains item indices outside the item range in 'sets'.", call. = FALSE)
    }
    res
  }

  position_input <- method == "Best-Worst"
  Best <- convert_response(data$Best, "Best", as_position = position_input)
  if (anyNA(Best)) {
    stop("'data$Best' cannot contain missing values for choice MDU.", call. = FALSE)
  }
  N <- nrow(Best)
  Worst <- NULL
  if (method == "Best-Worst") {
    Worst <- convert_response(data$Worst, "Worst", as_position = TRUE)
    if (anyNA(Worst)) {
      stop("'data$Worst' cannot contain missing values for Best-Worst MDU.", call. = FALSE)
    }
    if (!identical(dim(Best), dim(Worst))) {
      stop("'data$Best' and 'data$Worst' must have the same dimensions.", call. = FALSE)
    }
  }

  Y_best <- matrix(NA_integer_, N, P)
  Y_dif <- if (method == "Best-Worst") make_ydif_from_bw(Best, Worst, S) else matrix(NA_integer_, N, P)
  score_mat <- matrix(0, N, M)

  for (n in seq_len(N)) {
    for (p in seq_len(P)) {
      b_rel <- Best[n, p]
      if (is.na(b_rel)) next
      Y_best[n, p] <- b_rel
      b <- S[p, b_rel]
      score_mat[n, b] <- score_mat[n, b] + 1

      if (method == "Best-Worst") {
        w <- S[p, Worst[n, p]]
        score_mat[n, w] <- score_mat[n, w] - 1
      }
    }
  }

  list(
    Y_best = Y_best,
    Y_dif = Y_dif,
    S = S,
    Best = Best,
    Worst = Worst,
    score_mat = score_mat,
    N = N,
    M = M,
    P = P,
    C = C,
    item_names = item_names,
    input = if (method == "Best-Worst") "Best-Worst" else "Best"
  )
}

.make_init_mdu_choice <- function(score_mat, D, min_lambda = 1e-4) {
  N <- nrow(score_mat)
  M <- ncol(score_mat)
  if (D >= M) {
    stop("'ndim' must be smaller than the number of items.", call. = FALSE)
  }

  alpha_score <- colSums(score_mat, na.rm = TRUE)
  alpha_init <- as.numeric(scale(alpha_score, center = TRUE, scale = FALSE)) * 0.1

  score_centered <- scale(score_mat, center = TRUE, scale = FALSE)
  score_centered <- t(scale(t(score_centered), center = TRUE, scale = FALSE))
  score_centered[!is.finite(score_centered)] <- 0

  svd_res <- svd(score_centered, nu = D, nv = D)
  sv <- pmax(svd_res$d[seq_len(D)], 0)
  theta_temp <- sweep(svd_res$u[, seq_len(D), drop = FALSE], 2, sqrt(sv), "*")
  delta_temp <- sweep(svd_res$v[, seq_len(D), drop = FALSE], 2, sqrt(sv), "*")

  qr_obj <- qr(t(delta_temp))
  Q_rot <- qr.Q(qr_obj)
  R_rot <- qr.R(qr_obj)
  sign_diag <- sign(diag(R_rot))
  sign_diag[sign_diag == 0] <- 1
  Q_rot <- sweep(Q_rot, 2, sign_diag, "*")
  delta_rot <- delta_temp %*% Q_rot
  theta_rot <- theta_temp %*% Q_rot

  delta_sd <- apply(delta_rot, 2, stats::sd, na.rm = TRUE)
  delta_sd[!is.finite(delta_sd) | delta_sd <= 1e-8] <- 1
  delta_scaled <- sweep(delta_rot, 2, delta_sd, "/")
  theta_scaled <- sweep(theta_rot, 2, delta_sd, "/")

  init_delta <- matrix(0, nrow = M, ncol = D)
  for (d in seq_len(D)) {
    vals <- delta_scaled[d:M, d]
    vals <- vals - mean(vals)
    init_delta[d:M, d] <- vals
  }
  init_theta <- sweep(theta_scaled, 2, colMeans(theta_scaled), "-")

  lambda_init <- 1
  if (!is.finite(lambda_init) || lambda_init < min_lambda) lambda_init <- min_lambda

  list(
    alpha_raw = alpha_init,
    delta = init_delta,
    theta = init_theta,
    sigma_alpha = 1,
    lambda = lambda_init
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
#' @param data Numeric matrix or data frame (N rows x M items) for
#'   `method = "rating"`. For choice methods, a list containing `Best` and,
#'   for `method = "Best-Worst"`, `Worst`. For `method = "Best-Worst"`, a
#'   single matrix/data frame is treated as pre-coded `Y_dif` pair indices.
#' @param ndim Number of unfolding dimensions.
#' @param distance Character; `"squared"` uses squared Euclidean distance
#'   (default, often easier for optimization), while `"euclidean"` uses Euclidean
#'   distance.
#' @param alpha Character; `"random"` estimates item-specific alpha values as
#'   random effects (default), while `"fix"` estimates a single common alpha.
#' @param method Character; `"rating"` for continuous ratings, `"Best"` for
#'   best-only choice tasks, `"Best-Worst"` for best-worst choice tasks, or
#'   `"MDS"` for fitting a multidimensional scaling model to a distance matrix.
#' @param sets Matrix or data frame of presented item sets (P tasks x C items)
#'   for choice methods.
#' @param prior Prior configuration: `prior_flat()`, `prior_normal()`, or
#'   `prior_weak()`. `prior_flat()` creates a maximum-likelihood model suitable
#'   for `classic()`. The latent coordinates `delta` and `theta` are always
#'   treated as random effects with normal scale priors, similarly to IRT ability
#'   parameters.
#' @param y_range Optional response range for `method = "rating"`. If supplied
#'   with the default flat prior, `prior_weak()` is used. For choice methods
#'   (`"Best"` and `"Best-Worst"`), `prior_weak()` behaves like
#'   `prior_normal()` because there is no observed rating scale.
#' @param init Optional named list of initial values.
#' @param fixed Optional named list of parameter values to fix.
#' @param view Character vector of parameter names to prioritize in summaries.
#' @param distance_eps Small positive constant added to the distance.
#' @param missing Missing value handling strategy: "listwise" (default) or "fiml" (Full Information Maximum Likelihood).
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#'
#' @return An `RTMB_Model` object.
#' @example inst/examples/ex_mdu.R
#' @export
rtmb_mdu <- function(data, ndim = 2,
                     distance = c("squared", "euclidean"),
                     alpha = c("random", "fix"),
                     method = c("rating", "Best", "Best-Worst", "MDS"),
                     sets = NULL,
                     prior = prior_flat(), y_range = NULL,
                     init = NULL, fixed = NULL, view = NULL,
                     distance_eps = 1e-4,
                     missing = c("listwise", "fiml"), WAIC = FALSE) {

  missing_arg <- match.arg(missing)
  distance_missing <- missing(distance)
  method <- match.arg(method)
  if (method == "MDS" && distance_missing) {
    distance <- "euclidean"
    message("method = 'MDS' uses distance = 'euclidean' by default.")
  } else {
    distance <- match.arg(distance)
  }
  alpha_type <- match.arg(alpha)
  if (isTRUE(WAIC) && method == "rating" && missing_arg != "listwise") {
    stop("WAIC = TRUE is currently supported for rating MDU only with missing = 'listwise'.", call. = FALSE)
  }

  if (!is.numeric(ndim) || length(ndim) != 1L || ndim < 1) {
    stop("'ndim' must be a positive integer.", call. = FALSE)
  }
  D <- as.integer(ndim)

  if (method == "MDS") {
    Y <- .validate_mds_distance_matrix(data)
    N <- nrow(Y)
    M <- N
    if (D >= M) {
      stop("'ndim' must be smaller than the size of the distance matrix.", call. = FALSE)
    }
    item_names <- colnames(Y)
    if (is.null(item_names)) item_names <- rownames(Y)
    if (is.null(item_names)) item_names <- paste0("Item", seq_len(M))
    row_names <- rownames(Y)
    if (!is.null(row_names) && !identical(as.character(row_names), as.character(item_names))) {
      item_names <- row_names
    }
    person_names <- character(0)
  } else if (method == "rating") {
    if (is.data.frame(data)) {
      if (!all(sapply(data, is.numeric))) {
        stop("For rating MDU, all variables in the data must be numeric. Character or factor variables are not supported.", call. = FALSE)
      }
    } else if (!is.numeric(data) && !is.logical(data)) {
      stop("For rating MDU, the data matrix must be numeric.", call. = FALSE)
    }
    Y <- as.matrix(data)
    storage.mode(Y) <- "double"
    if (missing_arg == "listwise") {
      Y <- na.omit(Y)
    } else if (anyNA(Y) && missing_arg != "fiml") {
      stop("rtmb_mdu() requires complete data unless missing = 'fiml'.", call. = FALSE)
    }
    N <- nrow(Y)
    M <- ncol(Y)
    if (D >= M) {
      stop("'ndim' must be smaller than the number of items/columns.", call. = FALSE)
    }
    item_names <- colnames(Y)
    if (is.null(item_names)) item_names <- paste0("Item", seq_len(M))
    person_names <- rownames(Y)
    if (is.null(person_names)) person_names <- paste0("Person", seq_len(N))
  } else {
    choice_data <- .normalize_mdu_choice_data(data, sets, method)
    N <- choice_data$N
    M <- choice_data$M
    if (D >= M) {
      stop("'ndim' must be smaller than the number of items.", call. = FALSE)
    }
    item_names <- choice_data$item_names
    person_source <- if (is.list(data)) data$Best else data
    person_names <- rownames(.as_mdu_integer_matrix(person_source, "data"))
    if (is.null(person_names)) person_names <- paste0("Person", seq_len(N))
  }
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
  if (prior_type == "weak" && method == "rating" && is.null(y_range)) {
    stop("When using prior_weak() with method = 'rating', please specify 'y_range' (e.g., y_range = c(1, 7)).", call. = FALSE)
  }

  mu_alpha_init <- if (method == "rating") mean(Y, na.rm = TRUE) else 0
  alpha_sd <- prior$alpha_sd %||% prior$mu_sd %||% 10
  sigma_alpha_rate <- prior$sigma_alpha_rate %||% (1 / 2.5)
  beta_rate <- prior$beta_rate %||% 0.1
  sigma_rate <- prior$sigma_rate %||% 0.1
  delta_sd <- prior$delta_sd %||% prior$tau_sd %||% 1
  theta_sd <- prior$theta_sd %||% 1
  lambda_rate <- prior$lambda_rate %||% (1 / 2.5)

  if (prior_type == "weak" && method == "rating") {
    y_range <- as.numeric(y_range)
    if (length(y_range) != 2L || any(!is.finite(y_range)) || y_range[1] >= y_range[2]) {
      stop("'y_range' must be an increasing numeric vector of length 2.", call. = FALSE)
    }
    sd_ratio <- prior$sd_ratio %||% 0.5
    max_beta <- prior$max_beta %||% 1
    if (!is.finite(sd_ratio) || sd_ratio <= 0) {
      stop("'prior$sd_ratio' must be a positive finite value.", call. = FALSE)
    }
    if (!is.finite(max_beta) || max_beta <= 0) {
      stop("'prior$max_beta' must be a positive finite value.", call. = FALSE)
    }
    delta_sd <- prior$delta_sd %||% 1
    theta_sd <- prior$theta_sd %||% 1
  }

  if (method == "MDS") {
    sigma_rate_mds <- prior$sigma_rate %||% 1
    tau_rate_mds <- prior$tau_rate %||% prior$lambda_rate %||% 1
    dat_mdu <- list(
      Y = Y,
      N = N,
      ndim = D,
      sigma_rate = sigma_rate_mds,
      tau_rate = tau_rate_mds
    )

    setup_ast <- quote({
      N <- nrow(Y)
      D <- ndim
    })

    param_ast <- if (prior_type == "flat") {
      quote({
        delta <- Dim(c(N, D), type = "centered_tri")
        sigma <- Dim(lower = 0)
      })
    } else {
      quote({
        delta <- Dim(c(N, D), type = "centered_tri", random = TRUE)
        tau <- Dim(lower = 0)
        sigma <- Dim(lower = 0)
      })
    }

    loop_ast <- if (distance == "squared") {
      quote(for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
          d_ij <- squared_distance(delta[i, ], delta[j, ])
          Y[i, j] ~ normal(d_ij, sigma)
        }
      })
    } else {
      quote(for (i in 1:(N - 1)) {
        for (j in (i + 1):N) {
          d_ij <- distance(delta[i, ], delta[j, ])
          Y[i, j] ~ normal(d_ij, sigma)
        }
      })
    }

    model_exprs <- list(loop_ast)
    if (prior_type != "flat") {
      model_exprs[[length(model_exprs) + 1]] <- quote(delta ~ centered_tri_multi_normal(tau))
      model_exprs[[length(model_exprs) + 1]] <- quote(tau ~ exponential(tau_rate))
      model_exprs[[length(model_exprs) + 1]] <- quote(sigma ~ exponential(sigma_rate))
    }
    model_ast <- as.call(c(list(as.name("{")), model_exprs))

    code_obj <- list(setup = setup_ast, parameters = param_ast, model = model_ast)
    if (isTRUE(WAIC)) {
      pair_ll_ast <- if (distance == "squared") {
        quote({
          log_lik <- numeric(N * (N - 1) / 2)
          idx_ll <- 1
          for (i in 1:(N - 1)) {
            for (j in (i + 1):N) {
              d_ij <- squared_distance(delta[i, ], delta[j, ])
              log_lik[idx_ll] <- normal_lpdf(Y[i, j], d_ij, sigma)
              idx_ll <- idx_ll + 1
            }
          }
        })
      } else {
        quote({
          log_lik <- numeric(N * (N - 1) / 2)
          idx_ll <- 1
          for (i in 1:(N - 1)) {
            for (j in (i + 1):N) {
              d_ij <- distance(delta[i, ], delta[j, ])
              log_lik[idx_ll] <- normal_lpdf(Y[i, j], d_ij, sigma)
              idx_ll <- idx_ll + 1
            }
          }
        })
      }
      code_obj$generate <- .rtmb_waic_generate_ast(NULL, pair_ll_ast)
    }
    code_obj$env <- parent.frame()

    if (is.null(init)) {
      init <- .make_init_mds(Y, D, distance = distance)
      if (prior_type == "flat") init$tau <- NULL
    }

    par_names_list <- list(
      delta = list(item_names, dim_names)
    )
    if (is.null(view)) view <- if (prior_type == "flat") c("delta", "sigma") else c("delta", "tau", "sigma")
  } else if (method == "rating") {
    dat_mdu <- list(
      Y = Y,
      ndim = D,
      distance_eps = distance_eps,
      mu_alpha_init = mu_alpha_init,
      alpha_sd = alpha_sd,
      beta_rate = beta_rate,
      sigma_rate = sigma_rate,
      delta_sd = delta_sd,
      theta_sd = theta_sd,
      sigma_alpha_rate = sigma_alpha_rate
    )
    if (prior_type == "weak") {
      dat_mdu$y_range <- y_range
      dat_mdu$sd_ratio <- sd_ratio
      dat_mdu$max_beta <- max_beta
    }

    setup_ast <- if (prior_type == "weak") {
      quote({
        N <- nrow(Y)
        M <- ncol(Y)
        D <- ndim
        half_range <- diff(y_range) / 2
        mu_alpha_init <- mean(y_range)
        alpha_sd <- half_range
        sigma_rate <- 1 / (half_range * sd_ratio)
        beta_rate <- 1 / max_beta
      })
    } else {
      quote({
        N <- nrow(Y)
        M <- ncol(Y)
        D <- ndim
      })
    }

    param_ast <- if (alpha_type == "random") {
      bquote({
        mu_alpha <- Dim()
        alpha_raw <- Dim(M, random = TRUE)
        sigma_alpha <- Dim(1, lower = 0)
        beta <- Dim(lower = 0)
        sigma <- Dim(M, lower = 0)
        delta <- Dim(c(M, D), type = "centered_tri", random = TRUE)
        theta <- Dim(c(N, D), random = TRUE)
      })
    } else {
      bquote({
        alpha <- Dim()
        beta <- Dim(lower = 0)
        sigma <- Dim(M, lower = 0)
        delta <- Dim(c(M, D), type = "centered_tri", random = TRUE)
        theta <- Dim(c(N, D), random = TRUE)
      })
    }

    transform_ast <- if (alpha_type == "random") {
      quote({
        alpha <- mu_alpha + alpha_raw * sigma_alpha
      })
    } else {
      NULL
    }

    loop_ast <- if (distance == "squared" && alpha_type == "random") {
      quote(for (i in 1:N) {
        for (j in 1:M) {
          if (!is.na(Y[i, j])) {
            d_ij <- squared_distance(theta[i, ], delta[j, ]) + distance_eps
            mu_ij <- beta * (alpha[j] - d_ij)
            Y[i, j] ~ normal(mu_ij, sigma[j])
          }
        }
      })
    } else if (distance == "squared") {
      quote(for (i in 1:N) {
        for (j in 1:M) {
          if (!is.na(Y[i, j])) {
            d_ij <- squared_distance(theta[i, ], delta[j, ]) + distance_eps
            mu_ij <- beta * (alpha - d_ij)
            Y[i, j] ~ normal(mu_ij, sigma[j])
          }
        }
      })
    } else if (alpha_type == "random") {
      quote(for (i in 1:N) {
        for (j in 1:M) {
          if (!is.na(Y[i, j])) {
            d_ij <- distance(theta[i, ], delta[j, ]) + distance_eps
            mu_ij <- beta * (alpha[j] - d_ij)
            Y[i, j] ~ normal(mu_ij, sigma[j])
          }
        }
      })
    } else {
      quote(for (i in 1:N) {
        for (j in 1:M) {
          if (!is.na(Y[i, j])) {
            d_ij <- distance(theta[i, ], delta[j, ]) + distance_eps
            mu_ij <- beta * (alpha - d_ij)
            Y[i, j] ~ normal(mu_ij, sigma[j])
          }
        }
      })
    }

    model_exprs <- list(loop_ast)
    model_exprs[[length(model_exprs) + 1]] <- quote(delta ~ centered_tri_multi_normal(delta_sd))
    model_exprs[[length(model_exprs) + 1]] <- quote(for (d in 1:D) theta[, d] ~ normal(0, theta_sd))
    if (alpha_type == "random") {
      model_exprs[[length(model_exprs) + 1]] <- quote(alpha_raw ~ normal(0, 1))
      model_exprs[[length(model_exprs) + 1]] <- quote(sigma_alpha ~ exponential(sigma_alpha_rate))
    }
    if (prior_type %in% c("normal", "weak")) {
      if (alpha_type == "fix") {
        model_exprs[[length(model_exprs) + 1]] <- quote(alpha ~ normal(mu_alpha_init, alpha_sd))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- quote(mu_alpha ~ normal(mu_alpha_init, alpha_sd))
      }
      model_exprs[[length(model_exprs) + 1]] <- quote(beta ~ exponential(beta_rate))
      model_exprs[[length(model_exprs) + 1]] <- quote(sigma ~ exponential(sigma_rate))
    }
    model_ast <- as.call(c(list(as.name("{")), model_exprs))
    code_obj <- list(setup = setup_ast, parameters = param_ast)
    if (!is.null(transform_ast)) code_obj$transform <- transform_ast
    code_obj$model <- model_ast
    if (isTRUE(WAIC)) {
      rating_ll_ast <- if (distance == "squared" && alpha_type == "random") {
        quote({
          log_lik <- numeric(N * M)
          idx_ll <- 1
          for (i in 1:N) for (j in 1:M) {
            d_ij <- squared_distance(theta[i, ], delta[j, ]) + distance_eps
            mu_ij <- beta * (alpha[j] - d_ij)
            log_lik[idx_ll] <- normal_lpdf(Y[i, j], mu_ij, sigma[j])
            idx_ll <- idx_ll + 1
          }
        })
      } else if (distance == "squared") {
        quote({
          log_lik <- numeric(N * M)
          idx_ll <- 1
          for (i in 1:N) for (j in 1:M) {
            d_ij <- squared_distance(theta[i, ], delta[j, ]) + distance_eps
            mu_ij <- beta * (alpha - d_ij)
            log_lik[idx_ll] <- normal_lpdf(Y[i, j], mu_ij, sigma[j])
            idx_ll <- idx_ll + 1
          }
        })
      } else if (alpha_type == "random") {
        quote({
          log_lik <- numeric(N * M)
          idx_ll <- 1
          for (i in 1:N) for (j in 1:M) {
            d_ij <- distance(theta[i, ], delta[j, ]) + distance_eps
            mu_ij <- beta * (alpha[j] - d_ij)
            log_lik[idx_ll] <- normal_lpdf(Y[i, j], mu_ij, sigma[j])
            idx_ll <- idx_ll + 1
          }
        })
      } else {
        quote({
          log_lik <- numeric(N * M)
          idx_ll <- 1
          for (i in 1:N) for (j in 1:M) {
            d_ij <- distance(theta[i, ], delta[j, ]) + distance_eps
            mu_ij <- beta * (alpha - d_ij)
            log_lik[idx_ll] <- normal_lpdf(Y[i, j], mu_ij, sigma[j])
            idx_ll <- idx_ll + 1
          }
        })
      }
      code_obj$generate <- .rtmb_waic_generate_ast(NULL, rating_ll_ast)
    }
    code_obj$env <- parent.frame()

    if (is.null(init)) {
      init <- make_init_mdu(Y, D, distance = distance, alpha = alpha_type, distance_eps = distance_eps)
      if (alpha_type == "random") {
        alpha_init <- init$alpha
        sigma_alpha_init <- stats::sd(alpha_init, na.rm = TRUE)
        if (!is.finite(sigma_alpha_init) || sigma_alpha_init <= 1e-8) sigma_alpha_init <- 1
        init$alpha_raw <- (alpha_init - mean(alpha_init, na.rm = TRUE)) / sigma_alpha_init
        init$mu_alpha <- mean(alpha_init, na.rm = TRUE)
        init$sigma_alpha <- sigma_alpha_init
        init$alpha <- NULL
      }
    }

    par_names_list <- list(
      sigma = item_names,
      delta = list(item_names, dim_names),
      theta = list(person_names, dim_names)
    )
    if (alpha_type == "random") par_names_list$alpha_raw <- item_names

    if (is.null(view)) view <- c("alpha", "beta", "sigma", "delta")
  } else {
    dat_mdu <- list(
      sets = choice_data$S,
      ndim = D,
      distance_eps = distance_eps,
      mu_alpha_init = mu_alpha_init,
      alpha_sd = alpha_sd,
      delta_sd = delta_sd,
      theta_sd = theta_sd,
      sigma_alpha_rate = sigma_alpha_rate,
      lambda_rate = lambda_rate
    )
    if (method == "Best") {
      dat_mdu$Y_best <- choice_data$Y_best
    } else if (identical(choice_data$input, "Y_dif")) {
      dat_mdu$Y_dif <- choice_data$Y_dif
    } else {
      dat_mdu$Best <- choice_data$Best
      dat_mdu$Worst <- choice_data$Worst
    }

    setup_ast <- if (method == "Best") {
      quote({
        S <- sets
        N <- nrow(Y_best)
        P <- nrow(S)
        C <- ncol(S)
        M <- max(S)
        D <- ndim
      })
    } else if (identical(choice_data$input, "Y_dif")) {
      quote({
        S <- sets
        N <- nrow(Y_dif)
        P <- nrow(S)
        C <- ncol(S)
        M <- max(S)
        D <- ndim
      })
    } else {
      quote({
        S <- sets
        N <- nrow(Best)
        P <- nrow(S)
        C <- ncol(S)
        M <- max(S)
        D <- ndim
        Y_dif <- make_ydif_from_bw(Best, Worst, S)
      })
    }

    param_ast <- if (alpha_type == "random") {
      bquote({
        delta <- Dim(c(M, D), type = "centered_tri", random = TRUE)
        theta <- Dim(c(N, D), random = TRUE)
        alpha_raw <- Dim(M, type = "centered", random = TRUE)
        sigma_alpha <- Dim(1, lower = 0)
        lambda <- Dim(1, lower = 0)
      })
    } else {
      bquote({
        delta <- Dim(c(M, D), type = "centered_tri", random = TRUE)
        theta <- Dim(c(N, D), random = TRUE)
        alpha <- Dim()
        lambda <- Dim(1, lower = 0)
      })
    }

    transform_ast <- if (alpha_type == "random") {
      quote({
        alpha <- alpha_raw * sigma_alpha
      })
    } else {
      NULL
    }

    utility_ast <- if (distance == "squared" && alpha_type == "random") {
      quote(for (m in 1:M) {
        U[m] <- alpha[m] - squared_distance(theta[n, ], delta[m, ]) - distance_eps
      })
    } else if (distance == "squared") {
      quote(for (m in 1:M) {
        U[m] <- alpha - squared_distance(theta[n, ], delta[m, ]) - distance_eps
      })
    } else if (alpha_type == "random") {
      quote(for (m in 1:M) {
        U[m] <- alpha[m] - distance(theta[n, ], delta[m, ]) - distance_eps
      })
    } else {
      quote(for (m in 1:M) {
        U[m] <- alpha - distance(theta[n, ], delta[m, ]) - distance_eps
      })
    }

    likelihood_ast <- if (method == "Best") {
      bquote(for (n in 1:N) {
        U <- rep(alpha[1] * 0, M)
        .(utility_ast)
        for (p in 1:P) {
          U_set <- U[S[p, ]]
          Y_best[n, p] ~ categorical_logit(lambda * U_set)
        }
      })
    } else {
      bquote(for (n in 1:N) {
        U <- rep(alpha[1] * 0, M)
        .(utility_ast)
        for (p in 1:P) {
          U_set <- U[S[p, ]]
          Y_dif[n, p] ~ bw_categorical_logit(U_set, lambda)
        }
      })
    }

    model_exprs <- list(likelihood_ast)
    model_exprs[[length(model_exprs) + 1]] <- quote(delta ~ centered_tri_multi_normal(delta_sd))
    model_exprs[[length(model_exprs) + 1]] <- quote(for (d in 1:D) theta[, d] ~ normal(0, theta_sd))
    if (alpha_type == "random") {
      model_exprs[[length(model_exprs) + 1]] <- quote(alpha_raw ~ centered_multi_normal(1))
      model_exprs[[length(model_exprs) + 1]] <- quote(sigma_alpha ~ exponential(sigma_alpha_rate))
    } else if (prior_type %in% c("normal", "weak")) {
      model_exprs[[length(model_exprs) + 1]] <- quote(alpha ~ normal(mu_alpha_init, alpha_sd))
    }
    if (prior_type %in% c("normal", "weak")) {
      model_exprs[[length(model_exprs) + 1]] <- quote(lambda ~ exponential(lambda_rate))
    }
    model_ast <- as.call(c(list(as.name("{")), model_exprs))
    code_obj <- list(setup = setup_ast, parameters = param_ast)
    if (!is.null(transform_ast)) code_obj$transform <- transform_ast
    code_obj$model <- model_ast
    if (isTRUE(WAIC)) {
      choice_ll_ast <- if (method == "Best") {
        bquote({
          log_lik <- numeric(N * P)
          idx_ll <- 1
          for (n in 1:N) {
            U <- rep(alpha[1] * 0, M)
            .(utility_ast)
            for (p in 1:P) {
              U_set <- U[S[p, ]]
              log_lik[idx_ll] <- categorical_logit_lpmf(Y_best[n, p], lambda * U_set)
              idx_ll <- idx_ll + 1
            }
          }
        })
      } else {
        bquote({
          log_lik <- numeric(N * P)
          idx_ll <- 1
          for (n in 1:N) {
            U <- rep(alpha[1] * 0, M)
            .(utility_ast)
            for (p in 1:P) {
              U_set <- U[S[p, ]]
              log_lik[idx_ll] <- bw_categorical_logit_lpmf(Y_dif[n, p], U_set, lambda)
              idx_ll <- idx_ll + 1
            }
          }
        })
      }
      code_obj$generate <- .rtmb_waic_generate_ast(NULL, choice_ll_ast)
    }
    code_obj$env <- parent.frame()

    if (is.null(init)) {
      init <- .make_init_mdu_choice(choice_data$score_mat, D)
      if (alpha_type == "fix") init <- init[setdiff(names(init), c("alpha_raw", "sigma_alpha"))]
      if (alpha_type == "fix") init$alpha <- 0
    }

    par_names_list <- list(
      delta = list(item_names, dim_names),
      theta = list(person_names, dim_names)
    )
    if (alpha_type == "random") par_names_list$alpha_raw <- item_names

    if (is.null(view)) {
      view <- if (alpha_type == "random") {
        c("lambda", "sigma_alpha", "alpha", "delta")
      } else {
        c("lambda", "alpha", "delta")
      }
    }
  }

  if (!is.null(init)) {
    if (!is.null(init$phi_raw) && is.null(init$alpha_raw)) {
      init$alpha_raw <- init$phi_raw
      init$phi_raw <- NULL
    }
    if (!is.null(init$sigma_phi) && is.null(init$sigma_alpha)) {
      init$sigma_alpha <- init$sigma_phi
      init$sigma_phi <- NULL
    }
    if (alpha_type == "fix") {
      init$alpha_raw <- NULL
      init$sigma_alpha <- NULL
    }
  }

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
    method = method,
    alpha = alpha_type,
    input = if (method == "rating") "rating" else if (method == "MDS") "distance" else choice_data$input,
    distance = distance,
    marginal = if (method == "MDS") "delta" else c("delta", "theta"),
    sets = if (method %in% c("rating", "MDS")) NULL else choice_data$S
  )
  obj
}
