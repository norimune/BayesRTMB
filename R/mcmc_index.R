# Standalone implementation of rank-normalized R-hat and ESS
# Based on Vehtari et al. (2021) and the monitornew.R implementation.

#' Compute optimal FFT size
#'
#' @importFrom stats median quantile var density
#' @keywords internal
#' @noRd
fft_next_good_size <- function(N) {
  if (N <= 2) return(2)
  while (TRUE) {
    m <- N
    while ((m %% 2) == 0) m <- m / 2
    while ((m %% 3) == 0) m <- m / 3
    while ((m %% 5) == 0) m <- m / 5
    if (m <= 1) return(N)
    N <- N + 1
  }
}

#' @keywords internal
#' @noRd
autocovariance <- function(y) {
  N <- length(y)
  M <- fft_next_good_size(N)
  Mt2 <- 2 * M
  yc <- y - mean(y)
  yc <- c(yc, rep.int(0, Mt2 - N))
  transform <- stats::fft(yc)
  ac <- stats::fft(Conj(transform) * transform, inverse = TRUE)
  ac <- Re(ac)[1:N] / (N^2 * 2)
  return(ac)
}

#' @keywords internal
#' @noRd
backtransform_ranks <- function(r, c = 3/8) {
  S <- length(r)
  (r - c) / (S - 2 * c + 1)
}

#' @keywords internal
#' @noRd
z_scale <- function(x) {
  r <- rank(x, ties.method = "average")
  z <- stats::qnorm(backtransform_ranks(r, c = 3/8))
  z[is.na(x)] <- NA
  if (!is.null(dim(x))) {
    z <- array(z, dim = dim(x), dimnames = dimnames(x))
  }
  return(z)
}

#' @keywords internal
#' @noRd
split_chains <- function(sims) {
  if (is.vector(sims)) {
    dim(sims) <- c(length(sims), 1)
  }
  niter <- dim(sims)[1]
  if (niter < 2) return(sims)
  half <- niter / 2
  cbind(sims[1:floor(half), , drop = FALSE], sims[ceiling(half + 1):niter, , drop = FALSE])
}

#' @keywords internal
#' @noRd
is_constant <- function(x, tol = .Machine$double.eps) {
  if (all(is.na(x))) return(TRUE)
  abs(max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) < tol
}

#' @keywords internal
#' @noRd
rhat_rfun <- function(sims) {
  if (anyNA(sims)) return(NA)
  if (any(!is.finite(sims))) return(NaN)
  if (is_constant(sims)) return(1.0)
  if (is.vector(sims)) dim(sims) <- c(length(sims), 1)
  
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  if (chains < 2) return(NA)
  
  chain_mean <- colMeans(sims)
  chain_var <- apply(sims, 2, var)
  
  var_between <- n_samples * var(chain_mean)
  var_within <- mean(chain_var)
  
  if (var_within < .Machine$double.eps) return(1.0)
  sqrt((var_between / var_within + n_samples - 1) / n_samples)
}

#' @keywords internal
#' @noRd
ess_rfun <- function(sims) {
  if (is.vector(sims)) dim(sims) <- c(length(sims), 1)
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  if (n_samples < 3L || anyNA(sims)) return(NA)
  if (any(!is.finite(sims))) return(NaN)
  if (is_constant(sims)) return(NA)
  
  acov <- lapply(seq_len(chains), function(i) autocovariance(sims[, i]))
  acov <- do.call(cbind, acov) * n_samples / (n_samples - 1)
  chain_mean <- colMeans(sims)
  mean_var <- mean(acov[1, ])
  var_plus <- mean_var * (n_samples - 1) / n_samples
  if (chains > 1) var_plus <- var_plus + var(chain_mean)
  
  rho_hat_t <- rep.int(0, n_samples)
  t <- 0
  rho_hat_even <- 1
  rho_hat_t[t + 1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
  rho_hat_t[t + 2] <- rho_hat_odd
  
  while (t < nrow(acov) - 5 && !is.nan(rho_hat_even + rho_hat_odd) && (rho_hat_even + rho_hat_odd > 0)) {
    t <- t + 2
    rho_hat_even <- 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
    rho_hat_odd <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat_t[t + 1] <- rho_hat_even
      rho_hat_t[t + 2] <- rho_hat_odd
    } else {
      break
    }
  }
  max_t <- t
  if (rho_hat_even > 0) rho_hat_t[max_t + 1] <- rho_hat_even
  
  t <- 0
  while (t <= max_t - 4) {
    t <- t + 2
    if (rho_hat_t[t + 1] + rho_hat_t[t + 2] > rho_hat_t[t - 1] + rho_hat_t[t]) {
      rho_hat_t[t + 1] <- (rho_hat_t[t - 1] + rho_hat_t[t]) / 2
      rho_hat_t[t + 2] <- rho_hat_t[t + 1]
    }
  }
  
  S <- chains * n_samples
  tau_hat <- -1 + 2 * sum(rho_hat_t[1:max_t]) + rho_hat_t[max_t + 1]
  tau_hat <- max(tau_hat, 1 / log10(S))
  return(S / tau_hat)
}

# --- Public Interface Functions ---

#' Calculate Rank-normalized Split-R-hat
#' @param sims A matrix of samples (iterations x chains).
#' @return A numeric value.
#' @export
r_hat <- function(sims) {
  if (is_constant(sims)) return(1.0)
  bulk_rhat <- rhat_rfun(z_scale(split_chains(sims)))
  sims_folded <- abs(sims - median(sims))
  tail_rhat <- rhat_rfun(z_scale(split_chains(sims_folded)))
  max(bulk_rhat, tail_rhat, na.rm = TRUE)
}

#' Calculate Bulk Effective Sample Size
#' @param sims A matrix of samples (iterations x chains).
#' @return A numeric value.
#' @export
ess_bulk <- function(sims) {
  ess_rfun(z_scale(split_chains(sims)))
}

#' Calculate Tail Effective Sample Size (at 2.5\% and 97.5\% quantiles)
#' @param sims A matrix of samples (iterations x chains).
#' @return A numeric value.
#' @export
ess_tail95 <- function(sims) {
  q025 <- (sims <= quantile(sims, 0.025)) + 0
  q025_ess <- ess_rfun(split_chains(q025))
  q975 <- (sims <= quantile(sims, 0.975)) + 0
  q975_ess <- ess_rfun(split_chains(q975))
  min(q025_ess, q975_ess, na.rm = TRUE)
}

#' Basic Effective Sample Size for a single chain or pooled chains
#' @param sims A numeric vector or matrix of samples.
#' @return A numeric value.
#' @export
ess_basic <- function(sims) {
  ess_rfun(split_chains(sims))
}

#' Maximum A Posteriori (MAP) Estimate
#' @param z A numeric vector of samples.
#' @return A numeric value.
#' @export
map_est <- function(z) {
  if (abs(max(z) - min(z)) < 1e-10) return(z[1])
  d <- density(z)
  d$x[which.max(d$y)]
}

#' Calculate 95\% Quantiles
#' @param x A numeric vector.
#' @return A numeric vector of length 2.
#' @export
quantile95 <- function(x) {
  quantile(x, probs = c(0.025, 0.975), names = TRUE)
}
