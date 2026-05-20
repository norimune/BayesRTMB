#' RTMB-based Factor Analysis Wrapper
#'
#' @description
#' Performs exploratory factor analysis (EFA) using RTMB. It supports standard
#' rotation methods (e.g., varimax, promax) as well as regularized factor analysis
#' using a Spike-and-Slab Prior (SSP) for estimating sparse loading matrices.
#'
#' @param data Observation data frame or matrix (N x J).
#' @param nfactors Number of factors (K).
#' @param rotate String specifying the rotation method (e.g., "varimax", "promax", "ssp"). If NULL, no rotation is applied. Specifying "ssp" performs regularized factor analysis.
#' @param score Logical; if TRUE, factor scores are calculated in the generate block (default is FALSE).
#' @param prior Prior configuration: `prior_flat()`, `prior_normal()`, or `prior_weak()`.
#'   Hyperparameters can be specified within these functions (e.g., `prior_normal(mean_sd = 10, sd_rate = 10)`).
#'   Available parameters for FA: `mean_sd`, `sd_rate`, `loadings_sd`, and `ssp_ratio` (if `rotate = "ssp"`).
#' @param y_range A numeric vector of length 2 specifying the theoretical min and max values of the items.
#'   Used to construct weakly informative priors when `prior = prior_weak()`.
#' @param init List of initial values.
#' @param fixed A named list of parameter values to fix (optional).
#' @param missing Missing value handling strategy: "listwise" (default) or "fiml" (Full Information Maximum Likelihood).
#' @param WAIC Logical; if TRUE, add pointwise `log_lik` to the generate block for WAIC.
#' @example inst/examples/ex_fa.R
#' @export
rtmb_fa <- function(data, nfactors = 1, rotate = NULL, score = FALSE,
                    prior = prior_flat(),
                    y_range = NULL,
                    init = NULL, fixed = NULL,
                    missing = c("listwise", "fiml"), WAIC = FALSE) {
  missing <- match.arg(missing)
  if (isTRUE(WAIC) && missing != "listwise") {
    stop("WAIC = TRUE is currently supported for rtmb_fa() only with missing = 'listwise'.", call. = FALSE)
  }
  if (!is.numeric(nfactors) || length(nfactors) != 1L || is.na(nfactors) ||
      nfactors < 1 || nfactors != as.integer(nfactors)) {
    stop("'nfactors' must be a positive integer.", call. = FALSE)
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
  K <- as.integer(nfactors)
  J <- ncol(Y)
  N <- nrow(Y)

  if (is.null(J) || J < 2L) {
    stop("Factor analysis requires at least two observed variables.", call. = FALSE)
  }
  if (is.null(N) || N < 2L) {
    stop("Factor analysis requires at least two complete observations.", call. = FALSE)
  }
  if (K >= J) stop("The number of factors (K) must be less than the number of observed variables (J).")

  var_names <- colnames(data)
  if (is.null(var_names)) var_names <- paste0("V", 1:J)
  factor_names <- paste0("Factor", 1:K)
  colnames(Y) <- var_names

  # --- 1. Pre-calculate Sufficient Statistics in R (Handling NAs) ---
  if (any(is.na(Y))) {
    Y_bar <- apply(Y, 2, mean, na.rm = TRUE)
    S_Y <- cov(Y, use = "pairwise.complete.obs") * (N - 1)
  } else {
    Y_bar <- apply(Y, 2, mean)
    S_Y <- cov(Y) * (N - 1)
  }

  # --- 2. Prior Handling (Extracting from prior object) ---
  if (is.null(prior)) prior <- prior_flat()

  # Automatically switch to prior_weak() if y_range is provided and prior is default flat
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
  
  # Error handling: if weak is requested but y_range is missing
  if (prior_type == "weak" && is.null(y_range)) {
    stop("When using 'prior_weak()', you must specify 'y_range' (e.g., y_range = c(1, 5)) to define the scaling of the priors.", call. = FALSE)
  }

  # Default/Initial settings
  prior_mean_center <- 0
  prior_mean_sd <- if (prior_type == "normal") prior$mu_sd else prior$mean_sd
  prior_sd_rate <- if (prior_type == "normal") prior$sigma_rate else prior$sd_rate
  prior_loadings_sd <- if (prior_type == "normal") prior$b_sd else prior$loadings_sd
  ssp_ratio <- if (!is.null(prior$ssp_ratio)) prior$ssp_ratio else 0.25

  # Weak Information Prior Logic (GLMER style)
  if (prior_type == "weak") {
    sd_ratio <- if (!is.null(prior$sd_ratio)) prior$sd_ratio else 0.5
    half_d_y <- diff(y_range) / 2
    base_scale <- half_d_y * sd_ratio
    
    prior_mean_center <- mean(y_range)
    prior_mean_sd <- half_d_y
    prior_sd_rate <- 1.0 / base_scale
    prior_loadings_sd <- base_scale
  }

  dat_fa <- list(
    Y = Y,
    nfactors = K,
    prior_mean_center = prior_mean_center,
    prior_mean_sd = prior_mean_sd,
    prior_sd_rate = prior_sd_rate,
    prior_loadings_sd = prior_loadings_sd
  )

  # Determine if SSP model is used
  is_ssp <- !is.null(rotate) && rotate == "ssp"

  # --- 3. Simplified Setup AST ---
  if (missing == "listwise") {
    setup_ast <- quote({
      N <- nrow(Y)
      J <- ncol(Y)
      K <- nfactors
      Y_bar <- colMeans(Y)
      S_Y <- cov(Y) * (N - 1)
    })
  } else {
    setup_ast <- quote({
      N <- nrow(Y)
      J <- ncol(Y)
      K <- nfactors
    })
  }

  if (is_ssp) {
    dat_fa$ssp_ratio <- ssp_ratio
    if (score || isTRUE(WAIC)) dat_fa$Y <- Y

    # Automatic generation of initial values for SSP
    if (is.null(init)) {
      tryCatch({
        eig <- eigen(S_Y / (N - 1))
        L_pca <- eig$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(eig$values[1:K], 0.01)), nrow = K, ncol = K)
        var_Y <- diag(S_Y / (N - 1))

        if (K > 1) {
          pm <- stats::promax(L_pca, m = 4)
          init_Lambda <- unclass(pm$loadings)
          T_mat <- pm$rotmat
          Phi <- cov2cor(solve(t(T_mat) %*% T_mat))
          init_CF_Omega <- t(chol(Phi))
        } else {
          init_Lambda <- L_pca
          init_CF_Omega <- matrix(1, nrow = 1, ncol = 1)
          Phi <- matrix(1, nrow = 1, ncol = 1)
        }
        h2_raw <- rowSums(init_Lambda * (init_Lambda %*% Phi))
        init_sd <- sqrt(pmax(var_Y - h2_raw, 0.01 * var_Y))
        L_std <- init_Lambda / sqrt(var_Y)
        init_r <- ifelse(abs(L_std) > 0.2, 0.9, 0.1)

        init <- list(mean = Y_bar, Lambda_star = init_Lambda, sd = init_sd, CF_Omega = init_CF_Omega, r = init_r, tau = matrix(1.0, nrow = J, ncol = K))
      }, error = function(e) { init <- NULL })
    }

    param_ast <- quote({
      mean = Dim(dim = J)
      Lambda_star = Dim(c(J, K))
      r = Dim(c(J, K), lower = 0.001, upper = 0.999)
      tau = Dim(c(J, K), lower = 0)
      sd = Dim(J, lower = 0)
      CF_Omega = Dim(c(K, K), type = "CF_corr")
    })

    tran_ast <- quote({
      Lambda <- Lambda_star * r * tau
      h2 <- rowSums(Lambda * (Lambda %*% CF_Omega))
      var_Y <- h2 + sd^2
      sd_Y <- sqrt(var_Y)
      L <- Lambda / sd_Y
      fa_cor <- CF_Omega %*% t(CF_Omega)
    })

    model_exprs <- list()
    # In SSP case, Lambda is a transformed variable. 
    # S_Y ~ sufficient_multi_normal_fa(N, Y_bar, mean, sd, Lambda %*% CF_Omega)
    # Note: Lambda here is the loadings *before* correlation matrix rotation.
    if (missing == "listwise") {
      model_exprs[[1]] <- quote(S_Y ~ sufficient_multi_normal_fa(N, Y_bar, mean, sd, Lambda %*% CF_Omega))
    } else {
      model_exprs[[1]] <- quote(Y ~ fiml_multi_normal_fa(mean, Lambda %*% CF_Omega, sd))
    }
    if (prior_type == "normal") {
       if (!is.null(prior_mean_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(0, .(prior_mean_sd)))
       if (!is.null(prior_sd_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(.(prior_sd_rate)))
    } else if (prior_type == "weak") {
       if (!is.null(prior_mean_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(prior_mean_center, prior_mean_sd))
       if (!is.null(prior_sd_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(prior_sd_rate))
    }
    model_exprs[[length(model_exprs) + 1]] <- quote(CF_Omega ~ lkj_CF_corr(1))
    model_exprs[[length(model_exprs) + 1]] <- bquote(mu_r <- log(ssp_ratio / (1 - ssp_ratio)))
    model_exprs[[length(model_exprs) + 1]] <- quote(r ~ logit_normal(mu_r, 3))
    model_exprs[[length(model_exprs) + 1]] <- quote(tau ~ exponential(1))
    model_exprs[[length(model_exprs) + 1]] <- quote(Lambda_star ~ laplace(0, 1))
    model_ast <- as.call(c(list(as.name("{")), model_exprs))

    base_gq <- quote({
      Sigma <- Lambda %*% fa_cor %*% t(Lambda) + diag(sd^2)
      var_total <- diag(Sigma)
      var_common <- rowSums(Lambda * (Lambda %*% CF_Omega))
      communality <- var_common / var_total
      out <- list(communality = communality)
    })
    waic_expr <- if (isTRUE(WAIC)) {
      quote({ out$log_lik <- multi_normal_lpdf(Y, mean, Sigma, sum = FALSE) })
    } else {
      quote({})
    }
    score_expr <- if (score) {
      quote({
        # Use explicit matrix creation to ensure AD type is preserved
        Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE)
        out$score <- Y_c %*% solve(Sigma, Lambda %*% fa_cor)
      })
    } else quote({})
    ret_expr <- quote({ return(out) })
    gq_ast <- as.call(c(list(as.name("{")), as.list(base_gq)[-1], as.list(waic_expr)[-1], as.list(score_expr)[-1], as.list(ret_expr)[-1]))
    
    code_obj <- list(setup = setup_ast, parameters = param_ast, transform = tran_ast, model = model_ast, generate = gq_ast, env = parent.frame())
    p_names <- list(
      mean = var_names,
      Lambda_star = list(var_names, factor_names),
      r = list(var_names, factor_names),
      tau = list(var_names, factor_names),
      sd = var_names,
      CF_Omega = factor_names,
      Lambda = list(var_names, factor_names),
      L = list(var_names, factor_names),
      fa_cor = factor_names,
      communality = var_names
    )
    if (score) {
      ind_names <- rownames(data); if (is.null(ind_names)) ind_names <- paste0("Id", 1:N)
      p_names[["score"]] <- list(ind_names, factor_names)
    }
    obj <- rtmb_model(data = dat_fa, code = code_obj, par_names = p_names, init = init, view = c("L", "sd", "fa_cor"), fixed = fixed)
    obj$type <- "fa"
    obj$extra <- list(
      source = "wrapper",
      prior_type = if (is_ssp) "ssp" else prior$type,
      marginal = "mean"
    )
    return(obj)

  } else {
    # --- Standard rotation logic ---
    if (score || isTRUE(WAIC)) dat_fa$Y <- Y
    if (is.null(init)) {
      tryCatch({
        eig <- eigen(S_Y / (N - 1))
        L_pca <- eig$vectors[, 1:K, drop = FALSE] %*% diag(sqrt(pmax(eig$values[1:K], 0.01)), nrow=K, ncol=K)
        L_top <- L_pca[1:K, , drop = FALSE]; qr_res <- qr(t(L_top)); Q_mat <- qr.Q(qr_res)
        init_loadings <- L_pca %*% Q_mat
        init_sd <- pmax(diag(S_Y / (N - 1)) - rowSums(init_loadings^2), 0.01)^0.5
        for (j in 1:J) for (k in 1:K) if (j < k) init_loadings[j, k] <- 0
        init <- list(mean = Y_bar, L_raw = init_loadings, sd = init_sd)
      }, error = function(e) { init <- NULL })
    }

    param_ast <- quote({
      mean <- Dim(dim = J)
      L_raw <- Dim(dim = c(J, K), type = "lower_tri")
      sd <- Dim(dim = J, lower = 0)
    })

    tran_ast <- quote({
      h2 <- rowSums(L_raw^2)
      var_Y <- h2 + sd^2
      sd_Y <- sqrt(var_Y)
      L <- L_raw/sd_Y
    })

    model_exprs <- list()
    # Standard FA case: L_raw is a parameter (lower triangular matrix)
    if (missing == "listwise") {
      model_exprs[[1]] <- quote(S_Y ~ sufficient_multi_normal_fa(N, Y_bar, mean, sd, L_raw))
    } else {
      model_exprs[[1]] <- quote(Y ~ fiml_multi_normal_fa(mean, L_raw, sd))
    }
    if (prior_type == "normal") {
       if (!is.null(prior_mean_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(0, .(prior_mean_sd)))
       if (!is.null(prior_sd_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(.(prior_sd_rate)))
       if (!is.null(prior_loadings_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(L_raw ~ lower_tri_normal(0, .(prior_loadings_sd)))
    } else if (prior_type == "weak") {
       if (!is.null(prior_mean_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(0, prior_mean_sd))
       if (!is.null(prior_sd_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(prior_sd_rate))
       if (!is.null(prior_loadings_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(L_raw ~ lower_tri_normal(0, prior_loadings_sd))
    }
    model_ast <- as.call(c(list(as.name("{")), model_exprs))

    base_gq <- quote({
      Sigma <- L_raw %*% t(L_raw) + diag(sd^2)
      var_total <- diag(Sigma); var_common <- rowSums(L_raw^2); communality <- var_common / var_total
      out <- list(communality = communality)
    })
    waic_expr <- if (isTRUE(WAIC)) {
      quote({ out$log_lik <- multi_normal_lpdf(Y, mean, Sigma, sum = FALSE) })
    } else {
      quote({})
    }

    if (!is.null(rotate)) {
      rot_loadings_name <- paste0("L_", rotate)
      if (exists(rotate, mode = "function")) { rot_fn <- match.fun(rotate); fn_call <- as.name(rotate) }
      else if (requireNamespace("GPArotation", quietly = TRUE) && exists(rotate, where = asNamespace("GPArotation"), mode = "function")) {
        rot_fn <- getFromNamespace(rotate, "GPArotation"); fn_call <- call("::", as.name("GPArotation"), as.name(rotate))
      } else stop("Rotation function not found: ", rotate)
      
      dummy_L <- matrix(rnorm(J * K), J, K); test_rot <- rot_fn(dummy_L)
      is_matrix_rot <- is.matrix(test_rot); has_phi <- !is_matrix_rot && !is.null(test_rot$Phi)
      if (is_matrix_rot) {
        rot_expr <- bquote({
          rot_obj <- .(fn_call)(L)
          out[[.(rot_loadings_name)]] <- unclass(rot_obj)
          dimnames(out[[.(rot_loadings_name)]]) <- list(.(var_names), .(factor_names))
          class(out[[.(rot_loadings_name)]]) <- c("rtmb_estimate_matrix", "matrix", "array")
        })
        score_expr <- if (score) bquote({
          Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE)
          rot_raw <- unclass(.(fn_call)(L_raw)); if (!is.matrix(rot_raw)) rot_raw <- unclass(rot_raw$loadings)
          out$score <- Y_c %*% solve(Sigma, rot_raw)
        }) else quote({})
      } else {
        if (has_phi) {
          rot_expr <- bquote({
            rot_obj <- .(fn_call)(L)
            out$fa_cor <- rot_obj$Phi
            out[[.(rot_loadings_name)]] <- unclass(rot_obj$loadings)
            dimnames(out[[.(rot_loadings_name)]]) <- list(.(var_names), .(factor_names))
            class(out[[.(rot_loadings_name)]]) <- c("rtmb_estimate_matrix", "matrix", "array")
          })
          score_expr <- if (score) bquote({
            Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE)
            rot_raw_obj <- .(fn_call)(L_raw); out$score <- Y_c %*% solve(Sigma, unclass(rot_raw_obj$loadings) %*% rot_raw_obj$Phi)
          }) else quote({})
        } else {
          rot_expr <- bquote({
            rot_obj <- .(fn_call)(L)
            out[[.(rot_loadings_name)]] <- unclass(rot_obj$loadings)
            dimnames(out[[.(rot_loadings_name)]]) <- list(.(var_names), .(factor_names))
            class(out[[.(rot_loadings_name)]]) <- c("rtmb_estimate_matrix", "matrix", "array")
          })
          score_expr <- if (score) bquote({
            Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE)
            out$score <- Y_c %*% solve(Sigma, unclass(.(fn_call)(L_raw)$loadings))
          }) else quote({})
        }
      }
    } else { has_phi <- FALSE; rot_expr <- quote({}); score_expr <- if (score) quote({ Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE); out$score <- Y_c %*% solve(Sigma, L_raw) }) else quote({}) }

    ret_expr <- quote({ return(out) })
    gq_ast <- as.call(c(list(as.name("{")), as.list(base_gq)[-1], as.list(waic_expr)[-1], as.list(rot_expr)[-1], as.list(score_expr)[-1], as.list(ret_expr)[-1]))
    
    code_obj <- list(setup = setup_ast, parameters = param_ast, transform = tran_ast, model = model_ast, generate = gq_ast, env = parent.frame())
    p_names <- list(
      mean = var_names,
      L_raw = list(var_names, factor_names),
      sd = var_names,
      L = list(var_names, factor_names),
      communality = var_names
    )
    if (!is.null(rotate)) {
      p_names[[paste0("L_", rotate)]] <- list(var_names, factor_names)
      if (has_phi) p_names[["fa_cor"]] <- factor_names
    }
    if (score) { ind_names <- rownames(data); if (is.null(ind_names)) ind_names <- paste0("Id", 1:N); p_names[["score"]] <- list(ind_names, factor_names) }

    target_view <- if (!is.null(rotate)) c(paste0("L_", rotate), "sd", "fa_cor") else c("L", "sd", "fa_cor")
    obj <- rtmb_model(data = dat_fa, code = code_obj, par_names = p_names, init = init, view = target_view, fixed = fixed)
    obj$type <- "fa"
    obj$extra <- list(
      source = "wrapper",
      prior_type = prior$type,
      marginal = "mean"
    )
    return(obj)
  }
}
