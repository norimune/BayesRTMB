#' Base class for RTMB Fit objects
#'
#' An R6 base class providing common methods for Bayesian and MAP inference results.
#'
#' @field model The `RTMB_Model` object used for estimation.
#'
#' @importFrom R6 R6Class
#' @export
RTMB_Fit_Base <- R6::R6Class(
  classname = "RTMB_Fit_Base",
  
  public = list(
    model = NULL,

    #' @description Abstract method to get a point estimate for a target parameter.
    #' @param target Character string specifying the target parameter name.
    #' @return A numeric array or matrix of the point estimate.
    get_point_estimate = function(target) {
      stop("get_point_estimate must be implemented by subclasses.")
    },

    #' @description Calculate Expected A Posteriori (EAP) estimates from posterior samples.
    #' @param pars Character vector specifying the names of parameters to extract.
    #'        Use "parameters" for only model parameters, "all" for all variables
    #'        including transformed and generated quantities, or a character vector
    #'        of specific variable names. Default is "parameters".
    #' @param chains Numeric vector specifying the chains to use. Default is NULL (all chains).
    #' @param best_chains Integer; number of best chains to retain based on mean log-posterior (lp) or ELBO. Default is NULL.
    #' @return A named list of EAP estimates structured for use as `init`.
    EAP = function(pars = "parameters", chains = NULL, best_chains = NULL) {
      inc_tran <- TRUE
      inc_gen <- TRUE

      if (identical(pars, "parameters")) {
        target_vars <- names(self$model$par_list)
        inc_tran <- FALSE
        inc_gen <- FALSE
      } else if (identical(pars, "all")) {
        target_vars <- c(names(self$model$par_list),
                         if(!is.null(self$transform_dims)) names(self$transform_dims) else NULL,
                         if(!is.null(self$generate_dims)) names(self$generate_dims) else NULL)
      } else {
        target_vars <- pars
        inc_tran <- any(target_vars %in% names(self$transform_dims))
        inc_gen <- any(target_vars %in% names(self$generate_dims))
      }

      samps <- self$draws(pars = target_vars, chains = chains, best_chains = best_chains,
                          inc_random = TRUE, inc_transform = inc_tran, inc_generate = inc_gen)

      if (dim(samps)[3] == 0) stop("No matching parameters found.")

      flat_names <- dimnames(samps)[[3]]
      if (any(duplicated(flat_names))) {
        # 同じ名前が複数ある場合は最後（最新）のインデックスを残す
        keep_idx <- !rev(duplicated(rev(flat_names)))
        samps <- samps[, , keep_idx, drop = FALSE]
        flat_names <- dimnames(samps)[[3]]
      }

      # Calculate mean over iterations and chains
      eap_flat <- apply(samps, 3, mean, na.rm = TRUE)

      res <- list()
      for (v in target_vars) {
        if (v %in% names(self$model$par_list)) {
          v_dim <- self$model$par_list[[v]]$dim
        } else if (!is.null(self$transform_dims) && v %in% names(self$transform_dims)) {
          v_dim <- self$transform_dims[[v]]
        } else if (!is.null(self$generate_dims) && v %in% names(self$generate_dims)) {
          v_dim <- self$generate_dims[[v]]
        } else {
          next
        }

        pattern <- paste0("^", v, "(\\[.*\\])?$")
        match_idx <- grep(pattern, flat_names)

        if (length(match_idx) > 0) {
          val <- unname(eap_flat[match_idx])
          if (length(v_dim) > 1) {
            dim(val) <- v_dim
          }
          res[[v]] <- val
        }
      }
      return(res)
    },

    #' @description Calculate Maximum A Posteriori (MAP) estimates from posterior samples.
    #' @param pars Character vector specifying the names of parameters to extract.
    #'        Use "parameters" for only model parameters, "all" for all variables
    #'        including transformed and generated quantities, or a character vector
    #'        of specific variable names. Default is "parameters".
    #' @param chains Numeric vector specifying the chains to use. Default is NULL (all chains).
    #' @param best_chains Integer; number of best chains to retain based on mean log-posterior (lp) or ELBO. Default is NULL.
    #' @param type Character string specifying the type of MAP estimate.
    #'        "marginal" (default) calculates the peak of the marginal posterior density for each parameter.
    #'        "joint" returns the parameter values from the iteration with the highest log-posterior.
    #' @return A named list of MAP estimates structured for use as `init`.
    MAP = function(pars = "parameters", chains = NULL, best_chains = NULL, type = c("marginal", "joint")) {
      type <- match.arg(type)

      inc_tran <- TRUE
      inc_gen <- TRUE

      if (identical(pars, "parameters")) {
        target_vars <- names(self$model$par_list)
        inc_tran <- FALSE
        inc_gen <- FALSE
      } else if (identical(pars, "all")) {
        target_vars <- c(names(self$model$par_list),
                         if(!is.null(self$transform_dims)) names(self$transform_dims) else NULL,
                         if(!is.null(self$generate_dims)) names(self$generate_dims) else NULL)
      } else {
        target_vars <- pars
        inc_tran <- any(target_vars %in% names(self$transform_dims))
        inc_gen <- any(target_vars %in% names(self$generate_dims))
      }

      samps <- self$draws(pars = target_vars, chains = chains, best_chains = best_chains,
                          inc_random = TRUE, inc_transform = inc_tran, inc_generate = inc_gen)

      if (dim(samps)[3] == 0) stop("No matching parameters found.")

      flat_names <- dimnames(samps)[[3]]
      if (any(duplicated(flat_names))) {
        keep_idx <- !rev(duplicated(rev(flat_names)))
        samps <- samps[, , keep_idx, drop = FALSE]
        flat_names <- dimnames(samps)[[3]]
      }

      if (type == "joint") {
        # 同時MAP推定 (lpが最大のサンプルを採用)
        lp_samps <- self$draws(pars = "lp", chains = chains, best_chains = best_chains,
                               inc_random = FALSE, inc_transform = FALSE, inc_generate = FALSE)
        if (dim(lp_samps)[3] == 0) stop("Log-probability ('lp') not found. Cannot determine joint MAP.")

        max_idx <- which(lp_samps == max(lp_samps, na.rm = TRUE), arr.ind = TRUE)
        best_i <- max_idx[1, 1]
        best_c <- max_idx[1, 2]

        map_flat <- samps[best_i, best_c, ]

      } else {
        # 周辺MAP推定 (各パラメータの分布の頂点を採用: デフォルト)
        map_flat <- apply(samps, 3, function(z) {
          valid_z <- z[is.finite(z)]
          if (length(valid_z) == 0) return(NA)
          if (abs(max(valid_z) - min(valid_z)) < 1e-10) return(valid_z[1])
          d <- density(valid_z)
          return(d$x[which.max(d$y)])
        })
      }

      res <- list()
      for (v in target_vars) {
        if (v %in% names(self$model$par_list)) {
          v_dim <- self$model$par_list[[v]]$dim
        } else if (!is.null(self$transform_dims) && v %in% names(self$transform_dims)) {
          v_dim <- self$transform_dims[[v]]
        } else if (!is.null(self$generate_dims) && v %in% names(self$generate_dims)) {
          v_dim <- self$generate_dims[[v]]
        } else {
          next
        }

        pattern <- paste0("^", v, "(\\[.*\\])?$")
        match_idx <- grep(pattern, flat_names)

        if (length(match_idx) > 0) {
          val <- unname(map_flat[match_idx])
          if (length(v_dim) > 1) {
            dim(val) <- v_dim
          }
          res[[v]] <- val
        }
      }
      return(res)
    },

    #' @description Rotate sampled parameters.
    #' @param target Character string specifying the target variable to base the rotation on.
    #' @param reference Matrix to rotate towards. If NULL, the target's point estimate is used.
    #' @param linked Character vector of variable names to be rotated in the same direction.
    #' @param overwrite Logical; whether to overwrite the stored draws. If FALSE, adds to generated quantities. Default is FALSE.
    #' @param suffix Character string to append to the rotated variable names when overwrite is FALSE. Default is "rot".
    #' @param ... Additional arguments passed to the rotation function.
    #' @return The updated object invisibly.
    rotate = function(target, reference = NULL, linked = NULL) {
      cat("Applying orthogonal Procrustes rotation (Saving to generate as _rot)...\n")

      target_map <- self$get_point_estimate(target)

      ref_name <- paste0(target, "_ref")
      if (!is.null(reference)) {
        self$model$data[[ref_name]] <- reference
      } else {
        self$model$data[[ref_name]] <- target_map
      }

      exprs <- list()
      exprs[[length(exprs) + 1]] <- bquote(svd_res <- svd(t(.(as.name(target))) %*% .(as.name(ref_name))))
      exprs[[length(exprs) + 1]] <- quote(Q <- svd_res$u %*% t(svd_res$v))

      target_rot <- paste0(target, "_rot")
      exprs[[length(exprs) + 1]] <- bquote(.(as.name(target_rot)) <- .(as.name(target)) %*% Q)

      if (!is.null(linked)) {
        for (l_var in linked) {
          l_rot <- paste0(l_var, "_rot")
          exprs[[length(exprs) + 1]] <- bquote(.(as.name(l_rot)) <- .(as.name(l_var)) %*% Q)
        }
      }

      ret_list <- list()
      ret_list[[target_rot]] <- as.name(target_rot)
      if (!is.null(linked)) {
        for (l_var in linked) {
          l_rot <- paste0(l_var, "_rot")
          ret_list[[l_rot]] <- as.name(l_rot)
        }
      }
      ret_list_call <- as.call(c(list(as.name("list")), ret_list))
      exprs[[length(exprs) + 1]] <- bquote(return(.(ret_list_call)))

      gen_ast <- as.call(c(list(as.name("{")), exprs))

      self$generated_quantities(gen_ast)

      # MCMC and VB need to rename dimensions in generate_fit. MAP_Fit does not have generate_fit.
      if (!is.null(self$generate_fit)) {
        all_names <- c(
          if(!is.null(self$fit)) dimnames(self$fit)[[3]] else character(0),
          if(!is.null(self$transform_fit)) dimnames(self$transform_fit)[[3]] else character(0),
          dimnames(self$generate_fit)[[3]]
        )

        get_renamed_dimnames <- function(orig_var, new_var) {
          orig_pattern <- paste0("^", orig_var, "(\\[.*\\])?$")
          orig_flat_names <- grep(orig_pattern, all_names, value = TRUE)
          if (length(orig_flat_names) > 0) {
            return(sub(paste0("^", orig_var), new_var, orig_flat_names))
          }
          return(NULL)
        }

        gen_names <- dimnames(self$generate_fit)[[3]]

        new_target_names <- get_renamed_dimnames(target, target_rot)
        if (!is.null(new_target_names)) {
          target_gen_pattern <- paste0("^", target_rot, "(\\[.*\\])?$")
          idx <- grep(target_gen_pattern, gen_names)
          if (length(idx) == length(new_target_names)) {
            gen_names[idx] <- new_target_names
          }
        }

        if (!is.null(linked)) {
          for (l_var in linked) {
            l_rot <- paste0(l_var, "_rot")
            new_linked_names <- get_renamed_dimnames(l_var, l_rot)
            if (!is.null(new_linked_names)) {
              l_gen_pattern <- paste0("^", l_rot, "(\\[.*\\])?$")
              idx <- grep(l_gen_pattern, gen_names)
              if (length(idx) == length(new_linked_names)) {
                gen_names[idx] <- new_linked_names
              }
            }
          }
        }
        dimnames(self$generate_fit)[[3]] <- gen_names
      }

      invisible(self)
    },

    #' @description Rotate factor loadings and optional factor scores.
    #' @param target Character string specifying the target variable to base the rotation on.
    #' @param linked Character vector of variable names to be rotated in the same direction.
    #' @param scores Character vector of variable names to be rotated as factor scores (inverse direction).
    #' @param rotate Character string specifying the rotation method.
    #' @param ... Additional arguments passed to the rotation function.
    #' @return The updated object invisibly.
    fa_rotate = function(target = "L", linked = NULL, scores = NULL, rotate = "promax", ...) {
      cat(sprintf("Applying %s rotation to %s (Saving to generate as _%s)...\n", rotate, target, rotate))

      if (exists(rotate, mode = "function")) {
        rot_fn <- match.fun(rotate)
        fn_call <- as.name(rotate)
      } else if (requireNamespace("GPArotation", quietly = TRUE) &&
                 exists(rotate, where = asNamespace("GPArotation"), mode = "function")) {
        rot_fn <- getFromNamespace(rotate, "GPArotation")
        fn_call <- call("::", as.name("GPArotation"), as.name(rotate))
      } else {
        stop("Rotation function not found: ", rotate)
      }

      target_map <- self$get_point_estimate(target)
      
      if (length(dim(target_map)) != 2) stop("fa_rotate は行列(2次元)パラメータに対してのみ適用可能です。")

      test_rot <- rot_fn(target_map, ...)
      is_matrix_rot <- is.matrix(test_rot)
      has_phi <- !is_matrix_rot && !is.null(test_rot$Phi)

      exprs <- list()
      rot_name <- paste0(target, "_", rotate)

      exprs[[length(exprs) + 1]] <- bquote(rot_obj <- .(fn_call)(.(as.name(target))))
      ret_list <- list()

      if (is_matrix_rot) {
        exprs[[length(exprs) + 1]] <- quote(rot_mat <- unclass(rot_obj))
        ret_list[[rot_name]] <- as.name("rot_mat")
      } else {
        exprs[[length(exprs) + 1]] <- quote(rot_mat <- unclass(rot_obj$loadings))
        ret_list[[rot_name]] <- as.name("rot_mat")
        if (has_phi) {
          ret_list[["fa_cor"]] <- quote(rot_obj$Phi)
        }
      }

      if (!is.null(linked) || !is.null(scores)) {
        if (is_matrix_rot) {
          warning("指定された回転関数は行列のみを返すため、回転行列が取得できません。linked / scores は回転されません。")
        } else {
          exprs[[length(exprs) + 1]] <- quote(rot_Th <- if (!is.null(rot_obj$Th)) rot_obj$Th else rot_obj$rotmat)

          if (!is.null(linked)) {
            for (l_var in linked) {
              l_rot <- paste0(l_var, "_", rotate)
              exprs[[length(exprs) + 1]] <- bquote(.(as.name(l_rot)) <- .(as.name(l_var)) %*% rot_Th)
              ret_list[[l_rot]] <- as.name(l_rot)
            }
          }

          if (!is.null(scores)) {
            exprs[[length(exprs) + 1]] <- quote(rot_Th_inv <- solve(t(rot_Th)))
            for (s_var in scores) {
              s_rot <- paste0(s_var, "_", rotate)
              exprs[[length(exprs) + 1]] <- bquote(.(as.name(s_rot)) <- .(as.name(s_var)) %*% rot_Th_inv)
              ret_list[[s_rot]] <- as.name(s_rot)
            }
          }
        }
      }

      ret_list_call <- as.call(c(list(as.name("list")), ret_list))
      exprs[[length(exprs) + 1]] <- bquote(return(.(ret_list_call)))

      gen_ast <- as.call(c(list(as.name("{")), exprs))
      self$generated_quantities(gen_ast)

      if (!is.null(self$generate_fit)) {
        all_names <- c(
          if(!is.null(self$fit)) dimnames(self$fit)[[3]] else character(0),
          if(!is.null(self$transform_fit)) dimnames(self$transform_fit)[[3]] else character(0),
          dimnames(self$generate_fit)[[3]]
        )

        get_renamed_dimnames <- function(orig_var, new_var) {
          orig_pattern <- paste0("^", orig_var, "(\\[.*\\])?$")
          orig_flat_names <- grep(orig_pattern, all_names, value = TRUE)
          if (length(orig_flat_names) > 0) {
            return(sub(paste0("^", orig_var), new_var, orig_flat_names))
          }
          return(NULL)
        }

        gen_names <- dimnames(self$generate_fit)[[3]]

        new_target_names <- get_renamed_dimnames(target, rot_name)
        if (!is.null(new_target_names)) {
          target_gen_pattern <- paste0("^", rot_name, "(\\[.*\\])?$")
          idx <- grep(target_gen_pattern, gen_names)
          if (length(idx) == length(new_target_names)) {
            gen_names[idx] <- new_target_names
          }
        }

        if (!is.null(linked)) {
          for (l_var in linked) {
            l_rot <- paste0(l_var, "_", rotate)
            new_linked_names <- get_renamed_dimnames(l_var, l_rot)
            if (!is.null(new_linked_names)) {
              l_gen_pattern <- paste0("^", l_rot, "(\\[.*\\])?$")
              idx <- grep(l_gen_pattern, gen_names)
              if (length(idx) == length(new_linked_names)) {
                gen_names[idx] <- new_linked_names
              }
            }
          }
        }

        if (!is.null(scores)) {
          for (s_var in scores) {
            s_rot <- paste0(s_var, "_", rotate)
            new_score_names <- get_renamed_dimnames(s_var, s_rot)
            if (!is.null(new_score_names)) {
              s_gen_pattern <- paste0("^", s_rot, "(\\[.*\\])?$")
              idx <- grep(s_gen_pattern, gen_names)
              if (length(idx) == length(new_score_names)) {
                gen_names[idx] <- new_score_names
              }
            }
          }
        }

        dimnames(self$generate_fit)[[3]] <- gen_names
      }

      invisible(self)
    }
  )
)
