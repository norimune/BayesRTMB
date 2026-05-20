`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @keywords internal
#' @noRd
.names0 <- function(x) {
  n <- names(x)
  if (is.null(n)) character(0L) else n
}

#' @keywords internal
#' @noRd
.estimate_component_names <- function(fit, component) {
  if (component == "parameters") return(.names0(fit$model$par_list))

  if (component == "transform") {
    n <- .names0(fit$transform_dims)
    if (length(n) == 0L) n <- .names0(fit$transform)
    return(n)
  }

  if (component == "generate") {
    n <- .names0(fit$generate_dims)
    if (length(n) == 0L) n <- .names0(fit$generate)
    return(n)
  }

  stop("Unknown component: ", component, call. = FALSE)
}

#' @keywords internal
#' @noRd
.get_marginal_mode <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  if (length(unique(x)) == 1L) return(x[1L])
  # Avoid density estimation warning when values are nearly constant
  if (diff(range(x)) < 1e-7) return(mean(x))
  d <- density(x)
  d$x[which.max(d$y)]
}

#' @export
print.rtmb_estimate_matrix <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  x_print <- formatC(unclass(x), format = "f", digits = digits)
  dimnames(x_print) <- dimnames(x)
  rn <- rownames(x_print)
  cn <- colnames(x_print)
  if (is.null(rn)) rn <- rep("", nrow(x_print))
  if (is.null(cn)) cn <- paste0("[,", seq_len(ncol(x_print)), "]")

  row_width <- max(nchar(rn), 0L)
  col_widths <- pmax(nchar(cn), apply(nchar(x_print), 2L, max))
  header <- paste0(
    strrep(" ", row_width + 1L),
    paste(format(cn, width = col_widths, justify = "right"), collapse = " ")
  )
  rows <- vapply(seq_len(nrow(x_print)), function(i) {
    paste0(
      format(rn[i], width = row_width, justify = "left"),
      " ",
      paste(format(x_print[i, ], width = col_widths, justify = "right"), collapse = " ")
    )
  }, character(1L))
  cat(paste(c(header, rows), collapse = "\n"), "\n", sep = "")
  invisible(x)
}

#' @keywords internal
#' @noRd
.apply_estimate_dimnames <- function(fit, name, val, dims = NULL) {
  if (is.null(dims)) {
    if (!is.null(fit$model$par_list[[name]])) {
      dims <- fit$model$par_list[[name]]$dim
    } else if (!is.null(fit$transform_dims[[name]])) {
      dims <- fit$transform_dims[[name]]
    } else if (!is.null(fit$generate_dims[[name]])) {
      dims <- fit$generate_dims[[name]]
    }
  }
  if (is.null(dims) || length(dims) == 0L) return(val)

  names_def <- fit$model$par_names[[name]]
  if (length(dims) > 1L && is.null(dim(val))) dim(val) <- dims
  if (is.null(names_def) && length(dims) == 2L) {
    y_names <- colnames(fit$model$data$Y)
    if (!is.null(y_names) && length(y_names) == dims[1L]) {
      names_def <- list(y_names, paste0("Factor", seq_len(dims[2L])))
    }
  }

  if (!is.null(names_def)) {
    if (is.list(names_def) && length(names_def) == length(dims)) {
      dimnames(val) <- names_def
    } else if (is.atomic(names_def) && length(dims) == 1L && length(names_def) == dims[1L]) {
      names(val) <- names_def
    } else if (is.atomic(names_def) && length(dims) > 1L && length(names_def) == dims[1L]) {
      dn <- vector("list", length(dims))
      dn[[1L]] <- names_def
      if (length(dims) == 2L && dims[1L] == dims[2L]) dn[[2L]] <- names_def
      dimnames(val) <- dn
    }
  }
  if (length(dim(val)) == 2L && is.numeric(val)) {
    class(val) <- unique(c("rtmb_estimate_matrix", class(val)))
  }
  val
}

#' @keywords internal
#' @noRd
.collect_point_estimates <- function(fit, type = "EAP", chains = NULL, best_chains = NULL, ...) {
  if (inherits(fit, c("Classic_Fit", "map_fit"))) {
    # Combine par, transform, and generate fields from optimization result
    res <- c(fit$par %||% list(), fit$transform %||% list(), fit$generate %||% list())
    for (name in names(res)) {
      res[[name]] <- .apply_estimate_dimnames(fit, name, res[[name]])
    }
    return(res)
  }

  if (inherits(fit, c("mcmc_fit", "advi_fit"))) {
    # Fetch draws for all variables (respecting chains/best_chains)
    samps <- fit$draws(chains = chains, best_chains = best_chains, 
                       inc_random = TRUE, inc_transform = TRUE, inc_generate = TRUE)
    
    if (type == "joint_map") {
      # Returns the sample with the highest log-posterior (joint MAP)
      # Find the iteration/chain with max lp
      lp_idx <- which(dimnames(samps)[[3]] == "lp")
      if (length(lp_idx) == 0) lp_idx <- 1
      
      lp_vals <- samps[, , lp_idx]
      max_idx <- which(lp_vals == max(lp_vals, na.rm = TRUE), arr.ind = TRUE)
      iter_idx <- max_idx[1, 1]
      chain_idx <- max_idx[1, 2]
      
      # This is our point estimate vector (flattened)
      flat_ests <- samps[iter_idx, chain_idx, ]
      names(flat_ests) <- dimnames(samps)[[3]]
    } else {
      # EAP or marginal_map
      d <- dim(samps)
      samps_flat <- matrix(samps, nrow = d[1] * d[2], ncol = d[3])
      colnames(samps_flat) <- dimnames(samps)[[3]]
      
      if (type == "EAP") {
        flat_ests <- colMeans(samps_flat, na.rm = TRUE)
      } else {
        flat_ests <- apply(samps_flat, 2, .get_marginal_mode)
      }
    }
    
    # Group back into structured list
    res <- list()
    
    # Helper to group by base name
    group_by_name <- function(base_names, dims_list = NULL) {
      for (name in base_names) {
        idx <- grep(paste0("^", name, "(\\[|$)"), names(flat_ests))
        if (length(idx) > 0) {
          val <- flat_ests[idx]
          dims <- if (!is.null(dims_list)) dims_list[[name]] else fit$model$par_list[[name]]$dim
          if (length(dims) > 1) dim(val) <- dims
          val <- .apply_estimate_dimnames(fit, name, val, dims)
          res[[name]] <<- val
        }
      }
    }
    
    group_by_name(names(fit$model$par_list))
    group_by_name(names(fit$transform_dims), fit$transform_dims)
    group_by_name(names(fit$generate_dims), fit$generate_dims)
    
    return(res)
  }
  stop("Unsupported fit object type.")
}

#' @keywords internal
#' @noRd
.select_estimates <- function(fit, est_list, pars = NULL, component = "all", drop = FALSE) {
  component <- match.arg(component, c("all", "parameters", "transform", "generate"))

  # 1. Parameter selection keywords
  if (identical(pars, "all")) {
    pars <- NULL
  } else if (identical(pars, "parameters")) {
    pars <- .estimate_component_names(fit, "parameters")
  } else if (identical(pars, "transform")) {
    pars <- .estimate_component_names(fit, "transform")
  } else if (identical(pars, "generate")) {
    pars <- .estimate_component_names(fit, "generate")
  }

  res <- select_parameters(est_list, pars)

  # 2. Component filtering
  if (component != "all") {
    keep <- .estimate_component_names(fit, component)
    res <- res[names(res) %in% keep]
  }

  if (length(res) == 0L) return(list())
  if (isTRUE(drop) && length(res) == 1L) return(res[[1L]])
  return(res)
}

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
    #' @param ... Additional arguments.
    #' @return A numeric array or matrix of the point estimate.
    get_point_estimate = function(target, ...) {
      stop("get_point_estimate must be implemented by subclasses.")
    },

    #' @description Get point estimates for parameters, transformed parameters, and generated quantities.
    #' @param pars Optional character or numeric vector of parameter names or indices to extract.
    #'        Supports special keywords: "parameters", "transform", "generate", and "all".
    #' @param type Character string specifying the estimation type.
    #' @param component Character string specifying the component to filter by.
    #' @param chains Numeric vector of chains to include.
    #' @param best_chains Number of best chains to include.
    #' @param drop Logical; if TRUE and only one parameter is selected, return the value directly instead of a list.
    #' @param ... Additional arguments passed to draws().
    #' @return A named list of point estimates, or a single value if `drop = TRUE`.
    estimate = function(pars = NULL, 
                        type = c("mean", "EAP", "marginal_map", "joint_map", "MAP"), 
                        component = c("all", "parameters", "transform", "generate"), 
                        chains = NULL, 
                        best_chains = NULL, 
                        drop = TRUE, 
                        ...) {
      names0 <- function(x) {
        n <- names(x)
        if (is.null(n)) character(0L) else n
      }
      coalesce_null <- function(a, b) {
        if (!is.null(a)) a else b
      }
      apply_estimate_dimnames <- function(fit, name, val, dims = NULL) {
        if (is.null(dims)) {
          if (!is.null(fit$model$par_list[[name]])) {
            dims <- fit$model$par_list[[name]]$dim
          } else if (!is.null(fit$transform_dims[[name]])) {
            dims <- fit$transform_dims[[name]]
          } else if (!is.null(fit$generate_dims[[name]])) {
            dims <- fit$generate_dims[[name]]
          }
        }
        if (is.null(dims) || length(dims) == 0L) return(val)

        names_def <- fit$model$par_names[[name]]
        if (length(dims) > 1L && is.null(dim(val))) dim(val) <- dims
        if (is.null(names_def) && length(dims) == 2L) {
          y_names <- colnames(fit$model$data$Y)
          if (!is.null(y_names) && length(y_names) == dims[1L]) {
            names_def <- list(y_names, paste0("Factor", seq_len(dims[2L])))
          }
        }

        if (!is.null(names_def)) {
          if (is.list(names_def) && length(names_def) == length(dims)) {
            dimnames(val) <- names_def
          } else if (is.atomic(names_def) && length(dims) == 1L && length(names_def) == dims[1L]) {
            names(val) <- names_def
          } else if (is.atomic(names_def) && length(dims) > 1L && length(names_def) == dims[1L]) {
            dn <- vector("list", length(dims))
            dn[[1L]] <- names_def
            if (length(dims) == 2L && dims[1L] == dims[2L]) dn[[2L]] <- names_def
            dimnames(val) <- dn
          }
        }
        if (length(dim(val)) == 2L && is.numeric(val)) {
          class(val) <- unique(c("rtmb_estimate_matrix", class(val)))
        }
        val
      }
      select_parameters_local <- function(all_list, pars) {
        if (is.null(pars)) return(all_list)
        if (identical(pars, "parameters") || identical(pars, "all")) return(all_list)

        if (is.numeric(pars)) {
          res <- tryCatch(all_list[pars], error = function(e) NULL)
          if (is.null(res)) return(all_list[0])
          return(res[names(res) != ""])
        }

        if (is.character(pars)) {
          is_negative <- grepl("^-", pars)
          if (any(is_negative)) {
            if (!all(is_negative)) {
              stop("Cannot mix positive and negative parameter names in selection.")
            }
            exclude_names <- gsub("^-", "", pars)
            return(all_list[!(names(all_list) %in% exclude_names)])
          }
          return(all_list[names(all_list) %in% pars])
        }

        return(all_list)
      }
      estimate_component_names <- function(fit, component) {
        if (component == "parameters") return(names0(fit$model$par_list))

        if (component == "transform") {
          n <- names0(fit$transform_dims)
          if (length(n) == 0L) n <- names0(fit$transform)
          return(n)
        }

        if (component == "generate") {
          n <- names0(fit$generate_dims)
          if (length(n) == 0L) n <- names0(fit$generate)
          return(n)
        }

        stop("Unknown component: ", component, call. = FALSE)
      }
      marginal_mode <- function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0L) return(NA_real_)
        if (length(unique(x)) == 1L) return(x[1L])
        # Avoid density estimation warning when values are nearly constant
        if (diff(range(x)) < 1e-7) return(mean(x))
        d <- density(x)
        d$x[which.max(d$y)]
      }
      collect_point_estimates <- function(fit, type = "EAP", chains = NULL, best_chains = NULL, ...) {
        if (inherits(fit, c("Classic_Fit", "map_fit"))) {
          res <- c(coalesce_null(fit$par, list()),
                   coalesce_null(fit$transform, list()),
                   coalesce_null(fit$generate, list()))
          for (name in names(res)) {
            res[[name]] <- apply_estimate_dimnames(fit, name, res[[name]])
          }
          return(res)
        }

        if (inherits(fit, c("mcmc_fit", "advi_fit"))) {
          samps <- fit$draws(chains = chains, best_chains = best_chains,
                             inc_random = TRUE, inc_transform = TRUE, inc_generate = TRUE)

          if (type == "joint_map") {
            lp_idx <- which(dimnames(samps)[[3]] == "lp")
            if (length(lp_idx) == 0) lp_idx <- 1

            lp_vals <- samps[, , lp_idx]
            max_idx <- which(lp_vals == max(lp_vals, na.rm = TRUE), arr.ind = TRUE)
            iter_idx <- max_idx[1, 1]
            chain_idx <- max_idx[1, 2]

            flat_ests <- samps[iter_idx, chain_idx, ]
            names(flat_ests) <- dimnames(samps)[[3]]
          } else {
            d <- dim(samps)
            samps_flat <- matrix(samps, nrow = d[1] * d[2], ncol = d[3])
            colnames(samps_flat) <- dimnames(samps)[[3]]

            if (type == "EAP") {
              flat_ests <- colMeans(samps_flat, na.rm = TRUE)
            } else {
              flat_ests <- apply(samps_flat, 2, marginal_mode)
            }
          }

          res <- list()
          group_by_name <- function(base_names, dims_list = NULL) {
            for (name in base_names) {
              idx <- grep(paste0("^", name, "(\\[|$)"), names(flat_ests))
              if (length(idx) > 0) {
                val <- flat_ests[idx]
                dims <- if (!is.null(dims_list)) dims_list[[name]] else fit$model$par_list[[name]]$dim
                if (length(dims) > 1) dim(val) <- dims
                val <- apply_estimate_dimnames(fit, name, val, dims)
                res[[name]] <<- val
              }
            }
          }

          group_by_name(names(fit$model$par_list))
          group_by_name(names(fit$transform_dims), fit$transform_dims)
          group_by_name(names(fit$generate_dims), fit$generate_dims)

          return(res)
        }
        stop("Unsupported fit object type.")
      }
      select_estimates <- function(fit, est_list, pars = NULL, component = "all", drop = FALSE) {
        component <- match.arg(component, c("all", "parameters", "transform", "generate"))

        if (identical(pars, "all")) {
          pars <- NULL
        } else if (identical(pars, "parameters")) {
          pars <- estimate_component_names(fit, "parameters")
        } else if (identical(pars, "transform")) {
          pars <- estimate_component_names(fit, "transform")
        } else if (identical(pars, "generate")) {
          pars <- estimate_component_names(fit, "generate")
        }

        res <- select_parameters_local(est_list, pars)

        if (component != "all") {
          keep <- estimate_component_names(fit, component)
          res <- res[names(res) %in% keep]
        }

        if (length(res) == 0L) return(list())
        if (isTRUE(drop) && length(res) == 1L) return(res[[1L]])
        return(res)
      }

      type <- match.arg(type)
      component <- match.arg(component)
      
      # Normalization
      if (type == "mean") type <- "EAP"
      if (type == "MAP") type <- "marginal_map"

      # Determine default pars if NULL
      if (is.null(pars)) {
        if (!identical(component, "all")) {
          pars <- component
        } else if (inherits(self, c("mcmc_fit", "advi_fit"))) {
          pars <- "parameters"
        } else {
          pars <- "all"
        }
      }

      est_all <- collect_point_estimates(self, type = type, chains = chains, best_chains = best_chains, ...)
      
      # Validate pars existence if character and not keywords/regex/exclusion
      if (is.character(pars) && !any(grepl("^-", pars)) && 
          !all(pars %in% c("parameters", "transform", "generate", "all"))) {
        missing <- setdiff(pars, names(est_all))
        if (length(missing) > 0L) {
          stop(
            "Unknown estimate variable(s): ", paste(missing, collapse = ", "),
            "\nAvailable variables: ", paste(names(est_all), collapse = ", "),
            call. = FALSE
          )
        }
      }
      
      res <- select_estimates(self, est_all, pars = pars, component = component, drop = drop)
      return(res)
    },

    #' @description Calculate Expected A Posteriori (EAP) estimates from posterior samples.
    #' @param pars Optional character vector of parameter names to extract.
    #' @param chains Numeric vector of chains to include.
    #' @param best_chains Number of best chains to include.
    #' @param drop Logical; whether to drop the list if only one parameter is selected.
    #' @param ... Additional arguments passed to `estimate()`.
    #' @return A named list of EAP estimates.
    EAP = function(pars = "parameters", chains = NULL, best_chains = NULL, drop = FALSE, ...) {
      return(self$estimate(pars = pars, type = "EAP", chains = chains, best_chains = best_chains, drop = drop, ...))
    },

    #' @description Calculate Maximum A Posteriori (MAP) estimates.
    #' @param pars Optional character vector of parameter names to extract.
    #' @param chains Numeric vector of chains to include.
    #' @param best_chains Number of best chains to include.
    #' @param type Character string; "marginal" or "joint" MAP.
    #' @param drop Logical; whether to drop the list if only one parameter is selected.
    #' @param ... Additional arguments passed to `estimate()`.
    #' @return A named list of MAP estimates.
    MAP = function(pars = "parameters", chains = NULL, best_chains = NULL, 
                   type = c("marginal", "joint"), drop = FALSE, ...) {
      type <- match.arg(type)
      est_type <- if (type == "joint") "joint_map" else "marginal_map"
      return(self$estimate(pars = pars, type = est_type, chains = chains, best_chains = best_chains, drop = drop, ...))
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
      message("Applying orthogonal Procrustes rotation (Saving to generate as _rot)...")

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
      message(sprintf("Applying %s rotation to %s (Saving to generate as _%s)...", rotate, target, rotate))

      if (exists(rotate, mode = "function")) {
        rot_fn <- match.fun(rotate)
        fn_call <- as.name(rotate)
      } else if (requireNamespace("GPArotation", quietly = TRUE) &&
                 exists(rotate, where = asNamespace("GPArotation"), mode = "function")) {
        rot_fn <- getFromNamespace(rotate, "GPArotation")
        fn_call <- call("::", as.name("GPArotation"), as.name(rotate))
      } else {
        stop("Rotation function not found: ", rotate, ". If this is from GPArotation, please install it using install.packages('GPArotation').")
      }

      target_map <- self$get_point_estimate(target)

      if (length(dim(target_map)) != 2) stop("fa_rotate is only applicable to matrix (2D) parameters.")

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
          warning("The specified rotation function returns only a matrix, so the rotation matrix cannot be obtained. 'linked' and 'scores' will not be rotated.")
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
