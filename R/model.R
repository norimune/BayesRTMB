#' @import MASS
#' @importFrom grDevices col2rgb hcl.colors rgb
#' @importFrom graphics abline axis legend lines matplot mtext pairs par points polygon segments strwidth text
#' @importFrom stats Gamma binomial coef cor cov delete.response density dunif gaussian glm median model.frame model.matrix model.offset model.response na.omit nlminb poisson quantile rnorm runif sd terms var
#' @importFrom utils capture.output getFromNamespace modifyList read.csv
NULL

.rtmb_closed_eval_env <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    ns <- asNamespace("BayesRTMB")
    out <- new.env(parent = baseenv())
    imports <- parent.env(ns)
    if (startsWith(environmentName(imports), "imports:")) {
      for (nm in ls(imports, all.names = TRUE)) {
        assign(nm, get(nm, envir = imports, inherits = FALSE), envir = out)
      }
    }
    for (nm in ls(ns, all.names = TRUE)) {
      assign(nm, get(nm, envir = ns, inherits = FALSE), envir = out)
    }
    lockEnvironment(out, bindings = FALSE)
    cache <<- out
    out
  }
})

.rtmb_assigned_vars <- function(expr) {
  vars <- character()
  if (is.call(expr)) {
    if (identical(expr[[1]], as.name("<-")) || identical(expr[[1]], as.name("="))) {
      if (is.name(expr[[2]])) vars <- c(vars, as.character(expr[[2]]))
    } else if (identical(expr[[1]], as.name("{")) && length(expr) > 1L) {
      for (i in seq_along(expr)[-1]) {
        vars <- c(vars, .rtmb_assigned_vars(expr[[i]]))
      }
    } else if (length(expr) > 1L) {
      for (i in seq_along(expr)[-1]) {
        vars <- c(vars, .rtmb_assigned_vars(expr[[i]]))
      }
    }
  }
  unique(vars)
}

.rtmb_setup_env_list <- function(x) {
  if (is.null(x)) return(list())
  if (is.environment(x)) x <- as.list(x, all.names = TRUE)
  if (!is.list(x)) {
    stop("code$setup_env must be a named list or environment.", call. = FALSE)
  }
  if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
    stop("code$setup_env must be a named list.", call. = FALSE)
  }
  x
}

.rtmb_setup_env <- function(env = parent.frame(), setup = NULL, exclude = character()) {
  vars <- if (is.null(setup)) {
    ls(env, all.names = TRUE)
  } else {
    all.vars(setup, functions = FALSE)
  }
  vars <- setdiff(unique(vars), exclude)
  vars <- vars[vapply(vars, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(vars) == 0L) return(list())
  mget(vars, envir = env, inherits = FALSE)
}


#' Create an RTMB_Model Object
#'
#' @description
#' The `rtmb_model` function acts as the core constructor for compiling and combining
#' user-defined data with the model code (defined via `rtmb_code`). It generates an
#' `RTMB_Model` (R6 class) instance, which serves as the foundation for performing
#' Bayesian inference, including Maximum A Posteriori (MAP) estimation, Variational
#' Inference (ADVI), and Markov Chain Monte Carlo (MCMC) sampling.
#'
#' @details
#' \strong{Model Compilation and Pre-checking:}
#' When this function is called, it evaluates the provided data and model blocks.
#' It performs a "sandbox execution" (pre-check) using dummy initial values to dynamically
#' detect common structural errors, such as undefined variables, out-of-bounds indices,
#' or incompatible matrix operations, before proceeding to the computationally expensive
#' Automatic Differentiation (MakeADFun) phase. Cryptic backend errors are caught and
#' translated into user-friendly hints.
#'
#' \strong{Writing AD-Compatible Code (Important):}
#' To ensure the model is differentiable, you must follow specific syntax rules
#' when writing code within \code{rtmb_code}. Avoid discrete branching (\code{if}, \code{ifelse})
#' based on parameters and use numerically stable functions.
#' See \code{\link{rtmb_syntax}} for a detailed guide on Automatic Differentiation requirements.
#'
#' \strong{Initial Values (\code{init}):}
#' Initial parameter values can be specified as a flat numeric vector or a named list.
#' If a partial list is provided, the unspecified parameters are automatically initialized
#' with random values drawn from a uniform distribution on the unconstrained scale.
#' If \code{init = NULL}, all parameters are initialized randomly.
#'
#' \strong{Parameter Labeling (\code{par_names}):}
#' By default, vector or matrix parameters are displayed with numeric indices (e.g., \code{mu[1]}).
#' You can pass a named list of character vectors to \code{par_names} to assign meaningful
#' labels to specific dimensions (e.g., \code{mu[Control]}), vastly improving the readability
#' of summary outputs and trace plots.
#'
#' \strong{Available Inference Methods on the Returned Object:}
#' The returned \code{RTMB_Model} instance provides the following core methods:
#' \itemize{
#'   \item \code{$optimize(...)}: Performs Maximum A Posteriori (MAP) or Maximum Likelihood estimation. Returns a \code{MAP_Fit} object.
#'   \item \code{$sample(...)}: Draws posterior samples using the NUTS (No-U-Turn Sampler) algorithm. Returns an \code{MCMC_Fit} object.
#'   \item \code{$variational(...)}: Performs Mean-field or Full-rank Automatic Differentiation Variational Inference (ADVI). Returns a \code{VB_Fit} object.

#' }
#'
#' @param data A named list containing observation data and constants (e.g., sample size, matrices) used in the model.
#' @param code A model definition block generated by \code{rtmb_code(...)} (must include \code{parameters} and \code{model} blocks).
#' @param par_names A named list of character vectors corresponding to the dimensions of specific parameters (optional).
#' @param init A list or numeric vector of initial values for parameters (optional). If not specified, initialized randomly.
#' @param view Character vector of parameter names to be displayed preferentially at the top when outputting results like \code{summary()} (optional).
#' @param fixed A named list of parameter values to fix (optional). Useful for scoring or plug-in estimation where some parameters (e.g., item parameters) are fixed to known values.
#'
#' @return An \code{RTMB_Model} class instance with a compiled and pre-tested automatic differentiation function.
#'
#' @examples
#' \donttest{
#' # Simulate data for 3 groups
#' set.seed(123)
#' N <- 60
#' group_idx <- sample(1:3, N, replace = TRUE)
#' group_names <- c("Control", "Treatment_A", "Treatment_B")
#'
#' # True group means: Control = 0, Treatment_A = 2, Treatment_B = -1
#' true_means <- c(0, 2, -1)
#' y <- true_means[group_idx] + rnorm(N, mean = 0, sd = 0.5)
#'
#' data_list <- list(N = N, K = 3, group_idx = group_idx, y = y)
#'
#' # Define the model using rtmb_code
#' model_code <- rtmb_code(
#'   parameters = {
#'     mu    = Dim(K)             # Vector of length K (group means)
#'     sigma = Dim(1, lower = 0)  # Scalar (residual standard deviation)
#'   },
#'   model = {
#'     # Priors
#'     for (k in 1:K) mu[k] ~ normal(0, 10)
#'     sigma ~ exponential(1)
#'
#'     # Likelihood
#'     for (i in 1:N) {
#'       y[i] ~ normal(mu[group_idx[i]], sigma)
#'     }
#'   }
#' )
#'
#' # --- 1. Basic Model Creation ---
#' # Create the RTMB_Model object
#' mod_basic <- rtmb_model(
#'   data = data_list,
#'   code = model_code
#' )
#'
#' # Perform Maximum A Posteriori (MAP) estimation
#' map_basic <- mod_basic$optimize()
#' # The summary displays default parameter names: mu[1], mu[2], mu[3]
#' map_basic$summary()
#'
#' # --- 2. Optional: Adding Custom Parameter Names and initial values ---
#' # You can optionally use 'par_names' to assign meaningful labels
#' # to vector or matrix elements for easier interpretation.
#' mod_named <- rtmb_model(
#'   data = data_list,
#'   code = model_code,
#'   init = list(mu = rep(0, 3), sigma = 1),
#'   par_names = list(mu = group_names)
#' )
#'
#' map_named <- mod_named$optimize()
#' # The summary now displays: mu[Control], mu[Treatment_A], mu[Treatment_B]
#' map_named$summary()
#' }
#' @param silent Logical; if TRUE, suppresses diagnostic messages during model creation. Default is FALSE.
#' @param gr_test Logical; if TRUE, evaluate the gradient once after the AD tape is built.
#'   By default, model creation checks that \code{MakeADFun()} can build the AD tape but skips
#'   this extra gradient evaluation for speed.
#'
#' @return An \code{RTMB_Model} class instance with a compiled and pre-tested automatic differentiation function.
#'
#' @export
rtmb_model <- function(data, code, par_names = list(), init = NULL, view = NULL, fixed = NULL, silent = FALSE, gr_test = FALSE) {

  if (is.null(silent)) silent <- FALSE
  if (getOption("BayesRTMB.silent", FALSE)) silent <- TRUE
  if (is.null(gr_test)) gr_test <- FALSE
  if (!is.logical(gr_test) || length(gr_test) != 1L || is.na(gr_test)) {
    stop("gr_test must be TRUE or FALSE.", call. = FALSE)
  }

  old_silent <- options(BayesRTMB.silent = silent)
  on.exit(options(old_silent), add = TRUE)

  if (!"parameters" %in% names(code)) stop("The 'parameters = { ... }' block is required in code.")
  if (!"model" %in% names(code)) stop("The 'model = { ... }' block is required in code.")

  data <- validate_data(data)

  # --- 1. Static validation of Data block (optional) ---
  if ("data" %in% names(code)) {
    data_vars <- all.vars(code$data)
    missing_data <- setdiff(data_vars, names(data))
    if (length(missing_data) > 0) {
      stop(sprintf("[Declaration error in 'data' block]\nThe declared variable '%s' does not exist in the provided data list (%s).",
                   missing_data[1], paste(names(data), collapse = ", ")), call. = FALSE)
    }
  }

  if ("setup" %in% names(code)) {
    # Extract data into a closed environment and evaluate setup code within it.
    # Helper functions used later should be defined here or called with pkg::fun.
    setup_env <- .rtmb_setup_env_list(code$setup_env)
    setup_env_names <- setdiff(names(setup_env), names(data))
    setup_input <- setup_env
    setup_input[names(data)] <- data
    dat_env <- list2env(setup_input, parent = .rtmb_closed_eval_env())
    tryCatch({
      eval(code$setup, envir = dat_env)
    }, error = function(e) {
      msg <- conditionMessage(e)
      missing_fun <- .rtmb_extract_missing_function(msg)
      if (!is.na(missing_fun)) {
        stop(
          sprintf(
            "[Error in 'setup' block] Function '%s' was not found.\n  * Define it inside setup = { ... } or call it as package::%s().",
            missing_fun, missing_fun
          ),
          call. = FALSE
        )
      }
      stop(sprintf("[Error in 'setup' block] An error occurred during data preprocessing:\n  %s", msg), call. = FALSE)
    })
    assigned_setup_vars <- .rtmb_assigned_vars(code$setup)
    drop_setup_env <- setdiff(setup_env_names, assigned_setup_vars)
    if (length(drop_setup_env) > 0L) {
      rm(list = drop_setup_env, envir = dat_env)
    }
    # Return all created/modified variables to the data list after evaluation
    data <- env_to_ordered_list(dat_env, data, code$setup)
  }

  # --- 2. Static validation and evaluation of Parameters block ---
  param_exprs <- as.list(code$parameters)[-1]
  evaluated_par_list <- list()
  eval_env <- list2env(data, parent = parent.frame())

  for (e in param_exprs) {
    if (is.null(e)) next
    if (is.call(e) && (identical(e[[1]], as.name("=")) || identical(e[[1]], as.name("<-")))) {
      v_name <- as.character(e[[2]])
      if (identical(v_name, "lp")) {
        stop(
          "[Error in 'parameters' block]\n",
          "'lp' is a reserved variable used internally for accumulating log probability and cannot be declared as a parameter.",
          call. = FALSE
        )
      }
      d_expr <- e[[3]]

      used_vars <- all.vars(d_expr)
      missing_vars <- setdiff(used_vars, names(data))

      if (length(missing_vars) > 0) {
        stop(sprintf("[Error in 'parameters' block]\nUndefined variable '%s' is used.\n  [Location]: %s\n  * Please check if it is included in the 'data' list.",
                     missing_vars[1], paste(deparse(e), collapse = " ")), call. = FALSE)
      }

      p_obj <- eval(d_expr, envir = eval_env)
      
      # Handle initial values from 'init' argument
      if (!is.null(init)) {
        if (is.list(init) && !is.null(init[[v_name]])) {
          p_obj$init <- init[[v_name]]
        } else if (is.numeric(init)) {
          # If init is a flat vector, it's more complex to handle here without to_constrained logic
          # For now, handle list-based init which is standard for wrappers
        }
      }
      
      evaluated_par_list[[v_name]] <- p_obj
    } else {
      stop("Format error in 'parameters' block: ", paste(deparse(e), collapse = " "))
    }
  }

  for (v_name in names(evaluated_par_list)) {
    p <- evaluated_par_list[[v_name]]
    if (!is.null(p$type) && p$type %in% c("CF_cov", "CF_corr", "cov_matrix", "corr_matrix")) {
      if (length(p$dim) == 3) {
        if (p$dim[2] != p$dim[3]) {
          stop(sprintf("[Parameter definition error] Variable '%s' is specified as type = '%s', but its matrix slices do not have square dimensions (Specified dimension: %s).",
                       v_name, p$type, paste(p$dim, collapse = " x ")), call. = FALSE)
        }
      } else if (length(p$dim) < 2 || p$dim[1] != p$dim[2]) {
        stop(sprintf("[Parameter definition error] Variable '%s' is specified as type = '%s', but it does not have square matrix dimensions (Specified dimension: %s).",
                     v_name, p$type, paste(p$dim, collapse = " x ")), call. = FALSE)
      }
    }
  }

  check_unsupported_funcs <- function(expr) {
    if (is.call(expr)) {
      func_name <- as.character(expr[[1]])
      if (length(func_name) == 1 && func_name %in% c("ifelse", "apply", "sapply", "lapply", "tapply", "vapply")) {
        stop(sprintf("[Syntax error in 'model' block] Function '%s' is used, which is not supported by automatic differentiation (AD).\n  * Workaround for conditionals: Substitute with arithmetic operations like '(condition) * A + (!condition) * B'.\n  * Workaround for loop functions: Rewrite as a standard 'for' loop.", func_name), call. = FALSE)
      }
      lapply(as.list(expr)[-1], check_unsupported_funcs)
    }
  }
  if ("model" %in% names(code)) {
    check_unsupported_funcs(code$model)
  }
  original_code <- code
  code <- inject_transform_adreports(code)

  # --- 3. Dynamic compilation (model, transform, generate) ---
  user_env <- .rtmb_closed_eval_env()

  exec_model_ast <- code$model

  seed_candidates <- names(evaluated_par_list)[vapply(evaluated_par_list, function(p) p$length > 0L, logical(1))]
  ad_seed_name <- if (length(seed_candidates) > 0L) seed_candidates[[1L]] else NULL

  comp_model <- eval(bquote(model_code(.(code$model), env = user_env, ad_seed_name = .(ad_seed_name))))

  comp_transform <- NULL
  if ("transform" %in% names(code)) {
    comp_transform <- eval(bquote(transform_code(.(code$transform), env = user_env, ad_seed_name = .(ad_seed_name))))
  }

  comp_generate <- NULL
  if ("generate" %in% names(code)) {
    comp_generate <- eval(bquote(transform_code(.(code$generate), env = user_env, ad_seed_name = .(ad_seed_name))))
  }

  # --- 4. Pre-check (Sandbox execution) ---
  # --- 4. Pre-check (Sandbox execution) ---
  if (!silent) {
    cat("Pre-checking model code...\n")
  }

  test_unc_list <- lapply(evaluated_par_list, function(p) {
    rnorm(p$unc_length, mean = 0, sd = 0.1)
  })
  test_para <- to_constrained(test_unc_list, evaluated_par_list)

  if (!is.null(comp_transform)) {
    test_tran <- with_rtmb_error_handling({ comp_transform(data, test_para) }, "transform")
    if (is.list(test_tran)) {
      # Check for name collisions between transform outputs and parameter names
      colliding_names <- intersect(names(test_tran), names(evaluated_par_list))
      if (length(colliding_names) > 0) {
        warning(sprintf(
          "[Warning in 'transform' block] The following variable names in 'transform' conflict with parameter names: %s. ",
          paste(colliding_names, collapse = ", ")
        ),
        "This may cause unexpected behavior during model evaluation. ",
        "Consider renaming the transform variables to avoid conflicts.",
        call. = FALSE)
      }
      test_para <- utils::modifyList(test_para, test_tran)
    }
  }

  with_rtmb_error_handling({ 
    comp_model(data, test_para) 
  }, "model")

  if (!is.null(comp_generate)) {
    with_rtmb_error_handling({ comp_generate(data, test_para) }, "generate")
  }

  # --- 5. Create safe wrappers for execution ---
  safe_model <- function(dat, para) { with_rtmb_error_handling({ comp_model(dat, para) }, "model") }
  environment(safe_model) <- list2env(
    list(comp_model = comp_model),
    parent = asNamespace("BayesRTMB")
  )


  safe_transformed <- NULL
  if (!is.null(comp_transform)) {
    safe_transformed <- function(dat, para) {
      res <- with_rtmb_error_handling({ comp_transform(dat, para) }, "transform")
      if (!is.list(res)) stop("The return value must be a 'named list'.")
      res
    }
    environment(safe_transformed) <- list2env(
      list(comp_transform = comp_transform),
      parent = asNamespace("BayesRTMB")
    )
    attr(safe_transformed, "raw_expr") <- code$transform
  }

  safe_generate <- NULL
  if (!is.null(comp_generate)) {
    safe_generate <- function(dat, para) {
      res <- with_rtmb_error_handling({ comp_generate(dat, para) }, "generate")
      if (!is.list(res)) stop("The return value must be a 'named list'.")
      res
    }
    environment(safe_generate) <- list2env(
      list(comp_generate = comp_generate),
      parent = asNamespace("BayesRTMB")
    )
    attr(safe_generate, "raw_expr") <- code$generate
  }

  # --- 6. Pass to R6 class ---
  obj <- RTMB_Model$new(
    data       = data,
    par_list   = evaluated_par_list,
    log_prob   = safe_model,
    transform  = safe_transformed,
    generate   = safe_generate,
    par_names  = par_names,
    init       = init,
    view       = view,
    code       = original_code,
    gr_test    = gr_test
  )


  if (!is.null(fixed)) {
    obj <- obj$fixed_model(fixed)
  }



  return(obj)
}

report <- function(...) invisible(NULL)

# --- Helper functions for AST exploration ---

.rtmb_call_name <- function(expr) {
  if (!is.call(expr)) return(NA_character_)
  head <- expr[[1]]
  if (is.name(head)) return(as.character(head))
  if (is.call(head) &&
      (identical(head[[1]], as.name("::")) || identical(head[[1]], as.name(":::"))) &&
      length(head) >= 3L &&
      is.name(head[[3]])) {
    return(as.character(head[[3]]))
  }
  NA_character_
}

.rtmb_is_call_to <- function(expr, func_name) {
  identical(.rtmb_call_name(expr), func_name)
}

# Check if a specific function call (e.g., "report") is included in the expression
has_function_call <- function(expr, func_name) {
  if (is.atomic(expr) || is.name(expr)) return(FALSE)
  if (is.call(expr)) {
    if (.rtmb_is_call_to(expr, func_name)) return(TRUE)
    return(any(sapply(as.list(expr)[-1], has_function_call, func_name = func_name)))
  }
  return(FALSE)
}

extract_report_vars <- function(expr) {
  vars <- character()
  if (is.atomic(expr) || is.name(expr)) return(vars)
  if (is.call(expr)) {
    if (.rtmb_is_call_to(expr, "report")) {
      for (i in seq_along(expr)[-1]) {
        if (is.name(expr[[i]])) vars <- c(vars, as.character(expr[[i]]))
      }
    } else {
      for (i in seq_along(expr)[-1]) {
        vars <- c(vars, extract_report_vars(expr[[i]]))
      }
    }
  }
  return(unique(vars))
}

# Remove report() calls from AST to avoid runtime errors
remove_report_calls <- function(expr) {
  if (is.atomic(expr) || is.name(expr)) return(expr)
  if (is.call(expr)) {
    if (.rtmb_is_call_to(expr, "report")) return(NULL)
    new_args <- lapply(as.list(expr), remove_report_calls)
    new_args <- Filter(Negate(is.null), new_args)
    if (length(new_args) > 0) return(as.call(new_args)) else return(NULL)
  }
  return(expr)
}

# Extract all assigned variable names
get_assigned_vars <- function(expr) {
  vars <- character()
  if (is.atomic(expr) || is.name(expr)) return(vars)
  if (is.call(expr)) {
    if (identical(expr[[1]], as.name("<-")) || identical(expr[[1]], as.name("="))) {
      if (is.name(expr[[2]])) vars <- c(vars, as.character(expr[[2]]))
    }
    for (i in seq_along(expr)[-1]) {
      vars <- c(vars, get_assigned_vars(expr[[i]]))
    }
  }
  return(unique(vars))
}

# Inject ADREPORT into the model block based on the transform content
inject_transform_adreports <- function(code) {
  if (is.null(code$transform)) return(code)

  # 1. Check if report() is explicitly specified
  explicit_vars <- extract_report_vars(code$transform)

  if (length(explicit_vars) > 0) {
    # If explicit, target only those variables. Remove the report() call from transform
    target_vars <- explicit_vars
    code$transform <- remove_report_calls(code$transform)
  } else {
    # If not explicitly specified, target all variables assigned in transform
    target_vars <- get_assigned_vars(code$transform)
  }

  # 2. Add ADREPORT(var) to the end of the model block
  if (length(target_vars) > 0 && !is.null(code$model)) {
    ad_exprs <- lapply(target_vars, function(v) call("ADREPORT", as.name(v)))
    if (is.call(code$model) && identical(code$model[[1]], as.name("{"))) {
      # Inject into { ... }
      code$model <- as.call(c(as.list(code$model), ad_exprs))
    }
  }
  return(code)
}

# Extract the function name from common "function not found" error messages.
.rtmb_extract_missing_function <- function(msg) {
  patterns <- c(
    "could not find function [\"']([^\"']+)[\"']",
    "[\"']([^\"']+)[\"'] .*could not find function",
    "\u95a2\u6570 [\"']([^\"']+)[\"'] .*\u898b\u3064"
  )
  for (pat in patterns) {
    hit <- regexec(pat, msg)
    val <- regmatches(msg, hit)[[1]]
    if (length(val) >= 2L) return(val[2])
  }
  NA_character_
}

with_rtmb_error_handling <- function(expr, block_name) {
  tryCatch({
    expr
  }, error = function(e) {
    msg <- conditionMessage(e)
    call_expr <- conditionCall(e)

    call_str <- if (!is.null(call_expr)) {
      code_snippet <- paste(deparse(call_expr), collapse = " ")
      if (nchar(code_snippet) > 100) paste0(substr(code_snippet, 1, 97), "...") else code_snippet
    } else {
      "Unknown location"
    }

    # 1. Use of non-existent variables
    if (grepl("object '.*' not found", msg)) {
      var_name <- sub(".*object '(.*)'.*", "\\1", msg)
      stop(sprintf("[Error in '%s' block] Undefined variable '%s' is used.\n  [Location]: %s\n  * Please check if it is included in 'data' or 'parameters'.",
                   block_name, var_name, call_str), call. = FALSE)
    }

    # 2. Use of non-existent functions
    func_name <- .rtmb_extract_missing_function(msg)
    if (!is.na(func_name)) {
      stop(sprintf("[Error in '%s' block] Function '%s' was not found.\n  [Location]: %s\n  * Define it inside setup = { ... } or call it as package::%s().",
                   block_name, func_name, call_str, func_name), call. = FALSE)
    }

    # 3. Index error
    if (grepl("subscript out of bounds", msg)) {
      stop(sprintf("[Error in '%s' block] Array or matrix subscript out of bounds.\n  [Location]: %s",
                   block_name, call_str), call. = FALSE)
    }

    if (grepl("categorical_logit\\(\\) category index is outside", msg)) {
      stop(sprintf("[Error in '%s' block] %s\n  [Location]: %s",
                   block_name, msg, call_str), call. = FALSE)
    }

    # 4. Destruction of advector / Failure of type conversion
    if (grepl("lost class attribute", msg) || grepl("Invalid argument to 'advector'", msg) ||
        grepl("non-numeric argument to mathematical function", msg)) {
      stop(sprintf("[Error in '%s' block] Automatic differentiation (AD) type was lost or an invalid operation was performed.\n  [Location]: %s\n  * Possible cause: Assigning values to an empty vector initialized with 'numeric(N)' (turns into a list type).\n  * Solutions:\n    1. Vectorize the calculation to avoid loops.\n    2. Initialize while preserving AD type, e.g., 'vec <- rep(theta[1]*0, N)'.\n    3. Use 'advector(numeric(N))'.",
                   block_name, call_str), call. = FALSE)
    }

    # 5. Matrix dimension mismatch
    if (grepl("non-conformable arguments", msg)) {
      stop(sprintf("[Error in '%s' block] Non-conformable arguments in matrix multiplication (%%*%%) or operations.\n  [Location]: %s\n  * Note: R's standard 'automatic recycling rule' for vectors of different lengths does not apply to RTMB matrix operations. Ensure dimensions match exactly or use rep() to align sizes explicitly.",
                   block_name, call_str), call. = FALSE)
    }

    # Other errors
    stop(sprintf("[Error in '%s' block]\n  [Message]: %s\n  [Location]: %s",
                 block_name, msg, call_str), call. = FALSE)
  })
}

# Format MakeADFun backend errors with BayesRTMB hints.
.rtmb_format_makeadfun_error <- function(msg, context = "MakeADFun") {
  hints <- character(0)

  missing_fun <- .rtmb_extract_missing_function(msg)
  if (!is.na(missing_fun)) {
    hints <- c(
      hints,
      sprintf("Function '%s' was not found.", missing_fun),
      sprintf("Define '%s' inside setup = { ... } or call it as package::%s().", missing_fun, missing_fun)
    )
  }

  if (grepl("not a valid 'advector'|illegal operation|lost class attribute|Invalid argument to 'advector'", msg)) {
    hints <- c(
      hints,
      "An automatic-differentiation (AD) object was used after an invalid operation.",
      "Common causes in rtmb_code():",
      "  1. Indexing a matrix as x[t] when you intended a row, x[t, ].",
      "  2. Passing a logit/probability vector with the wrong length to categorical_logit().",
      "  3. Assigning past the end of an AD vector, such as vec[t + 1] when vec has length Trial_t.",
      "  4. Growing or rebuilding an AD vector inside a loop; initialize with rep(theta[1] * 0, N) or a fixed-size matrix first."
    )
  }

  if (grepl("incorrect number of dimensions", msg)) {
    hints <- c(
      hints,
      "An object was indexed with the wrong number of dimensions.",
      "For example, use x[t, ] for a matrix row and x[t, c] for one matrix cell."
    )
  }

  if (grepl("subscript out of bounds|attempt to select", msg)) {
    hints <- c(
      hints,
      "A vector, matrix, or array index is outside its valid range.",
      "Check that category labels match the length of the probability/logit vector."
    )
  }

  if (grepl("Comparison is generally unsafe for AD types", msg)) {
    hints <- c(
      hints,
      "A comparison was applied to an AD value.",
      "Avoid pmin(), pmax(), ifelse(), or if-statements on parameters/transformed parameters; clip data before rtmb_model() when possible."
    )
  }

  hint_text <- if (length(hints) > 0L) {
    paste0("\n[Hint]:\n  ", paste(hints, collapse = "\n  "))
  } else {
    ""
  }

  paste0("Failed to setup ", context, ".\n[Error]: ", msg, hint_text)
}

#' Internal function to convert an environment to a list while maintaining original order
#' @keywords internal
env_to_ordered_list <- function(env, orig_list, code_ast = NULL) {
  env_list <- as.list(env)
  orig_names <- names(orig_list)

  if (!is.null(code_ast)) {
    extract_assigned_vars <- function(expr) {
      vars <- character()
      if (is.call(expr)) {
        if (identical(expr[[1]], as.name("<-")) || identical(expr[[1]], as.name("="))) {
          if (is.name(expr[[2]])) vars <- c(vars, as.character(expr[[2]]))
        } else if (identical(expr[[1]], as.name("{"))) {
          if (length(expr) > 1) {
            for (i in seq_along(expr)[-1]) {
              vars <- c(vars, extract_assigned_vars(expr[[i]]))
            }
          }
        }
      }
      return(unique(vars))
    }

    assigned_vars <- extract_assigned_vars(code_ast)
    assigned_vars <- intersect(assigned_vars, names(env_list))
    new_names <- intersect(assigned_vars, setdiff(names(env_list), orig_names))
    other_new_names <- setdiff(names(env_list), c(orig_names, new_names))
    new_names <- c(new_names, other_new_names)
  } else {
    new_names <- setdiff(names(env_list), orig_names)
  }

  return(env_list[c(intersect(orig_names, names(env_list)), new_names)])
}

#' Internal function to search AST and inject namespace to package functions
#' @param expr Unevaluated code (AST)
#' @param pkg Target package name
#' @keywords internal
inject_namespace <- function(expr, pkg = "BayesRTMB") {
  if (is.atomic(expr) || is.name(expr)) {
    return(expr)
  }

  if (is.call(expr)) {
    if (identical(expr[[1]], quote(`::`)) || identical(expr[[1]], quote(`:::`))) {
      return(expr)
    }

    func_name <- as.character(expr[[1]])

    if (length(func_name) == 1) {
      if (identical(func_name, "c")) {
        expr[[1]] <- call(":::", as.name(pkg), as.name(".rtmb_c"))
      } else if (exists(func_name, envir = asNamespace(pkg), inherits = FALSE)) {
        expr[[1]] <- call(":::", as.name(pkg), as.name(func_name))
      }
    }

    if (length(expr) > 1) {
      for (i in seq_along(expr)[-1]) {
        if (!identical(expr[[i]], quote(expr=))) {
          expr[[i]] <- inject_namespace(expr[[i]], pkg)
        }
      }
    }
  }

  return(expr)
}

#' Define an RTMB Model with Stan-like Syntax
#'
#' @description
#' \code{rtmb_code} is the primary interface for defining Bayesian or Frequentist models
#' within the RTMB framework. It allows you to organize model logic into functional
#' blocks while utilizing a Stan-like \code{~} operator for specifying priors and likelihoods.
#'
#' @details
#' \strong{Key Blocks and Logic:}
#' \itemize{
#'   \item \code{setup}: Pre-compilation data processing.
#'   \item \code{parameters}: Declaration of estimated parameters. See \code{\link{parameter_types}} for available constraints.
#'   \item \code{transform}: Definition of intermediate variables. Use \code{\link{math_functions}} for stable AD calculations.
#'   \item \code{model}: Likelihood and Priors. Refer to \code{\link{distributions}} for available sampling statements.
#'   \item \code{generate}: Post-hoc calculation of predictions or diagnostics.
#' }
#'
#' \strong{Automatic Differentiation (AD) Note:}
#' All code within \code{parameters}, \code{transform}, and \code{model} is
#' automatically differentiated by RTMB. To ensure numerical stability, use
#' provided utility functions like \code{\link{log_sum_exp}} or \code{\link{inv_logit}}
#' instead of raw algebraic implementations.
#'
#' @param ... Model definition blocks.
#' @return A list of unevaluated code blocks.
#' @seealso
#' \code{\link{rtmb_model}} for building the model object.
#' \code{\link{parameter_types}} for constraining parameters.
#' \code{\link{distributions}} for likelihood functions.
#' \code{\link{math_functions}} for numerical utilities.
#'
#' @examples
#' \donttest{
#' # Simulate data for a linear regression
#' set.seed(123)
#' N <- 100
#' x <- rnorm(N)
#' y <- 2.0 + 1.5 * x + rnorm(N, mean = 0, sd = 1)
#' data_list <- list(N = N, x = x, y = y)
#'
#' # Define the model using rtmb_code
#' code <- rtmb_code(
#'   setup = {
#'     # Center the predictor variable (executed once)
#'     x_centered <- x - mean(x)
#'   },
#'   parameters = {
#'     # Define parameters and their constraints
#'     alpha = Dim(1)
#'     beta  = Dim(1)
#'     sigma = Dim(1, lower = 0)
#'   },
#'   transform = {
#'     # Calculate the linear predictor
#'     mu <- alpha + beta * x_centered
#'   },
#'   model = {
#'     # Priors
#'     alpha ~ normal(0, 10)
#'     beta  ~ normal(0, 10)
#'     sigma ~ exponential(1)
#'
#'     # Likelihood (Vectorized)
#'     y ~ normal(mu, sigma)
#'   },
#'   generate = {
#'     # Calculate generated quantities
#'     y_pred <- mu
#'
#'     # Must return a named list
#'     list(y_pred = y_pred)
#'   }
#' )
#'
#' # Create the model object
#' mod <- rtmb_model(data = data_list, code = code)
#'
#' # Fit the model using MAP estimation
#' map_res <- mod$optimize()
#' map_res$summary(pars = c("alpha", "beta", "sigma"))
#'
#' # The generated quantity 'y_pred' can also be summarized
#' map_res$summary("y_pred", max_rows = 5)
#' }
#' @export
rtmb_code <- function(...) {
  args <- as.list(match.call())[-1]
  code_list <- list()
  valid_blocks <- c("setup", "data", "parameters", "transform", "model", "generate")

  if (!is.null(names(args))) {
    for (name in names(args)) {
      if (name == "") next
      if (name %in% valid_blocks) {
        code_list[[name]] <- args[[name]]
      } else {
        stop(sprintf("Unknown block '%s' is specified. Valid blocks are %s.",
                     name, paste(valid_blocks, collapse = ", ")), call. = FALSE)
      }
    }

  } else if (length(args) == 1 && is.null(names(args))) {
    ast <- args[[1]]
    if (is.call(ast) && identical(ast[[1]], as.name("{"))) {
      exprs <- as.list(ast)[-1]

      for (e in exprs) {
        if (is.null(e)) next

        if (is.call(e)) {
          block_name <- as.character(e[[1]])
          if (block_name %in% valid_blocks) {
            code_list[[block_name]] <- if (length(e) >= 2) e[[2]] else quote({})
          } else {
            stop(sprintf("Unknown block '%s' is specified. Valid blocks are %s.",
                         block_name, paste(valid_blocks, collapse = ", ")), call. = FALSE)
          }
        } else {
          stop(sprintf("Invalid syntax in rtmb_code: '%s'",
                       paste(deparse(e), collapse = " ")), call. = FALSE)
        }
      }
    } else {
      stop("Invalid syntax in rtmb_code(). Separate blocks with commas, or enclose the entire definition in {}.", call. = FALSE)
    }
  } else {
    stop("Invalid syntax in rtmb_code(). Separate blocks with commas, or enclose the entire definition in {}.", call. = FALSE)
  }

  code_list$env <- parent.frame()
  return(code_list)
}


# Build the AD seed expression injected into generated model-code functions.
.rtmb_ad_seed_expr <- function(ad_seed_name = NULL) {
  if (is.null(ad_seed_name) || !nzchar(ad_seed_name)) {
    return(quote(.rtmb_ad_seed <- NULL))
  }
  substitute(.rtmb_ad_seed <- SEED, list(SEED = as.name(ad_seed_name)))
}

model_code <- function(expr, env = parent.frame(), ad_seed_name = NULL) {
  raw_expr <- substitute(expr)

  if (is.name(raw_expr)) {
    evaluated <- eval(raw_expr, envir = env)
    if (is.language(evaluated)) {
      raw_expr <- evaluated
    }
  }

  # --- Inject namespace into standard math functions etc. ---
  raw_expr <- inject_namespace(raw_expr, pkg = "BayesRTMB")

  rewrite_formula <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("~"))) {
        target <- x[[2]]
        dist_call <- x[[3]]

        if (!is.call(dist_call)) {
          stop(sprintf("The right side of '~' must be a distribution function call. (e.g., %s( ... ))",
                       as.character(dist_call)), call. = FALSE)
        }

        dist_name_raw <- dist_call[[1]]
        if (is.call(dist_name_raw) && identical(dist_name_raw[[1]], as.name(":::"))) {
          dist_name <- as.character(dist_name_raw[[3]])
        } else {
          dist_name <- as.character(dist_name_raw)
        }
        dist_args <- as.list(dist_call[-1])

        name_lpdf <- paste0(dist_name, "_lpdf")
        name_lpmf <- paste0(dist_name, "_lpmf")

        exists_in_pkg <- function(fname) {
          if (length(fname) > 1) return(FALSE)
          requireNamespace("BayesRTMB", quietly = TRUE) &&
            exists(fname, mode = "function", envir = asNamespace("BayesRTMB"), inherits = FALSE)
        }

        use_pkg_namespace <- FALSE

        if (exists_in_pkg(name_lpdf)) {
          actual_name <- name_lpdf
          use_pkg_namespace <- TRUE
        } else if (exists_in_pkg(name_lpmf)) {
          actual_name <- name_lpmf
          use_pkg_namespace <- TRUE
        } else if (exists(name_lpdf, mode = "function", envir = env)) {
          actual_name <- name_lpdf
        } else if (exists(name_lpmf, mode = "function", envir = env)) {
          actual_name <- name_lpmf
        } else {
          actual_name <- name_lpdf # Default
        }

        if (use_pkg_namespace) {
          func_call <- call(":::", as.name("BayesRTMB"), as.name(actual_name))
        } else {
          func_call <- as.name(actual_name)
        }

        new_call <- as.call(c(
          as.name("<-"),
          as.name("lp"),
          as.call(c(
            as.name("+"),
            as.name("lp"),
            as.call(c(func_call, target, dist_args))
          ))
        ))
        return(new_call)
      }
      x[] <- lapply(x, rewrite_formula)
    }
    return(x)
  }

  processed_expr <- rewrite_formula(raw_expr)

  if (is.call(processed_expr) && identical(processed_expr[[1]], as.name("{"))) {
    expr_elements <- as.list(processed_expr)[-1]
  } else {
    expr_elements <- list(processed_expr)
  }

  body_list <- c(
    list(as.name("{")),
    quote(RTMB::getAll(dat, par)),
    list(.rtmb_ad_seed_expr(ad_seed_name)),
    quote(lp <- 0),
    unname(expr_elements),
    quote(return(lp))
  )
  body_expr <- as.call(body_list)

  args <- as.pairlist(alist(dat = , par = ))
  fn_expr <- call("function", args, body_expr)

  log_prob_fn <- eval(fn_expr, envir = env)
  attr(log_prob_fn, "raw_expr") <- raw_expr

  return(log_prob_fn)
}


#' Transformed Code Wrapper for RTMB
#'
#' @param expr A block of code containing calculations for transformed parameters.
#' @param env Environment to assign to the generated function.
#' @param ad_seed_name Optional parameter name used internally to seed RTMB's AD type.
#' @return A function taking (dat, par) that returns a named list.
transform_code <- function(expr, env = parent.frame(), ad_seed_name = NULL) {
  raw_expr <- substitute(expr)

  raw_expr <- inject_namespace(raw_expr, pkg = "BayesRTMB")

  if (is.call(raw_expr) && identical(raw_expr[[1]], as.name("{"))) {
    expr_elements <- as.list(raw_expr)[-1]
  } else {
    expr_elements <- list(raw_expr)
  }

  if (length(expr_elements) == 0) {
    body_list <- c(
      list(as.name("{")),
      quote(RTMB::getAll(dat, par)),
      list(.rtmb_ad_seed_expr(ad_seed_name)),
      list(quote(return(list())))
    )
    body_expr <- as.call(body_list)
    args <- as.pairlist(alist(dat = , par = ))
    fn_expr <- call("function", args, body_expr)
    fn <- eval(fn_expr, envir = env)
    attr(fn, "raw_expr") <- raw_expr
    return(fn)
  }

  last_elem <- expr_elements[[length(expr_elements)]]

  is_explicit_return <- FALSE
  if (is.call(last_elem)) {
    if (identical(last_elem[[1]], as.name("list"))) {
      is_explicit_return <- TRUE
      ret_call <- last_elem
    } else if (identical(last_elem[[1]], as.name("return"))) {
      is_explicit_return <- TRUE
      ret_call <- last_elem[[2]]
    }
  } else if (is.name(last_elem)) {
    is_explicit_return <- TRUE
    ret_call <- last_elem
  }

  # --- report() AST
  reported_vars <- character()

  strip_report <- function(e) {
    if (is.atomic(e) || is.name(e)) return(e)
    if (is.call(e)) {
      if (identical(e[[1]], as.name("{"))) {
        new_args <- list()
        for (i in seq_along(e)[-1]) {
          sub_e <- e[[i]]
          if (is.call(sub_e) && .rtmb_is_call_to(sub_e, "report")) {
            if (length(sub_e) > 1) {
              for (j in seq_along(sub_e)[-1]) {
                if (is.name(sub_e[[j]])) {
                  reported_vars <<- c(reported_vars, as.character(sub_e[[j]]))
                }
              }
            }
            next
          }
          new_args <- c(new_args, list(strip_report(sub_e)))
        }
        return(as.call(c(list(as.name("{")), new_args)))
      }
      for (i in seq_along(e)[-1]) {
        e[[i]] <- strip_report(e[[i]])
      }
    }
    return(e)
  }

  raw_expr_stripped <- strip_report(raw_expr)
  reported_vars <- unique(reported_vars)

  if (is.call(raw_expr_stripped) && identical(raw_expr_stripped[[1]], as.name("{"))) {
    expr_elements <- as.list(raw_expr_stripped)[-1]
  } else {
    expr_elements <- if (is.null(raw_expr_stripped)) list() else list(raw_expr_stripped)
  }

  if (is_explicit_return) {
    expr_elements <- expr_elements[-length(expr_elements)]
  } else {
    # report()
    if (length(reported_vars) > 0) {
      defined_vars <- reported_vars
    } else {
      find_assignments <- function(x) {
        if (is.call(x)) {
          if (identical(x[[1]], as.name("<-")) || identical(x[[1]], as.name("="))) {
            if (is.name(x[[2]])) return(as.character(x[[2]]))
          }
          return(unique(unlist(lapply(as.list(x), find_assignments))))
        }
        return(NULL)
      }
      defined_vars <- find_assignments(raw_expr)
    }

    ret_list_args <- lapply(defined_vars, as.name)
    names(ret_list_args) <- defined_vars
    ret_call <- as.call(c(list(as.name("list")), ret_list_args))
  }

  body_list <- c(
    list(as.name("{")),
    quote(RTMB::getAll(dat, par)),
    list(.rtmb_ad_seed_expr(ad_seed_name)),
    unname(expr_elements),
    list(as.call(list(as.name("return"), ret_call)))
  )
  body_expr <- as.call(body_list)

  args <- as.pairlist(alist(dat = , par = ))
  fn_expr <- call("function", args, body_expr)
  fn <- eval(fn_expr, envir = env)
  attr(fn, "raw_expr") <- raw_expr

  return(fn)
}

#' Code block for parameter definitions
#'
#' @description
#' Defines the list of parameters to be passed to `rtmb_model`.
#' By writing in the format `variable_name = Dim(...)` within a block enclosed in `{}`,
#' evaluation is delayed, enabling strict error checking during model construction.
#'
#' @param expr A parameter definition expression enclosed in `{}`.
#' @return A lazily evaluated list of parameters.
parameters_code <- function(expr) {
  ast <- substitute(expr)

  if (!is.call(ast) || !identical(ast[[1]], as.name("{"))) {
    stop("parameters_code must receive a block enclosed in {}.\nExample: parameters_code({ beta = Dim(K) })", call. = FALSE)
  }

  exprs <- as.list(ast)[-1]
  param_list <- list()

  for (e in exprs) {
    if (is.null(e)) next

    if (is.call(e) && (identical(e[[1]], as.name("=")) || identical(e[[1]], as.name("<-")))) {
      var_name <- as.character(e[[2]])
      param_list[[var_name]] <- e[[3]]
    } else {
      stop(sprintf("Invalid syntax in parameters_code: '%s'\nPlease describe in the format 'variable_name = Dim(...)'.",
                   paste(deparse(e), collapse = " ")), call. = FALSE)
    }
  }

  return(param_list)
}


#' Pre-validation of data and parameters
#'
#' @description
#' Before passing data to `rtmb_model`, it checks for the presence of R-specific data types (such as data.frame)
#' and missing values, and outputs errors or warnings accordingly.
#'
#' @param dat_list A list of data to be passed to the model (usually containing matrices or numeric vectors).
#'
#' @return Returns invisible `NULL` on success. Interrupts execution with `stop()` or issues `warning()` if issues are found.
#'
#' @details
#' RTMB's automatic differentiation engine requires pure `matrix` or `numeric` types.
#' If a user mistakenly passes a `data.frame`, incomprehensible errors occur during the construction of the computation graph.
#' This function catches such issues in advance.
validate_data <- function(dat_list) {
  for (name in names(dat_list)) {
    if (inherits(dat_list[[name]], "rtmb_setup_df")) next
    if (is.data.frame(dat_list[[name]])) {
      dat_list[[name]] <- as.matrix(dat_list[[name]])
    }
  }
  for (name in names(dat_list)) {
    x <- dat_list[[name]]
    if (inherits(x, "formula") || is.list(x) || is.function(x)) next
  }
  return(dat_list)
}


#' Safe RTMB model construction (with error message translation)
#'
#' @description
#' Wraps `rtmb_model` to translate cryptic error messages generated during `MakeADFun`
#' execution (originating from C++/RTMB) into user-friendly hints.
#'
#' @param data A list of data that has passed `validate_data`.
#' @param code Code blocks for likelihood and priors defined with `model_code({})`.
#' @param par_names A list of variable names corresponding to the dimensions of each parameter (optional).
#' @param init A list or numeric vector of initial values (optional).
#' @param view Character vector of parameter names to be displayed preferentially in summary outputs (optional).

#'
#' @return Returns an `rtmb_model` object (R6 class instance) upon successful compilation.
#'
#' @export
safe_rtmb_model <- function(data, code, par_names = list(), init = NULL, view = NULL) {
  data <- validate_data(data)

  tryCatch({
    obj <- rtmb_model(data = data, code = code, par_names = par_names,
                      init = init, view = view)
    return(obj)

  }, error = function(e) {
    err_msg <- conditionMessage(e)

    if (grepl("lost class attribute", err_msg) || grepl("Invalid argument to 'advector'", err_msg)) {
      stop("[Model construction error] Automatic differentiation (AD) type was lost during execution.\n",
           "  * Cause: You may be performing sequential assignments within a 'for' loop to an empty 'array()' or using 'rep()' inside 'model_code'.\n",
           "  * Solution: Use vectorized matrix calculations, or use 'numeric(M)' etc. for initialization.\n",
           "[Original error]: ", err_msg, call. = FALSE)
    }

    if (grepl("illegal operation", err_msg) || grepl("not a valid 'advector'", err_msg)) {
      stop("[Model construction error] An invalid operation or type mismatch occurred.\n",
           "  * Cause: Dimensions may not match around matrix multiplication (%*%), or type specification in Dim() (e.g., CF_cov, lower_tri) may contradict the actual matrix size.\n",
           "[Original error]: ", err_msg, call. = FALSE)
    }

    if (grepl("requires numeric/complex matrix/vector", err_msg)) {
      stop("[Data type error] A different type was passed where a numeric vector or matrix is required.\n",
           "  * Cause: Are observation data or indices passed as data.frame or list?\n",
           "[Original error]: ", err_msg, call. = FALSE)
    }

    stop("[Unexpected error]\n", err_msg, call. = FALSE)
  })
}
