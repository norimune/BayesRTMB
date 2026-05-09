#' Common Features and Arguments of RTMB Wrapper Functions
#'
#' @description
#' The RTMB wrapper functions (`rtmb_lm`, `rtmb_glm`, `rtmb_glmer`, `rtmb_fa`, etc.)
#' share a unified interface designed to make Bayesian and Frequentist inference
#' accessible through familiar R formulas and standard model specifications.
#'
#' @details
#' All wrapper functions in this package are built upon the same core engine.
#' This ensures that regardless of the model type, the workflow for estimation,
#' summary, and expansion remains consistent.
#'
#' \strong{1. Unified Inference Methods:}
#' Every model object returned by a wrapper function provides the following methods:
#' \itemize{
#'   \item \code{$optimize()}: Performs MAP estimation (comparable to MLE).
#'   \item \code{$sample()}: Performs MCMC sampling (NUTS) for full Bayesian inference.
#'   \item \code{$variational()}: Performs Variational Inference (ADVI) for fast posterior approximation.
#' }
#'
#' \strong{2. Regularization and Variable Selection:}
#' For rtmb_lm, rtmb_glm, and rtmb_glmer, you can handle high-dimensional predictors (where the number of variables is large relative to the sample size) using the penalty argument:
#' \itemize{
#'   \item \code{"none"}: Standard flat or weakly informative priors.
#'   \item \code{"rhs"}: Regularized Horseshoe prior for continuous shrinkage.
#'   \item \code{"ssp"}: Spike-and-Slab prior for sparse variable selection.
#' }
#' \emph{Note: When using regularization, you must specify \code{y_range = c(min, max)}
#' to let the model set appropriate global scales for the priors.}
#'
#' \strong{3. Weakly Informative Priors (\code{y_range}):}
#' By providing the theoretical range of your response variable via \code{y_range},
#' the wrappers automatically construct "Weakly Informative Priors". These priors
#' are designed to be broad enough to cover any reasonable value but narrow enough
#' to stabilize the estimation and prevent the sampler from wandering into
#' non-sensical parameter space.
#'
#' \strong{4. Fixed vs. Random Effects:}
#' For mixed-effect models (e.g., \code{rtmb_glmer}), random effects are marginalized
#' using the Laplace approximation during \code{$optimize(laplace = TRUE)}.
#' When using \code{$sample()}, random effects are treated as unknown parameters
#' and sampled alongside fixed effects.
#'
#' \strong{5. Null Model Creation (\code{null}):}
#' You can specify a \code{null} argument (e.g., \code{null = "x1 ~ normal(0, 0.1)"})
#' in the wrappers to simultaneously create a restricted version of your model.
#' This is particularly useful for computing Bayes Factors or performing
#' model comparisons.
#'
#' @name rtmb_wrappers
#' @family wrappers
NULL

# Internal helper to merge priors without removing NULL values
.merge_prior <- function(default, user) {
  if (is.null(user)) return(default)
  res <- user
  for (name in names(default)) {
    if (!(name %in% names(user))) {
      res[[name]] <- default[[name]]
    }
  }
  return(res)
}


#' Specify a uniform or manual prior
#'
#' @param Intercept_sd Standard deviation for the intercept prior (Normal). Default is NULL (flat).
#' @param b_sd Standard deviation for the coefficients prior (Normal). Default is NULL (flat).
#' @param sigma_rate Rate for the residual standard deviation prior (Exponential). Default is NULL (flat).
#' @param tau_rate Rate for the random effects standard deviation prior (Exponential). Default is NULL (flat).
#' @param ... Optional hyperparameters
#' @return A list with class "rtmb_prior"
#' @export
prior_uniform <- function(Intercept_sd = NULL, b_sd = NULL, mu_sd = NULL, sigma_rate = NULL, tau_rate = NULL, ...) {
  res <- list(type = "uniform", Intercept_sd = Intercept_sd, b_sd = b_sd, mu_sd = mu_sd, sigma_rate = sigma_rate, tau_rate = tau_rate, ...)
  class(res) <- "rtmb_prior"
  return(res)
}

#' Specify a weakly informative prior
#'
#' @param sd_ratio Ratio of the prior standard deviation to the half-range of the response variable. Default is 0.5.
#' @param max_beta Maximum expected effect size. Default is 1.0.
#' @param ... Optional hyperparameters for other distributions.
#' @return A list with class "rtmb_prior"
#' @export
prior_weak <- function(sd_ratio = 0.5, max_beta = 1.0, ...) {
  res <- list(type = "weak", sd_ratio = sd_ratio, max_beta = max_beta, ...)
  class(res) <- "rtmb_prior"
  return(res)
}

#' Specify a JZS (Jeffreys-Zellner-Siow) prior
#'
#' @param r Scale factor for the Cauchy prior on effect sizes. Default is 0.707 (sqrt(2)/2).
#' @param ... Optional hyperparameters
#' @return A list with class "rtmb_prior"
#' @export
prior_jzs <- function(r = 0.707, ...) {
  res <- list(type = "jzs", r = r, ...)
  class(res) <- "rtmb_prior"
  return(res)
}

#' Specify a Spike-and-Slab prior for variable selection
#'
#' @param ssp_ratio Prior probability of inclusion for each variable. Default is 0.25.
#' @param max_beta Maximum expected effect size. Default is 1.0.
#' @param ... Optional hyperparameters for other distributions.
#' @return A list with class "rtmb_prior"
#' @export
prior_ssp <- function(ssp_ratio = 0.25, max_beta = 1.0, ...) {
  res <- list(type = "ssp", ssp_ratio = ssp_ratio, max_beta = max_beta, ...)
  class(res) <- "rtmb_prior"
  return(res)
}

#' Specify a Regularized Horseshoe prior for continuous shrinkage
#'
#' @param expected_vars Expected number of non-zero variables. Default is 3.
#' @param slab_scale Scale parameter for the slab distribution. Default is 2.0.
#' @param slab_df Degrees of freedom for the slab distribution. Default is 4.0.
#' @param ... Optional hyperparameters for other distributions.
#' @return A list with class "rtmb_prior"
#' @export
prior_rhs <- function(expected_vars = 3, slab_scale = 2.0, slab_df = 4.0, ...) {
  res <- list(type = "rhs", expected_vars = expected_vars, slab_scale = slab_scale, slab_df = slab_df, ...)
  class(res) <- "rtmb_prior"
  return(res)
}

#' Specify a JZS (Jeffrey-Zellner-Siow) prior for t-tests
#'
#' @param r Numeric; scale parameter for the Cauchy prior on effect size delta. Default is 0.707.
#' @return A list with class "rtmb_prior"
#' @export
prior_jzs <- function(r = 0.707) {
  res <- list(type = "jzs", r = r)
  class(res) <- "rtmb_prior"
  return(res)
}

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
rtmb_table <- function(x, y = NULL, data = NULL, correct = TRUE, prior = prior_uniform(), fixed = NULL, ...) {

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
      report(mu)
      report(chisq_val)
    },
    model = {
      # Likelihood
      Y ~ multinomial(N, p)

      # Prior (Flat Dirichlet by default)
      p ~ dirichlet(rep(1, R * C))
    },
    generate = {
      # We only generate what is not already in parameters/report
      mu <- p * N
      chisq_val <- sum((Y - mu)^2 / mu)
      list(mu = mu, chisq_val = chisq_val)
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

#' RTMB-based Log-linear Analysis (Contingency Table)
#'
#' @description
#' `rtmb_loglinear` performs log-linear analysis of contingency tables using RTMB.
#' It models cell frequencies using a Poisson distribution and a log link function.
#' For Bayesian inference, it supports weakly informative priors (Stan-style).
#' For frequentist inference (`classic = TRUE`), it provides Wald-based ANOVA for independence tests.
#'
#' @param formula A formula describing the table structure (e.g., `~ A + B` for independence, `~ A * B` for dependence).
#' @param data A data frame, table, or matrix.
#' @param classic Logical; if TRUE, perform frequentist estimation and ANOVA tests.
#' @param fisher Logical; if TRUE, perform Fisher's exact test (for 2D tables).
#' @param prior Prior specification (e.g., `prior_uniform()` or `prior_weak()`). Default is `prior_uniform()`.
#' @param ... Additional arguments passed to the internal estimation engine.
#'
#' @return An `RTMB_Model`, `MCMC_Fit`, or `Classic_Fit` object depending on the settings.
#'
#' @examples
#' \donttest{
#' # 2x2 table analysis
#' df <- as.data.frame(Titanic)
#' fit <- rtmb_loglinear(Survived ~ Sex + Age, data = df, classic = TRUE)
#' fit$anova()
#'
#' # Bayesian 4-way table interaction analysis
#' fit_bayas <- rtmb_loglinear(~ Class * Sex * Age * Survived, data = Titanic, prior = prior_weak())
#' }
#' @export
rtmb_loglinear <- function(formula, data, prior = prior_uniform(), fixed = NULL, ...) {

  # 1. Data Preparation
  # Convert table or matrix to long data frame
  was_table <- FALSE
  if (inherits(data, c("table", "matrix"))) {
    data <- as.data.frame(as.table(data))
    was_table <- TRUE
  }

  # Identify variables and frequency column
  all_vars <- all.vars(formula)
  response_var <- NULL

  # Check if formula has a response
  if (length(formula) == 3) {
    response_var <- as.character(formula[[2]])
  }

  # If no response is specified
  if (is.null(response_var)) {
    if (was_table && "Freq" %in% names(data)) {
      # For tables/matrices, we automatically use the Freq column
      formula <- update(formula, Freq ~ .)
      response_var <- "Freq"
    } else {
      # For raw data frames, we aggregate
      # Check if variables exist in data
      missing_vars <- setdiff(all_vars, names(data))
      if (length(missing_vars) > 0) {
        stop(sprintf("Variables not found in data: %s", paste(missing_vars, collapse = ", ")))
      }

      # Aggregate raw observations into counts
      data_counts <- as.data.frame(table(data[all_vars]))
      # Update formula to include 'Freq' as response
      formula <- update(formula, Freq ~ .)
      data <- data_counts
      response_var <- "Freq"
    }
  }


  # Ensure response is numeric (counts)
  data[[response_var]] <- as.numeric(data[[response_var]])

  # 2. Invoke rtmb_glmer (or internal logic) with family = "poisson"
  # Since rtmb_glmer is already robust, we leverage it.
  # However, we override some defaults to better suit table analysis.

  # For rtmb_table, we default to prior_weak with Stan-style values if requested
  if (inherits(prior, "prior_weak")) {
    # Stan-style defaults for log-linear (Poisson)
    if (is.null(prior$sd_ratio)) prior$sd_ratio <- 5.0 # For Intercept (alpha_prior_sd)
    if (is.null(prior$max_beta)) prior$max_beta <- 2.5 # For coefficients
  }

  # Call the engine
  # We use rtmb_glmer which handles formulas, random effects, and priors.

  # 2. Determine centering and predictors (to keep generate block clean in print_code)
  args <- list(...)
  is_weak_triggered <- !is.null(args$y_range) && prior$type == "uniform"
  is_centered <- (!prior$type %in% c("flat", "uniform")) || is_weak_triggered

  # Identify if there are predictors
  X_tmp <- stats::model.matrix(nobars(formula), data[1, , drop = FALSE])
  K <- ncol(X_tmp) - (if (attr(stats::terms(nobars(formula)), "intercept")) 1 else 0)

  if (is_centered) {
    int_name <- as.name("Intercept_c")
    x_name <- as.name("X_c")
  } else {
    int_name <- as.name("Intercept")
    x_name <- as.name("X")
  }

  if (K > 0) {
    gen_block <- bquote({
      eta <- .(int_name) + .(x_name) %*% b
      mu <- exp(eta)
      list(mu = mu)
    })
  } else {
    gen_block <- bquote({
      eta <- rep(.(int_name), N)
      mu <- exp(eta)
      list(mu = mu)
    })
  }

  res <- rtmb_glmer(
    formula = formula,
    data = data,
    family = "poisson",
    prior = prior,
    generate = gen_block, fixed = fixed,
    ...
  )

  # Tag as loglinear
  res$type <- "loglinear"
  res$extra$obs_Y <- data[[response_var]]
  fit <- res

  # Tag as rtmb_loglinear for potential custom methods
  class(fit) <- c("rtmb_loglinear", class(fit))

  return(fit)
}

#' RTMB-based GLMM wrapper function
#'
#' @param formula lme4-style formula (e.g., Y ~ X + (1 | GID))
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param laplace Logical; whether to marginalize random effects using Laplace approximation
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Use prior_weak(), prior_rhs(), or prior_ssp(). Default is NULL (flat prior).
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Required when using weakly informative or regularized priors with continuous models.
#' @param init List of initial values (generated automatically based on glm if omitted)
#' @param null Character string specifying the target parameter for the null model.
#' @param gmc Character vector of variable names for Grand Mean Centering (GMC). If "all", all numeric variables are centered.
#' @param cwc List for Centering Within Cluster (CWC). Should contain \code{cluster} (group variable) and \code{pars} (variable names to center).
#' @param view Character vector of parameter names to prioritize in summary.

#' @param factors Character vector of variable names to be treated as factors.
#' @param contrasts Character string specifying the contrast type ("treatment" or "sum").
#' @param sigma_by Character vector specifying variables to group residual variance by (heteroscedasticity).
#' @param resid_corr Residual correlation structure: "ar1" (Autoregressive), "cs" (Compound Symmetry), "toep" (Toeplitz), or "un" (Unstructured).
#' @param resid_time Variable name for time points in residual correlation.
#' @param resid_group Variable name for grouping in residual correlation.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glmer <- function(formula, data, family = "gaussian", laplace = FALSE,
                       prior = prior_uniform(),
                       y_range = NULL,
                       init = NULL,
                       fixed = NULL,
                       null = NULL,
                       gmc = NULL,
                       cwc = NULL,
                       view = NULL,
                       within = NULL,
                       factors = NULL,
                       contrasts = "treatment",
                       sigma_by = NULL,
                       resid_corr = NULL,
                       resid_time = NULL,
                       resid_group = NULL,
                       generate = NULL,
                       .force_sum = FALSE) {

  # --- Handle Wide to Long conversion if needed ---
  if (family == "gaussian" && is.call(formula[[2]]) && identical(formula[[2]][[1]], as.name("cbind"))) {
    processed <- .handle_wide_to_long(formula, data, within, factors)
    formula <- processed$formula
    data <- processed$data
    factors <- processed$factors
  }

  # --- 0. Contrast Management (Automatic sum-to-zero) ---


   if (!is.null(resid_corr)) {
     valid_resid_corr <- c("ar1", "cs", "un", "toep")
     if (!(resid_corr %in% valid_resid_corr)) {
       stop(sprintf("Invalid 'resid_corr' value: '%s'. Valid options are: %s",
            resid_corr, paste(valid_resid_corr, collapse = ", ")))
     }
   }

  # --- 0. Data Preparation (Factors) ---
  if (!is.null(factors)) {
    data <- as.data.frame(data)
    for (f in factors) {
      if (f %in% names(data)) {
        data[[f]] <- as.factor(data[[f]])
      }
    }
  }

  # --- 0. Data Centering (GMC / CWC) ---
  if (!is.null(gmc) || !is.null(cwc)) {
    data_centered <- as.data.frame(data)

    # 1. Grand Mean Centering
    if (!is.null(gmc)) {
      target_gmc <- if (identical(gmc, "all")) {
        names(data_centered)[sapply(data_centered, is.numeric)]
      } else {
        gmc
      }
      for (v in target_gmc) {
        if (v %in% names(data_centered)) {
          data_centered[[v]] <- data_centered[[v]] - mean(data_centered[[v]], na.rm = TRUE)
        } else {
          warning(sprintf("Variable '%s' for GMC not found in data.", v))
        }
      }
    }

    # 2. Centering Within Cluster
    if (!is.null(cwc)) {
      # Robustly extract cluster and pars
      if ("pars" %in% names(cwc)) {
        cluster_var <- cwc$cluster
        target_pars <- cwc$pars
      } else if ("cluster" %in% names(cwc)) {
        cluster_var <- cwc$cluster
        # Any elements not named "cluster" are considered target_pars
        target_pars <- unlist(cwc[names(cwc) != "cluster" | names(cwc) == ""])
      } else if (length(cwc) >= 2) {
        # Assume first is cluster, rest are pars
        cluster_var <- cwc[[1]]
        target_pars <- unlist(cwc[-1])
      } else {
        stop("Invalid 'cwc' format. Use list(cluster = ID, pars = c('var1', 'var2')) or list(ID, 'var1').")
      }

      # Resolve group ID (character or symbol)
      group_id <- if (is.character(cluster_var)) {
        data_centered[[cluster_var]]
      } else {
        # cluster_var might be a symbol (e.g., group)
        # Try to evaluate it in the context of data_centered
        tryCatch(eval(cluster_var, data_centered, parent.frame()),
                 error = function(e) NULL)
      }

      if (is.null(group_id)) {
        # If still NULL, try to deparse and check if it's a column name
        group_var_name <- deparse(cluster_var)
        if (group_var_name %in% names(data_centered)) {
          group_id <- data_centered[[group_var_name]]
        } else {
           stop("Cluster variable for CWC not found.")
        }
      }

      for (v in target_pars) {
        if (v %in% names(data_centered)) {
          group_means <- tapply(data_centered[[v]], group_id, mean, na.rm = TRUE)
          data_centered[[v]] <- data_centered[[v]] - group_means[as.character(group_id)]
        } else {
          warning(sprintf("Variable '%s' for CWC not found in data.", v))
        }
      }
    }
    data <- data_centered
  }

  # --- 0.1 Handle sigma_by (Heteroscedasticity) ---
  sigma_idx <- NULL
  num_sigma_groups <- 1
  sigma_group_names <- "sigma"
  if (!is.null(sigma_by)) {
    if (!family %in% c("gaussian", "lognormal", "student_t")) {
      stop("'sigma_by' is only supported for continuous families (gaussian, lognormal, student_t).")
    } else {
      # Identify categorical factors in sigma_by
      if (identical(sigma_by, "all")) {
        # Extract all variables in fixed part of formula and check if they are factors
        vars_in_form <- all.vars(nobars(formula))
        # Exclude response (the first variable in all.vars is the response)
        vars_in_form <- vars_in_form[-1]
        is_cat <- sapply(data[vars_in_form], function(x) is.factor(x) || is.character(x))
        target_vars <- vars_in_form[is_cat]
      } else {
        target_vars <- sigma_by
      }

      if (length(target_vars) > 0) {
        sigma_factor <- if (length(target_vars) > 1) {
          interaction(data[target_vars], drop = TRUE, sep = ":")
        } else {
          as.factor(data[[target_vars]])
        }
        sigma_idx <- as.integer(sigma_factor)
        num_sigma_groups <- nlevels(sigma_factor)
        sigma_group_names <- levels(sigma_factor)
      }
    }
  }

  # --- 0.2 Default Initial Values ---
  if (is.null(init)) {
    init <- list()
  }
  if (is.list(init)) {
    # Set safe default for residual correlation (e.g. CS, AR1) if not specified by user.
    # This avoids non-positive definite matrices during the initial check.
    if (!is.null(resid_corr) && !("rho_resid" %in% names(init)) && resid_corr %in% c("ar1", "cs", "toep")) {
      init$rho_resid <- 0.1
    }
  }

  if (is.null(prior)) {
    prior <- prior_uniform()
  }

  # Automatically switch to prior_weak() if y_range is provided and prior is default uniform
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && prior$type == "uniform" &&
      is.null(prior$Intercept_sd) && is.null(prior$b_sd) &&
      is.null(prior$sigma_rate) && is.null(prior$tau_rate)) {
    prior <- prior_weak()
  }

  if (!inherits(prior, "rtmb_prior")) {
    stop("prior must be an object of class 'rtmb_prior'. Use prior_weak(), prior_rhs(), or prior_ssp().")
  }

  prior_type <- prior$type
  regularization <- if (prior_type %in% c("rhs", "ssp")) prior_type else "none"
  use_weak_info <- prior_type %in% c("weak", "rhs", "ssp", "jzs")
  has_random <- !is.null(findbars(formula))

  # If classic mode and random effects are present (or forced), we internally
  # force 'sum' contrasts for stable Type III ANOVA/DFs.
  # Otherwise, we use the user-requested contrast (defaulting to 'treatment').
  actual_contrasts <- if (.force_sum) "sum" else contrasts

  if (!is.null(actual_contrasts)) {
    if (actual_contrasts == "sum") {
      old_opts <- options(contrasts = c("contr.sum", "contr.poly"))
      on.exit(options(old_opts), add = TRUE)
    } else if (actual_contrasts == "treatment") {
      old_opts <- options(contrasts = c("contr.treatment", "contr.poly"))
      on.exit(options(old_opts), add = TRUE)
    }
  }

  default_prior <- list(Intercept_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, nu_rate = 0.1, cutpoint_sd = 2.5, shape_rate = 1.0, phi_rate = 1.0, lkj_eta = 1.0)
  prior <- .merge_prior(default_prior, prior)

  if (use_weak_info) {
    if (is.null(prior$max_beta)) prior$max_beta <- 1.0
    if (is.null(prior$sd_ratio)) prior$sd_ratio <- 0.5
    if (is.null(prior$expected_vars)) prior$expected_vars <- 3
    if (is.null(prior$slab_scale)) prior$slab_scale <- 2.0
    if (is.null(prior$slab_df)) prior$slab_df <- 4.0
    if (is.null(prior$ssp_ratio)) prior$ssp_ratio <- 0.25
  }

  families_requiring_yrange <- c("gaussian", "lognormal", "student_t")
  is_continuous <- family %in% families_requiring_yrange

  if (use_weak_info && is_continuous && is.null(y_range) && prior_type != "jzs") {
    stop(sprintf("Specifying 'y_range' is required when using prior_%s() with family = '%s'.", prior_type, family))
  }

  if (has_random) {
    all_vars_form <- subbars(formula)
    fixed_form <- nobars(formula)
    bars <- findbars(formula)

    # Expand bars with * or : using terms() expansion
    expanded_bars <- list()
    if (length(bars) > 0) {
      for (bar in bars) {
        grp_expr <- bar[[3]]
        # Use terms() to expand a*b into a, b, a:b
        term_labels <- attr(terms(as.formula(paste("~", deparse(grp_expr)))), "term.labels")
        for (label in term_labels) {
          new_bar <- bar
          new_bar[[3]] <- parse(text = label)[[1]]
          expanded_bars[[length(expanded_bars) + 1]] <- new_bar
        }
      }
    }
    bars <- expanded_bars
    num_bars <- length(bars)
    suffix <- function(b) if (num_bars > 1) paste0("_", b) else ""

    # Add residual correlation variables to formula for model.frame if needed
    mf_form <- all_vars_form
    if (!is.null(resid_group)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_group)))
    if (!is.null(resid_time)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_time)))

    mf <- model.frame(mf_form, data = data)
    Y <- model.response(mf)
    X <- model.matrix(fixed_form, mf)
    X_assign <- attr(X, "assign")
    X_terms <- attr(terms(fixed_form), "term.labels")
    offset <- model.offset(mf)
    N <- nrow(mf)

    dat_ranef <- list()
    num_ranef_list <- list()
    num_groups_list <- list()
    ranef_names_list <- list()
    group_labels_list <- list()

    for (b in 1:num_bars) {
      bar <- bars[[b]]
      re_form <- as.formula(paste("~", deparse(bar[[2]])))

      Z_mf <- model.matrix(re_form, mf)

      # Handle grouping factor (potentially with interactions)
      group_var_label <- deparse(bar[[3]])
      vars <- trimws(unlist(strsplit(group_var_label, ":")))

      if (length(vars) > 1) {
        flist <- interaction(mf[vars], sep = ":", drop = TRUE)
      } else {
        flist <- as.factor(mf[[vars]])
      }
      group_idx <- as.integer(flist)

      if (length(group_idx) == 0) {
        stop("Grouping factor '", group_var_label, "' resulted in an empty vector.")
      }
      if (any(is.na(group_idx))) {
        stop("NA values found in grouping factor '", group_var_label, "'. Check your data.")
      }

      num_groups <- length(levels(flist))
      num_ranef <- ncol(Z_mf)

      Zt <- matrix(0, nrow = num_groups * num_ranef, ncol = N)
      for (i in 1:N) {
        g <- group_idx[i]
        Zt[((g - 1) * num_ranef + 1):(g * num_ranef), i] <- Z_mf[i, ]
      }

      dat_ranef[[paste0("Zt", suffix(b))]] <- Zt
      dat_ranef[[paste0("group_idx", suffix(b))]] <- group_idx

      num_ranef_list[[b]] <- num_ranef
      num_groups_list[[b]] <- num_groups

      r_names <- colnames(Z_mf)
      r_names[r_names == "(Intercept)"] <- "Int"
      ranef_names_list[[b]] <- r_names
      group_labels_list[[b]] <- group_var_label
    }
  } else {
    # Add residual correlation variables to formula for model.frame if needed
    mf_form <- formula
    if (!is.null(resid_group)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_group)))
    if (!is.null(resid_time)) mf_form <- update(mf_form, as.formula(paste("~ . +", resid_time)))

    mf <- model.frame(mf_form, data)
    Y <- model.response(mf)
    X <- model.matrix(formula, mf)
    X_assign <- attr(X, "assign")
    X_terms <- attr(terms(formula), "term.labels")
    offset <- model.offset(mf)
    N <- nrow(mf)
    num_bars <- 0
    dat_ranef <- list()
    num_ranef_list <- list()
    num_groups_list <- list()
    ranef_names_list <- list()
    group_labels_list <- list()
  }
  
  jzs_V_val <- NULL

  has_intercept <- "(Intercept)" %in% colnames(X)
  if (family == "ordered") has_intercept <- FALSE
  fixed_colnames <- colnames(X)
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  fixed_names <- colnames(X)
  K_tmp <- ncol(X)

  if (is.matrix(Y)) {
    if (ncol(Y) == 2 && family == "binomial") {
      trials <- as.numeric(Y[, 1] + Y[, 2])
      Y <- as.numeric(Y[, 1])
    } else if (ncol(Y) == 1) {
      Y <- as.numeric(Y[, 1])
      trials <- rep(1, length(Y))
    } else stop("Invalid matrix format for the response variable.")
  } else {
    trials <- rep(1, length(Y))
    if (is.factor(Y)) {
      Y <- if (family %in% c("bernoulli", "binomial")) as.numeric(Y) - 1 else as.numeric(Y)
    } else Y <- as.numeric(Y)
  }

  # --- Automatic generation of initial values for SSP ---
  if (regularization == "ssp" && is.null(init) && K_tmp > 0) {
    rhs_mod <- NULL
    suppressMessages(capture.output({
      rhs_prior_obj <- prior_rhs(expected_vars = prior$expected_vars, slab_scale = prior$slab_scale, slab_df = prior$slab_df)
      rhs_mod <- rtmb_glmer(formula, data, family, laplace, prior = rhs_prior_obj, y_range = y_range, init = NULL)
    }))

    jac_target <- if (laplace) "random" else "none"

    ad_setup <- NULL
    suppressMessages(capture.output({
      ad_setup <- rhs_mod$build_ad_obj(laplace = laplace, jacobian_target = jac_target)
    }))

    ad_obj <- ad_setup$ad_obj
    best_obj <- Inf
    best_par <- ad_obj$par

    for (i in 1:4) {
      start_par <- if (i == 1) ad_obj$par else ad_obj$par + rnorm(length(ad_obj$par), 0, 0.5)
      opt <- suppressWarnings(try(nlminb(start = start_par, objective = ad_obj$fn, gradient = ad_obj$gr), silent = TRUE))
      if (!inherits(opt, "try-error") && !is.na(opt$objective) && opt$objective < best_obj) {
        best_obj <- opt$objective
        best_par <- opt$par
      }
    }

    ad_obj$fn(best_par)
    full_par <- if (!is.null(ad_obj$env$last.par.best)) ad_obj$env$last.par.best else ad_obj$env$last.par

    unc_est_list <- unconstrained_vector_to_list(full_par, rhs_mod$par_list)
    con_est_list <- to_constrained(unc_est_list, rhs_mod$par_list)
    if (!is.null(rhs_mod$transform)) {
      tran_res <- rhs_mod$transform(rhs_mod$data, con_est_list)
      con_est_list <- c(con_est_list, tran_res)
    }

    b_est <- as.numeric(con_est_list$b)
    if (length(b_est) == 0 || any(is.na(b_est))) b_est <- rep(0, K_tmp)

    r_init <- ifelse(abs(b_est) > mean(abs(b_est)), 0.9, 0.1)
    tau_init <- pmax(abs(b_est), 0.01)
    beta_raw_init <- b_est / (r_init * tau_init)

    init <- list(beta_raw = beta_raw_init, r = r_init, tau = tau_init)
    if (has_intercept && !is.null(con_est_list$Intercept_c)) {
      init$Intercept_c <- as.numeric(con_est_list$Intercept_c)[1]
    }
  } else if (is.null(init)) {
    tryCatch({
      fixed_formula <- nobars(formula)
      glm_fam <- switch(family, "gaussian" = gaussian(), "student_t" = gaussian(), "lognormal" = gaussian(link="log"), "bernoulli" = binomial(), "binomial" = binomial(), "poisson" = poisson(), "neg_binomial" = poisson(), "gamma" = Gamma(link="log"), gaussian())
      init_fit <- suppressWarnings(glm(fixed_formula, data = data, family = glm_fam))
      init_coef <- unname(coef(init_fit))
      init_coef[is.na(init_coef)] <- 0

      init <- list()
      if (has_intercept) {
        if (prior_type %in% c("flat", "uniform")) {
          init$Intercept <- init_coef[1]
        } else {
          init$Intercept_c <- init_coef[1]
        }
      }

      if (regularization == "none" && K_tmp > 0) {
        if (has_intercept) init$b <- init_coef[-1]
        else init$b <- init_coef
      }

      # Initial values for residual correlation
      if (!is.null(resid_corr)) {
        if (resid_corr %in% c("ar1", "cs")) init$rho_resid <- 0.1
        else if (resid_corr == "toep") init$rho_resid <- rep(0.1, dat$max_T_resid - 1)
        else if (resid_corr == "un") init$L_resid <- diag(dat$max_T_resid)
      }

      if (family %in% c("gaussian", "lognormal", "student_t")) {
        init$sigma <- if (num_sigma_groups > 1) rep(stats::sd(Y) * 0.5, num_sigma_groups) else stats::sd(Y) * 0.5
      }

    }, error = function(e) init <- NULL)
  }

  dat <- list(Y = Y, trials = trials, X = X,
              sigma_idx = sigma_idx, num_sigma_groups = num_sigma_groups)

  if (!is.null(resid_corr)) {
    if (family != "gaussian") stop("Residual correlation structures are currently only supported for Gaussian models.")
    resid_corr <- tolower(resid_corr)
    if (!resid_corr %in% c("ar1", "cs", "un", "toep")) stop("Currently 'ar1', 'cs', 'un', and 'toep' are supported for resid_corr.")

    # Identify group for residuals
    if (is.null(resid_group)) {
      if (has_random) {
        resid_group_var <- as.character(bars[[1]][[3]])
      } else {
        stop("Residual correlation requires 'resid_group' (e.g., resid_group = 'ID') if no random effects are specified.")
      }
    } else {
      resid_group_var <- resid_group
    }

    if (!(resid_group_var %in% names(mf))) stop(sprintf("Residual grouping variable '%s' not found in data.", resid_group_var))

    group_resid <- as.numeric(as.factor(mf[[resid_group_var]]))
    dat$group_resid <- group_resid
    dat$num_groups_resid <- max(group_resid)
    if (!is.null(resid_time)) {
      dat$time_resid <- as.numeric(as.factor(mf[[resid_time]])) - 1
    } else {
      # Use row order within group
      dat$time_resid <- ave(rep(1, N), group_resid, FUN = seq_along) - 1
    }
    dat$max_T_resid <- max(table(group_resid))
  }

  if (has_random) {
    dat <- c(dat, dat_ranef)
  }
  if (!is.null(sigma_idx)) {
    dat$sigma_idx <- sigma_idx
    dat$num_sigma_groups <- num_sigma_groups
  }
  if (!is.null(offset)) dat$offset <- offset
  if (family == "ordered") dat$num_categories <- length(unique(Y))

  # --- Setup AST ---
  setup_exprs <- list()
  setup_exprs[[1]] <- quote(N <- length(Y))
  setup_exprs[[2]] <- quote(K <- ncol(X))
  if (has_random) {
    for (b in 1:num_bars) {
      group_idx_name <- as.name(paste0("group_idx", suffix(b)))
      Zt_name <- as.name(paste0("Zt", suffix(b)))
      Z_mat_name <- as.name(paste0("Z_mat", suffix(b)))
      num_groups_name <- as.name(paste0("num_groups", suffix(b)))
      num_ranef_name <- as.name(paste0("num_ranef", suffix(b)))

      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(num_groups_name) <- max(.(group_idx_name)))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(num_ranef_name) <- nrow(.(Zt_name)) / .(num_groups_name))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(Z_mat_name) <- matrix(0, nrow = N, ncol = .(num_ranef_name)))

      loop_body <- bquote({
        g <- .(group_idx_name)[i]
        .(Z_mat_name)[i, ] <- .(Zt_name)[((g - 1) * .(num_ranef_name) + 1):(g * .(num_ranef_name)), i]
      })
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(for (i in 1:N) .(loop_body))
    }
  }

  if (use_weak_info) {
    if (family %in% c("bernoulli", "binomial", "ordered")) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- pi / sqrt(3))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(alpha_prior_sd <- base_scale * .(prior$max_beta))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- 0)
    } else if (family %in% c("poisson", "neg_binomial", "gamma")) {
      # For Poisson/log-link models, Stan often uses Normal(log(mean(y)), 2.5 or 5) for intercept
      # and Normal(0, 2.5) for coefficients.
      setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- log(mean(Y) + 0.5))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(alpha_prior_sd <- .(prior$sd_ratio))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- 1.0)
    } else {
      if (prior_type == "jzs") {
         # For JZS, we use defaults if y_range is missing
         if (!is.null(y_range)) {
           setup_exprs[[length(setup_exprs) + 1]] <- quote(half_d_y <- diff(y_range) / 2)
           setup_exprs[[length(setup_exprs) + 1]] <- bquote(base_scale <- half_d_y * .(prior$sd_ratio))
           setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- half_d_y)
           setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mean(y_range))
         } else {
           # JZS defaults for alpha/sigma part if y_range is missing (though we use flat priors anyway)
           setup_exprs[[length(setup_exprs) + 1]] <- quote(base_scale <- 1.0)
           setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- 10.0)
           setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- 0)
         }
         setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate <- 1.0 / base_scale)
         setup_exprs[[length(setup_exprs) + 1]] <- bquote(jzs_r <- .(prior$r))
      } else {
         setup_exprs[[length(setup_exprs) + 1]] <- quote(half_d_y <- diff(y_range) / 2)
         setup_exprs[[length(setup_exprs) + 1]] <- bquote(base_scale <- half_d_y * .(prior$sd_ratio))
         setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- half_d_y)
         setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mean(y_range))
         setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate <- 1.0 / base_scale)
      }
    }

    setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_rate <- 1.0 / base_scale)

    if (K_tmp > 0) {
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_mean <- apply(X, 2, mean))
      setup_exprs[[length(setup_exprs) + 1]] <- quote(X_sd <- apply(X, 2, sd))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(beta_prior_sd <- .(prior$max_beta) * base_scale / X_sd)
      if (has_intercept && !prior_type %in% c("flat", "uniform")) {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(X_c <- X - rep(1, N) %*% t(X_mean))
      } else if (prior_type == "jzs") {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(X_c <- X - rep(1, N) %*% t(X_mean))
      }

      if (prior_type == "jzs") {
        # Calculate (X_c' X_c / N)^-1 for multivariate JZS
        # X_c is centered X.
        X_mean_val <- apply(X, 2, mean)
        X_c_val <- X - rep(1, N) %*% t(X_mean_val)
        
        X_for_jzs <- X_c_val
        
        if (ncol(X_for_jzs) > 0) {
          # Calculate (X'X/N)^-1 with numerical safety
          XtX_inv <- solve(t(X_for_jzs) %*% X_for_jzs / N)
          # Force symmetry and add a tiny jitter to ensure PD
          jzs_V_val <- 0.5 * (XtX_inv + t(XtX_inv))
          diag(jzs_V_val) <- diag(jzs_V_val) + 1e-10
          
          jzs_scales <- sqrt(diag(jzs_V_val))
          
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(jzs_V <- .(jzs_V_val))
        } else {
          # Null model case: no predictors for JZS
          setup_exprs[[length(setup_exprs) + 1]] <- quote(jzs_V <- matrix(0, 0, 0))
        }
      }

      if (regularization == "rhs") {
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(p0 <- min(.(prior$expected_vars), K - 1))
        setup_exprs[[length(setup_exprs) + 1]] <- quote(if (p0 < 1) p0 <- 1)
        setup_exprs[[length(setup_exprs) + 1]] <- quote(tau0 <- (p0 / (K - p0)) * (base_scale / sqrt(N)))
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(half_slab_df <- .(prior$slab_df) / 2.0)
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(half_slab_scale2 <- (.(prior$slab_df) * .(prior$slab_scale)^2) / 2.0)
      } else if (regularization == "ssp") {
        setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_scale <- base_scale / X_sd)
      }
    }
  } else {
    # Minimal calculations for fixed prior distributions
    # No centering or scaling for flat/uniform priors
    NULL
  }
  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(dat)
  eval(setup_ast, tmp_env)
  N <- tmp_env$N; K <- tmp_env$K
  num_sigma_groups <- if (!is.null(tmp_env$num_sigma_groups)) tmp_env$num_sigma_groups else 1

  # --- Parameters AST ---
  param_exprs <- list()
  if (has_intercept) {
    if (prior_type %in% c("flat", "uniform")) {
      param_exprs[[length(param_exprs) + 1]] <- quote(Intercept <- Dim(1))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- quote(Intercept_c <- Dim(1))
    }
  }
  if (K > 0) {
    if (prior_type == "jzs") {
      param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(K))
    } else if (regularization == "none") {
      param_exprs[[length(param_exprs) + 1]] <- quote(b <- Dim(K))
    } else if (regularization == "rhs") {
      param_exprs[[length(param_exprs) + 1]] <- quote(z <- Dim(K))
      param_exprs[[length(param_exprs) + 1]] <- quote(lambda <- Dim(K, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(w_lambda <- Dim(K, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(tau_hs <- Dim(1, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(w_tau <- Dim(1, lower = 0))
      param_exprs[[length(param_exprs) + 1]] <- quote(c2 <- Dim(1, lower = 0))
    } else if (regularization == "ssp") {
      param_exprs[[length(param_exprs) + 1]] <- quote(beta_raw <- Dim(K))
      param_exprs[[length(param_exprs) + 1]] <- quote(r <- Dim(K, lower = 0.001, upper = 0.999))
      param_exprs[[length(param_exprs) + 1]] <- quote(tau <- Dim(K, lower = 0))
    }
  }

  if (has_random) {
    for (b in 1:num_bars) {
      num_ranef_val <- num_ranef_list[[b]]
      num_groups_val <- num_groups_list[[b]]
      sd_name <- as.name(paste0("sd", suffix(b)))
      r_re_name <- as.name(paste0("r_re", suffix(b)))
      CF_corr_name <- as.name(paste0("CF_corr", suffix(b)))

      if (num_ranef_val == 1) {
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(sd_name) <- Dim(.(num_ranef_val), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(r_re_name) <- Dim(.(num_groups_val), random = TRUE))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(sd_name) <- Dim(.(num_ranef_val), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(CF_corr_name) <- Dim(c(.(num_ranef_val), .(num_ranef_val)), type = "CF_corr"))
        param_exprs[[length(param_exprs) + 1]] <- bquote(.(r_re_name) <- Dim(c(.(num_groups_val), .(num_ranef_val)), random = TRUE))
      }
    }
  }

  if (family %in% c("gaussian", "lognormal", "student_t")) param_exprs[[length(param_exprs) + 1]] <- quote(sigma <- Dim(num_sigma_groups, lower = 0))
  if (!is.null(resid_corr)) {
    if (resid_corr %in% c("ar1", "cs")) {
      param_exprs[[length(param_exprs) + 1]] <- quote(rho_resid <- Dim(type = "interval", lower = -0.99, upper = 0.99))
    } else if (resid_corr == "un") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(L_resid <- Dim(c(.(dat$max_T_resid), .(dat$max_T_resid)), type = "CF_corr"))
    } else if (resid_corr == "toep") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(rho_resid <- Dim(.(dat$max_T_resid - 1), type = "interval", lower = -0.99, upper = 0.99))
    }
  }
  if (family == "student_t") param_exprs[[length(param_exprs) + 1]] <- quote(nu <- Dim(1, lower = 2))
  if (family == "gamma") param_exprs[[length(param_exprs) + 1]] <- quote(shape <- Dim(1, lower = 0))
  if (family == "neg_binomial") param_exprs[[length(param_exprs) + 1]] <- quote(phi <- Dim(1, lower = 0))
  if (family == "ordered") param_exprs[[length(param_exprs) + 1]] <- quote(cutpoints <- Dim(num_categories - 1, type = "ordered"))
  param_ast <- as.call(c(list(as.name("{")), param_exprs))

  # --- Transform AST ---
  tran_exprs <- list()
  if (K > 0) {
    if (regularization == "rhs") {
      tran_exprs[[length(tran_exprs) + 1]] <- quote(lambda_sq <- lambda^2)
      tran_exprs[[length(tran_exprs) + 1]] <- quote(tau_sq <- tau_hs^2)
      tran_exprs[[length(tran_exprs) + 1]] <- quote(lambda_tilde <- sqrt((c2 * lambda_sq) / (c2 + tau_sq * lambda_sq)))
      tran_exprs[[length(tran_exprs) + 1]] <- quote(b <- z * lambda_tilde * tau_hs)
    } else if (regularization == "ssp") {
      tran_exprs[[length(tran_exprs) + 1]] <- quote(b <- beta_raw * r * tau)
    }
  }
  if (has_intercept) {
    if (!prior_type %in% c("flat", "uniform")) {
      if (K > 0) tran_exprs[[length(tran_exprs) + 1]] <- quote(Intercept <- Intercept_c - sum(X_mean * b))
      else tran_exprs[[length(tran_exprs) + 1]] <- quote(Intercept <- Intercept_c)
    }
  }
  if (has_random) {
    for (b in 1:num_bars) {
      if (num_ranef_list[[b]] > 1) {
        CF_corr_name <- as.name(paste0("CF_corr", suffix(b)))
        corr_name <- as.name(paste0("corr", suffix(b)))
        tran_exprs[[length(tran_exprs) + 1]] <- bquote(.(corr_name) <- .(CF_corr_name) %*% t(.(CF_corr_name)))
      }
    }
  }
  tran_ast <- if (length(tran_exprs) > 0) as.call(c(list(as.name("{")), tran_exprs)) else NULL

  # --- Model AST ---
  transform_exprs <- list()
  if (has_intercept) {
    if (prior_type %in% c("flat", "uniform")) {
      if (K > 0) transform_exprs[[1]] <- quote(eta <- Intercept + X %*% b)
      else transform_exprs[[1]] <- quote(eta <- rep(Intercept, N))
    } else {
      if (K > 0) transform_exprs[[1]] <- quote(eta <- Intercept_c + X_c %*% b)
      else transform_exprs[[1]] <- quote(eta <- rep(Intercept_c, N))
    }
  } else {
    if (K > 0) transform_exprs[[1]] <- quote(eta <- X %*% b)
    else transform_exprs[[1]] <- quote(eta <- rep(0, N))
  }

  if (has_random) {
    for (b in 1:num_bars) {
      Z_mat_name <- as.name(paste0("Z_mat", suffix(b)))
      r_re_name <- as.name(paste0("r_re", suffix(b)))
      group_idx_name <- as.name(paste0("group_idx", suffix(b)))
      sd_name <- as.name(paste0("sd", suffix(b)))

      if (num_ranef_list[[b]] > 1) {
        transform_exprs[[length(transform_exprs) + 1]] <- bquote(for (i in 1:N) eta[i] <- eta[i] + sum(.(Z_mat_name)[i, ] * .(r_re_name)[.(group_idx_name)[i], ] * .(sd_name)))
      } else {
        transform_exprs[[length(transform_exprs) + 1]] <- bquote(eta <- eta + .(Z_mat_name)[,1] * .(r_re_name)[.(group_idx_name)] * .(sd_name))
      }
    }
  }
  if (!is.null(offset)) transform_exprs[[length(transform_exprs) + 1]] <- quote(eta <- eta + offset)

  ll_data_exprs <- list()
  if (!is.null(resid_corr)) {
    vi_inner_expr <- switch(resid_corr,
      "ar1" = quote(Vi[r,c] <- s_r * s_c * (rho_resid^abs(ti[r] - ti[c]))),
      "cs"  = quote(if (r == c) Vi[r,c] <- s_r^2 else Vi[r,c] <- s_r * s_c * rho_resid),
      "toep"= quote({
               lag <- abs(ti[r] - ti[c])
               if (lag == 0) Vi[r,c] <- s_r * s_c
               else if (lag <= length(rho_resid)) Vi[r,c] <- s_r * s_c * rho_resid[lag]
               else Vi[r,c] <- 0
            }),
      "un" = quote({
               if (r == c) Vi[r,c] <- s_r^2
               else {
                 idx_r <- ti[r] + 1
                 idx_c <- ti[c] + 1
                 Vi[r,c] <- s_r * s_c * R_mat_resid[idx_r, idx_c]
               }
            })
    )

    # Build model components list
    ll_data_exprs <- list()
    ll_data_exprs[[length(ll_data_exprs) + 1]] <- quote(has_sig_idx <- !is.null(sigma_idx))
    if (resid_corr == "un") {
      ll_data_exprs[[length(ll_data_exprs) + 1]] <- quote(R_mat_resid <- L_resid %*% t(L_resid))
    }

    ll_data_exprs[[length(ll_data_exprs) + 1]] <- bquote(
      for (i in 1:num_groups_resid) {
        idx <- which(group_resid == i)
        ni <- length(idx)
        ti <- time_resid[idx]
        mu_i <- eta[idx]

        Vi <- matrix(0, ni, ni)
        for (r in 1:ni) {
          s_r <- if (has_sig_idx) sigma[sigma_idx[idx[r]]] else sigma
          for (c in 1:ni) {
            s_c <- if (has_sig_idx) sigma[sigma_idx[idx[c]]] else sigma
            .(vi_inner_expr)
          }
        }
        diag(Vi) <- diag(Vi) + 1e-3
        Y[idx] ~ multi_normal(mu_i, Vi)
      }
    )
  } else {
    ll_data_exprs[[1]] <- switch(family,
                                 "gaussian" = if (!is.null(sigma_idx)) quote(Y ~ normal(eta, sigma[sigma_idx])) else quote(Y ~ normal(eta, sigma)),
                                 "lognormal" = if (!is.null(sigma_idx)) quote(Y ~ lognormal(eta, sigma[sigma_idx])) else quote(Y ~ lognormal(eta, sigma)),
                                 "student_t" = if (!is.null(sigma_idx)) quote(Y ~ student_t(nu, eta, sigma[sigma_idx])) else quote(Y ~ student_t(nu, eta, sigma)),
                                 "gamma" = quote(Y ~ gamma(shape, shape / exp(eta))),
                                 "bernoulli" = quote(Y ~ bernoulli_logit(eta)),
                                 "binomial" = quote(Y ~ binomial_logit(trials, eta)),
                                 "poisson" = quote(Y ~ poisson(exp(eta))),
                                 "neg_binomial" = quote(Y ~ neg_binomial_2(exp(eta), phi)),
                                 "ordered" = quote(Y ~ ordered_logistic(eta, cutpoints))
    )
  }

  ll_random_exprs <- list()
  if (has_random) {
    for (b in 1:num_bars) {
      r_re_name <- as.name(paste0("r_re", suffix(b)))
      if (num_ranef_list[[b]] > 1) {
        CF_corr_name <- as.name(paste0("CF_corr", suffix(b)))
        num_ranef_name <- as.name(paste0("num_ranef", suffix(b)))
        num_groups_name <- as.name(paste0("num_groups", suffix(b)))

        ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(.(CF_corr_name) ~ lkj_CF_corr(.(prior$lkj_eta)))
        ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(for (j in 1:.(num_groups_name)) .(r_re_name)[j, ] ~ multi_normal_CF(rep(0, .(num_ranef_name)), rep(1, .(num_ranef_name)), .(CF_corr_name)))
      } else {
        ll_random_exprs[[length(ll_random_exprs) + 1]] <- bquote(.(r_re_name) ~ normal(0, 1))
      }
    }
  }

  prior_exprs <- list()
  if (prior_type %in% c("weak", "rhs", "ssp", "jzs")) {
    if (family %in% c("gaussian", "lognormal", "student_t")) {
      if (use_weak_info && prior_type == "weak") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(sigma ~ exponential(sigma_rate))
      } else if (prior_type == "jzs") {
        # Jeffreys prior for variance (p(sigma) proportional to 1/sigma)
        if (num_sigma_groups > 1) {
          prior_exprs[[length(prior_exprs) + 1]] <- quote(lp <- lp - sum(log(sigma)))
        } else {
          prior_exprs[[length(prior_exprs) + 1]] <- quote(lp <- lp - log(sigma))
        }
      } else if (!is.null(prior$sigma_rate)) {
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
      }
    }

    if (family == "student_t" && !is.null(prior$nu_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(nu ~ exponential(.(prior$nu_rate)))
    if (family == "gamma" && !is.null(prior$shape_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(shape ~ exponential(.(prior$shape_rate)))
    if (family == "neg_binomial" && !is.null(prior$phi_rate)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(phi ~ exponential(.(prior$phi_rate)))
    if (family == "ordered" && !is.null(prior$cutpoint_sd)) prior_exprs[[length(prior_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))

    if (has_intercept) {
      if (use_weak_info && prior_type == "weak") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(Intercept_c ~ normal(mid_y, alpha_prior_sd))
      } else if (prior_type == "jzs") {
        # Flat prior for intercept (p(alpha) proportional to 1) - no term added
        prior_exprs[[length(prior_exprs) + 1]] <- quote("# Flat prior for Intercept (JZS)")
      } else if (!is.null(prior$Intercept_sd)) {
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(Intercept_c ~ normal(0, .(prior$Intercept_sd)))
      }
    }

    if (K > 0) {
      if (prior_type == "jzs") {
        # Multivariate JZS (Single-g prior): b ~ Multivariate-Cauchy(0, r^2 * sigma^2 * (X'X/N)^-1)
        # This matches the BayesFactor package and handles correlations correctly.
        # We use bquote to embed the pre-calculated dimension.
        if (!is.null(jzs_V_val) && ncol(jzs_V_val) > 0) {
          K_jzs <- ncol(jzs_V_val)
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(b ~ multi_cauchy(rep(0, .(K_jzs)), (jzs_r * sigma)^2 * jzs_V))
        }
      } else if (regularization == "none") {
        if (use_weak_info && prior_type == "weak") {
          prior_exprs[[length(prior_exprs) + 1]] <- quote(b ~ normal(0, beta_prior_sd))
        } else if (!is.null(prior$b_sd)) {
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
        }
      } else if (regularization == "rhs") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(z ~ normal(0, 1))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(w_lambda ~ gamma(0.5, 0.5))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(lambda ~ half_normal(1 / sqrt(w_lambda)))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(w_tau ~ gamma(0.5, 0.5))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(tau_hs ~ half_normal(tau0 / sqrt(w_tau)))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(c2 ~ inverse_gamma(half_slab_df, half_slab_scale2))
      } else if (regularization == "ssp") {
        prior_exprs[[length(prior_exprs) + 1]] <- quote(beta_raw ~ laplace(0, 0.5))
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(mu_r <- log(.(prior$ssp_ratio) / (1 - .(prior$ssp_ratio))))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(r ~ logit_normal(mu_r, 3))
        prior_exprs[[length(prior_exprs) + 1]] <- quote(tau ~ exponential(1 / tau_scale))
      }
    }

    if (has_random) {
      for (b in 1:num_bars) {
        sd_name <- as.name(paste0("sd", suffix(b)))
        if (use_weak_info && prior_type == "weak") {
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(.(sd_name) ~ exponential(tau_rate))
        } else if (!is.null(prior$sd_rate)) {
          prior_exprs[[length(prior_exprs) + 1]] <- bquote(.(sd_name) ~ exponential(.(prior$sd_rate)))
        }
      }
    }
  } else if (prior_type == "uniform") {
    if (!is.null(prior$sigma_rate)) {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    }
    if (has_intercept && !is.null(prior$Intercept_sd)) {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(Intercept ~ normal(0, .(prior$Intercept_sd)))
    }
    if (K > 0 && !is.null(prior$b_sd)) {
      prior_exprs[[length(prior_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
    }
    if (has_random && !is.null(prior$tau_rate)) {
      for (b in 1:num_bars) {
        sd_name <- as.name(paste0("sd", suffix(b)))
        prior_exprs[[length(prior_exprs) + 1]] <- bquote(.(sd_name) ~ exponential(.(prior$tau_rate)))
      }
    }
  }

  model_exprs <- list()
  model_exprs[[length(model_exprs) + 1]] <- "# Transform"
  model_exprs <- c(model_exprs, transform_exprs)
  model_exprs[[length(model_exprs) + 1]] <- "# Likelihood (Data)"
  model_exprs <- c(model_exprs, ll_data_exprs)
  if (has_random) {
    model_exprs[[length(model_exprs) + 1]] <- "# Likelihood (Random)"
    model_exprs <- c(model_exprs, ll_random_exprs)
  }
  model_exprs[[length(model_exprs) + 1]] <- "# Priors"
  model_exprs <- c(model_exprs, prior_exprs)

  model_ast <- as.call(c(list(as.name("{")), model_exprs))

  code_obj <- list(setup = setup_ast, parameters = param_ast)
  if (!is.null(tran_ast)) code_obj$transform <- tran_ast
  code_obj$model <- model_ast
  if (!is.null(generate)) code_obj$generate <- generate

  par_names_list <- list()
  if (num_sigma_groups > 1) {
    par_names_list$sigma <- sigma_group_names
  }
  if (K > 0) {
    par_names_list$b <- fixed_names
    if (regularization == "rhs") {
      par_names_list$z <- fixed_names; par_names_list$lambda <- fixed_names; par_names_list$w_lambda <- fixed_names
    } else if (regularization == "ssp") {
      par_names_list$beta_raw <- fixed_names; par_names_list$r <- fixed_names; par_names_list$tau <- fixed_names
    }
  }
  if (has_random) {
    for (b in 1:num_bars) {
      r_names <- ranef_names_list[[b]]
      if (num_bars > 1) {
        r_names <- paste0(r_names, ":", group_labels_list[[b]])
      }
      par_names_list[[paste0("sd", suffix(b))]] <- paste0(group_labels_list[[b]], ":", r_names)
      if (num_ranef_list[[b]] > 1) par_names_list[[paste0("corr", suffix(b))]] <- r_names
    }
  }

  view_vars <- c()
  if (has_intercept) view_vars <- c("Intercept")
  if (K > 0) view_vars <- c(view_vars, "b")
  view_vars <- c(view_vars, "sigma")
  if (has_random) {
    for (b in 1:num_bars) {
      view_vars <- c(view_vars, paste0("sd", suffix(b)))
      if (num_ranef_list[[b]] > 1) view_vars <- c(view_vars, paste0("corr", suffix(b)))
    }
  }
  if (!is.null(resid_corr)) {
    if (resid_corr %in% c("ar1", "cs", "toep")) view_vars <- c(view_vars, "rho_resid")
    if (resid_corr == "un") view_vars <- c(view_vars, "L_resid")
  }

  ordered_data <- env_to_ordered_list(tmp_env, dat, setup_ast)
  obj <- rtmb_model(data = ordered_data, code = code_obj, par_names = par_names_list, init = init,
                    fixed = fixed,
                    view = if (!is.null(view)) view else view_vars, silent = FALSE)
  obj$formula <- formula
  obj$raw_data <- data
  obj$family <- family
  # Store preferences
  obj$contrasts <- actual_contrasts
  obj$requested_contrasts <- contrasts

  if (!is.null(null)) {
    obj <- obj$null_model(pars = null)
  }

  # --- Metadata for classic inference and other methods ---
  obj$type <- if (has_random || !is.null(resid_corr)) {
    if (family %in% c("gaussian", "lognormal", "student_t")) "lmer" else "glmm"
  } else {
    if (family %in% c("gaussian", "lognormal", "student_t")) "lm" else "glm"
  }
  obj$extra <- list(
    X_assign = X_assign, 
    X_terms = X_terms, 
    X_colnames = fixed_colnames,
    factors = factors,
    within = within
  )

  fixed_effects <- if (K > 0) "b" else character(0)
  if (has_intercept) {
    fixed_effects <- c(ifelse(prior_type %in% c("flat", "uniform"), "Intercept", "Intercept_c"), fixed_effects)
  }
  obj$extra$df_pars <- fixed_effects

  return(obj)
}

#' RTMB-based Linear Mixed Model (LMM) wrapper function
#'
#' @param formula Formula
#' @param data Data frame
#' @param laplace Logical; whether to marginalize random effects using Laplace approximation
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Default is NULL (flat prior).
#' @param y_range Theoretical minimum and maximum values of the response variable
#' @param init Initial values
#' @param null Null model parameters
#' @param gmc Character vector of variable names for GMC
#' @param cwc List for CWC
#' @param view Character vector of parameter names to prioritize in summary.
#' @param factors Character vector of variable names to be treated as factors.
#' @param contrasts Character string specifying the contrast type ("treatment" or "sum").
#' @param sigma_by Character vector specifying variables to group residual variance by (heteroscedasticity).
#' @param resid_corr Residual correlation structure (e.g., "ar1", "cs", "un", "toep").
#' @param resid_time Variable name for time points in residual correlation.
#' @param resid_group Variable name for grouping in residual correlation.
#' @param within Optional list for wide-to-long conversion. For repeated measures data in wide format,
#' specify the factor names and their levels, e.g., \code{list(Time = 4)} or \code{list(A = 2, B = 3)}.
#' The total number of levels must match the number of columns in \code{cbind()} on the LHS.
#' If omitted and the LHS is \code{cbind()}, the within-factor name is inferred from RHS variables not present in the data.
#' @param ... Additional arguments passed to \code{rtmb_model()}.
#'
#' @return RTMB_Model object
#' @export
#' @example inst/examples/ex_lm.R
rtmb_lmer <- function(formula, data, laplace = TRUE,
                       prior = prior_uniform(),
                       y_range = NULL,
                       init = NULL,
                       fixed = NULL,
                       null = NULL,
                       gmc = NULL,
                       cwc = NULL,
                       view = NULL,
                       sigma_by = NULL,
                       factors = NULL,
                       contrasts = "treatment",
                       resid_corr = NULL,
                       resid_time = NULL,
                       resid_group = NULL,
                       within = NULL,
                       ...) {
  rtmb_glmer(formula = formula, data = data, family = "gaussian",
             laplace = laplace,
             prior = prior,
             y_range = y_range,
             init = init,
             null = null,
             gmc = gmc,
             cwc = cwc,
             view = view,
             sigma_by = sigma_by,
             factors = factors,
             contrasts = contrasts,
             resid_corr = resid_corr,
             resid_time = resid_time,
             resid_group = resid_group,
             within = within,
             .force_sum = TRUE)
}

#' RTMB-based GLM wrapper function (no random effects)
#'
#' @param formula Formula
#' @param data Data frame
#' @param family Character string of the distribution family (e.g., "gaussian", "binomial", "poisson")
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Use prior_weak(), prior_rhs(), or prior_ssp(). Default is NULL (flat prior).
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param init List of initial values
#' @param null Character string specifying the target parameter for the null model.
#' @param gmc Character vector of variable names for GMC
#' @example inst/examples/ex_lm.R
#' @export
rtmb_glm <- function(formula, data, family = "gaussian",
                       prior = prior_uniform(),
                       y_range = NULL,
                       init = NULL, fixed = NULL, null = NULL,
                       gmc = NULL,
                       factors = NULL,
                       contrasts = "treatment") {
  rtmb_glmer(formula = formula, data = data, family = family,
             laplace = FALSE,
             prior = prior,
             y_range = y_range,
             init = init,
             fixed = fixed,
             null = null,
             gmc = gmc,
             factors = factors,
             contrasts = contrasts)
}

#' RTMB-based Linear Regression wrapper function
#'
#' @param formula Formula (e.g., Y ~ X1 + X2)
#' @param data Data frame
#' @param prior An object of class "rtmb_prior" specifying the prior distribution. Use prior_weak(), prior_rhs(), or prior_ssp(). Default is NULL (flat prior).
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Specifying this automatically enables weakly informative priors.
#' @param init List of initial values.
#' @param null Character string specifying the target parameter for the null model.
#' @param gmc Character vector of variable names for GMC

#' @param factors Character vector of variable names to be treated as factors.
#' @example inst/examples/ex_lm.R
#' @export
rtmb_lm <- function(formula, data,
                    prior = prior_uniform(),
                    y_range = NULL,
                    init = NULL, fixed = NULL, null = NULL,
                    gmc = NULL,
                    factors = NULL,
                    contrasts = "treatment") {
  rtmb_glmer(formula = formula, data = data, family = "gaussian",
             laplace = FALSE,
             prior = prior,
             y_range = y_range,
             init = init,
             fixed = fixed,
             null = null,
             gmc = gmc,
             factors = factors,
             contrasts = contrasts)
}



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
#' @param prior Prior configuration: `prior_uniform()` (default) or `prior_weak()`.
#'   Hyperparameters can be specified within these functions (e.g., `prior_uniform(mean_sd = 10, sd_rate = 10)`).
#'   Available parameters for FA: `mean_sd`, `sd_rate`, `loadings_sd`, and `ssp_ratio` (if `rotate = "ssp"`).
#' @param y_range A numeric vector of length 2 specifying the theoretical min and max values of the items.
#'   Used to construct weakly informative priors when `prior = prior_weak()`.
#' @param init List of initial values.
#' @param fixed A named list of parameter values to fix (optional).
#' @example inst/examples/ex_fa.R
#' @export
rtmb_fa <- function(data, nfactors = 1, rotate = NULL, score = FALSE,
                    prior = prior_uniform(),
                    y_range = NULL,
                    init = NULL, fixed = NULL) {

  Y <- as.matrix(data)
  K <- nfactors
  J <- ncol(Y)
  N <- nrow(Y)

  if (K >= J) stop("The number of factors (K) must be less than the number of observed variables (J).")

  var_names <- colnames(data)
  if (is.null(var_names)) var_names <- paste0("V", 1:J)

  # --- 1. Pre-calculate Sufficient Statistics in R (Handling NAs) ---
  if (any(is.na(Y))) {
    Y_bar <- apply(Y, 2, mean, na.rm = TRUE)
    S_Y <- cov(Y, use = "pairwise.complete.obs") * (N - 1)
  } else {
    Y_bar <- apply(Y, 2, mean)
    S_Y <- cov(Y) * (N - 1)
  }

  # --- 2. Prior Handling (Extracting from prior object) ---
  if (is.null(prior)) prior <- prior_uniform()
  prior_type <- if (inherits(prior, "rtmb_prior")) prior$type else "weak"
  
  # New logic: If y_range is specified, automatically treat as "weak"
  if (!is.null(y_range)) {
    prior_type <- "weak"
  }
  
  # Error handling: if weak is requested but y_range is missing
  if (prior_type == "weak" && is.null(y_range)) {
    stop("When using 'prior_weak()', you must specify 'y_range' (e.g., y_range = c(1, 5)) to define the scaling of the priors.", call. = FALSE)
  }

  # Default/Initial settings
  prior_mean_center <- 0
  prior_mean_sd <- prior$mean_sd
  prior_sd_rate <- prior$sd_rate
  prior_loadings_sd <- prior$loadings_sd
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
    K_factors = K,
    prior_mean_center = prior_mean_center,
    prior_mean_sd = prior_mean_sd,
    prior_sd_rate = prior_sd_rate,
    prior_loadings_sd = prior_loadings_sd
  )

  # Determine if SSP model is used
  is_ssp <- !is.null(rotate) && rotate == "ssp"

  # --- 3. Simplified Setup AST ---
  setup_ast <- quote({
    N <- nrow(Y)
    J <- ncol(Y)
    # Note: Using colMeans and cov for efficiency in RTMB if possible, 
    # but follow user's requested syntax for setup block.
    Y_bar <- apply(Y, 2, mean)
    S_Y <- cov(Y) * (N - 1)
  })

  if (is_ssp) {
    dat_fa$ssp_ratio <- ssp_ratio
    if (score) dat_fa$Y <- Y

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
      Lambda_star = Dim(c(J, K_factors))
      r = Dim(c(J, K_factors), lower = 0.001, upper = 0.999)
      tau = Dim(c(J, K_factors), lower = 0)
      sd = Dim(J, lower = 0)
      CF_Omega = Dim(c(K_factors, K_factors), type = "CF_corr")
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
    model_exprs[[1]] <- quote(S_Y ~ sufficient_multi_normal_fa(N, Y_bar, mean, sd, Lambda %*% CF_Omega))
    if (!is.null(prior_mean_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(prior_mean_center, prior_mean_sd))
    if (!is.null(prior_sd_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(prior_sd_rate))
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
    score_expr <- if (score) {
      quote({
        # Use explicit matrix creation to ensure AD type is preserved
        Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE)
        out$score <- Y_c %*% solve(Sigma, Lambda %*% fa_cor)
      })
    } else quote({})
    ret_expr <- quote({ return(out) })
    gq_ast <- as.call(c(list(as.name("{")), as.list(base_gq)[-1], as.list(score_expr)[-1], as.list(ret_expr)[-1]))
    
    code_obj <- list(setup = setup_ast, parameters = param_ast, transform = tran_ast, model = model_ast, generate = gq_ast, env = parent.frame())
    p_names <- list(mean = var_names, Lambda_star = var_names, r = var_names, tau = var_names, sd = var_names, CF_Omega = paste0("Factor", 1:K), Lambda = var_names, L = var_names, fa_cor = paste0("Factor", 1:K), communality = var_names)
    if (score) {
      ind_names <- rownames(data); if (is.null(ind_names)) ind_names <- paste0("Id", 1:N)
      p_names[["score"]] <- list(ind_names, paste0("Factor", 1:K))
    }
    obj <- rtmb_model(data = dat_fa, code = code_obj, par_names = p_names, init = init, view = c("L", "sd", "fa_cor"), fixed = fixed)
    return(obj)

  } else {
    # --- Standard rotation logic ---
    if (score) dat_fa$Y <- Y
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
      L_raw <- Dim(dim = c(J, K_factors), type = "lower_tri")
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
    model_exprs[[1]] <- quote(S_Y ~ sufficient_multi_normal_fa(N, Y_bar, mean, sd, L_raw))
    if (!is.null(prior_mean_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(0, prior_mean_sd))
    if (!is.null(prior_sd_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(prior_sd_rate))
    if (!is.null(prior_loadings_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(L_raw ~ lower_tri_normal(0, prior_loadings_sd))
    model_ast <- as.call(c(list(as.name("{")), model_exprs))

    base_gq <- quote({
      Sigma <- L_raw %*% t(L_raw) + diag(sd^2)
      var_total <- diag(Sigma); var_common <- rowSums(L_raw^2); communality <- var_common / var_total
      out <- list(communality = communality)
    })

    if (!is.null(rotate)) {
      rot_loadings_name <- paste0("L_", rotate)
      if (exists(rotate, mode = "function")) { rot_fn <- match.fun(rotate); fn_call <- as.name(rotate) }
      else if (requireNamespace("GPArotation", quietly = TRUE) && exists(rotate, where = asNamespace("GPArotation"), mode = "function")) {
        rot_fn <- getFromNamespace(rotate, "GPArotation"); fn_call <- call("::", as.name("GPArotation"), as.name(rotate))
      } else stop("Rotation function not found: ", rotate)
      
      dummy_L <- matrix(rnorm(J * K), J, K); test_rot <- rot_fn(dummy_L)
      is_matrix_rot <- is.matrix(test_rot); has_phi <- !is_matrix_rot && !is.null(test_rot$Phi)
      if (is_matrix_rot) {
        rot_expr <- bquote({ rot_obj <- .(fn_call)(L); out[[.(rot_loadings_name)]] <- unclass(rot_obj) })
        score_expr <- if (score) bquote({
          Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE)
          rot_raw <- unclass(.(fn_call)(L_raw)); if (!is.matrix(rot_raw)) rot_raw <- unclass(rot_raw$loadings)
          out$score <- Y_c %*% solve(Sigma, rot_raw)
        }) else quote({})
      } else {
        if (has_phi) {
          rot_expr <- bquote({ rot_obj <- .(fn_call)(L); out$fa_cor <- rot_obj$Phi; out[[.(rot_loadings_name)]] <- unclass(rot_obj$loadings) })
          score_expr <- if (score) bquote({
            Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE)
            rot_raw_obj <- .(fn_call)(L_raw); out$score <- Y_c %*% solve(Sigma, unclass(rot_raw_obj$loadings) %*% rot_raw_obj$Phi)
          }) else quote({})
        } else {
          rot_expr <- bquote({ rot_obj <- .(fn_call)(L); out[[.(rot_loadings_name)]] <- unclass(rot_obj$loadings) })
          score_expr <- if (score) bquote({
            Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE)
            out$score <- Y_c %*% solve(Sigma, unclass(.(fn_call)(L_raw)$loadings))
          }) else quote({})
        }
      }
    } else { has_phi <- FALSE; rot_expr <- quote({}); score_expr <- if (score) quote({ Y_c <- Y - matrix(mean, nrow = N, ncol = J, byrow = TRUE); out$score <- Y_c %*% solve(Sigma, L_raw) }) else quote({}) }

    ret_expr <- quote({ return(out) })
    gq_ast <- as.call(c(list(as.name("{")), as.list(base_gq)[-1], as.list(rot_expr)[-1], as.list(score_expr)[-1], as.list(ret_expr)[-1]))
    
    code_obj <- list(setup = setup_ast, parameters = param_ast, transform = tran_ast, model = model_ast, generate = gq_ast, env = parent.frame())
    p_names <- list(mean = var_names, L_raw = var_names, sd = var_names, L = var_names, communality = var_names)
    if (!is.null(rotate)) { p_names[[paste0("L_", rotate)]] <- var_names; if (has_phi) p_names[["fa_cor"]] <- paste0("Factor", 1:K) }
    if (score) { ind_names <- rownames(data); if (is.null(ind_names)) ind_names <- paste0("Id", 1:N); p_names[["score"]] <- list(ind_names, paste0("Factor", 1:K)) }

    target_view <- if (!is.null(rotate)) c(paste0("L_", rotate), "sd", "fa_cor") else c("L", "sd", "fa_cor")
    obj <- rtmb_model(data = dat_fa, code = code_obj, par_names = p_names, init = init, view = target_view, fixed = fixed)
    return(obj)
  }
}

#' RTMB-based IRT (Item Response Theory) Wrapper
#'
#' @description
#' Performs Item Response Theory modeling (1PL, 2PL, 3PL, and Graded Response Model) using RTMB.
#' Missing values (NA) in the data are automatically removed internally for efficient computation.
#'
#' @param data A data frame or matrix of item responses (N persons x J items).
#' @param model Character string for the model type: "1PL", "2PL", or "3PL".
#' @param type Character string for the data type: "binary" or "ordered".
#' @param prior Prior configuration: `prior_uniform()` (default) or `prior_weak()`.
#'   Hyperparameters can be specified within these functions (e.g., `prior_weak(b_sd = 5)`).
#'   Available parameters for IRT: `a_rate` (discrimination), `b_sd` (difficulty), `c_alpha`/`c_beta` (guessing).
#' @param init List of initial values.
#' @param fixed A named list of parameter values to fix (optional).
#' @param view Character vector of parameter names to prioritize in summary.
#' @param ... Additional arguments passed to \code{rtmb_model()}.
#' @example inst/examples/ex_irt.R
#' @export
rtmb_irt <- function(data, model = c("2PL", "1PL", "3PL"), type = c("binary", "ordered"),
                     prior = prior_uniform(), 
                     init = NULL, fixed = NULL, view = NULL, ...) {

  model <- match.arg(model)
  type <- match.arg(type)

  if (model == "3PL" && type == "ordered") {
    stop("The 3PL model only supports 'binary' data.")
  }

  Y <- as.matrix(data)

  # Get row and column names
  item_names <- colnames(Y)
  if (is.null(item_names)) item_names <- paste0("Item", 1:ncol(Y))
  person_names <- rownames(Y)
  if (is.null(person_names)) person_names <- paste0("Person", 1:nrow(Y))

  # Prior Handling
  if (is.null(prior)) prior <- prior_uniform()
  prior_type <- if (inherits(prior, "rtmb_prior")) prior$type else "weak"
  
  # Extract hyperparameters from prior object if they exist
  # Use priority: 1. Explicitly provided in prior_weak(...), 2. Default weak values (if type is weak)
  a_rate  <- prior$a_rate
  b_sd    <- prior$b_sd
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
  setup_exprs[[5]] <- quote(N_persons <- nrow(Y)) # or max(person_idx)
  setup_exprs[[6]] <- quote(N_items <- ncol(Y))   # or max(item_idx)
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
  model_exprs[[length(model_exprs) + 1]] <- bquote(theta ~ normal(0, .(theta_sd)))

  model_ast <- as.call(c(list(as.name("{")), model_exprs))

  # --- Mapping of Parameter Names ---
  code_obj <- list(setup = setup_ast, parameters = param_ast, model = model_ast)

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

  obj <- rtmb_model(data = as.list(tmp_env), code = code_obj, par_names = par_names_list, 
                    init = init, fixed = fixed, view = view_vars)

  return(obj)
}

#' Fit a Correlation Model using RTMB
#'
#' @description
#' `rtmb_corr` fits a correlation model to estimate means, standard deviations, and correlation structures.
#' It supports simple correlation, multilevel correlation, and classical frequentist estimation.
#'
#' @param x A matrix, data frame, formula, or expression (e.g., \code{cbind(V1, V2)}) of response variables.
#' @param data An optional data frame containing the variables.
#' @param ID A character string or expression specifying the group ID variable for multilevel models.
#' @param prior Prior configuration object: \code{prior_uniform()} or \code{prior_weak()}. Default is \code{prior_uniform()}.
#' @param y_range Optional numeric vector or matrix defining the theoretical range (min, max) of response variables.
#' Required when using \code{prior_weak()}. Can be a vector of length 2 (applies to all variables) or a matrix/list of length P.
#' @param init Optional list of initial values.
#' @param null Optional list specifying parameters to fix to null values.

#' @param ... Additional arguments passed to \code{rtmb_model}.
#' @return A \code{RTMB_Model} object.
#' @export
rtmb_corr <- function(x = NULL, data = NULL, ID = NULL,
                      covariates = NULL,
                      prior = prior_uniform(), y_range = NULL,
                      init = NULL, fixed = NULL, null = NULL, ...) {

  x_expr <- substitute(x)
  id_expr <- substitute(ID)

  # Evaluation logic for response variables (Y_mat)
  if (!is.null(data)) {
    # Evaluate x in the context of data (NSE support for cbind(a, b))
    Y_mat <- eval(x_expr, data, parent.frame())
  } else {
    Y_mat <- x
  }

  # Handle if Y_mat is a formula
  if (inherits(Y_mat, "formula")) {
    formula <- Y_mat
    lhs_expr <- formula[[2]]

    # 1. Try to evaluate LHS directly (handles cbind(y1, y2) or Y_df)
    Y_mat <- try(eval(lhs_expr, data, parent.frame()), silent = TRUE)

    # 2. If eval(lhs) failed, try a safer model.frame call
    if (inherits(Y_mat, "try-error") || is.null(Y_mat)) {
      resp_formula <- bquote(.(lhs_expr) ~ 1)
      mf_y <- try(model.frame(resp_formula, data = data, na.action = na.pass), silent = TRUE)
      if (!inherits(mf_y, "try-error")) {
         Y_mat <- model.response(mf_y)
         if (is.null(Y_mat)) Y_mat <- mf_y[, 1, drop = FALSE]
      } else {
         # Last resort: check if it's a name in data
         id_name <- as.character(lhs_expr)
         if (!is.null(data) && id_name %in% names(data)) {
            Y_mat <- data[[id_name]]
         }
      }
    }

    # Final check if it's a list (not df)
    if (is.list(Y_mat) && !is.data.frame(Y_mat)) {
       Y_mat <- do.call(cbind, Y_mat)
    }

    # If formula has RHS, extract covariates
    if (length(formula) == 3 && is.null(covariates)) {
       rhs_expr <- formula[[3]]
       cov_formula <- bquote(~ .(rhs_expr))
       mf_x <- try(model.frame(cov_formula, data = data, na.action = na.pass), silent = TRUE)
       if (inherits(mf_x, "try-error")) {
          # Try to subset data to only RHS variables to avoid errors from other complex columns
          rhs_vars <- all.vars(rhs_expr)
          data_sub <- if (!is.null(data)) data[, intersect(rhs_vars, names(data)), drop = FALSE] else data
          mf_x <- model.frame(cov_formula, data = data_sub, na.action = na.pass)
       }
       covariates <- model.matrix(attr(mf_x, "terms"), mf_x)
       # Remove intercept
       if ("(Intercept)" %in% colnames(covariates)) {
         covariates <- covariates[, colnames(covariates) != "(Intercept)", drop = FALSE]
       }
    }
  } else {
    formula <- NULL
    # If Y_mat is a character vector of names, subset from data
    if (is.character(Y_mat) && length(Y_mat) > 1 && !is.null(data)) {
      Y_mat <- data[, Y_mat, drop = FALSE]
    }
    # Handle list of vectors
    if (is.list(Y_mat) && !is.data.frame(Y_mat)) {
       Y_mat <- do.call(cbind, Y_mat)
    }
  }

  # Handle covariates (if provided separately or extracted from formula)
  X_mat <- NULL
  if (!is.null(covariates)) {
     if (is.matrix(covariates) || is.data.frame(covariates)) {
        X_mat <- as.matrix(covariates)
     } else if (inherits(covariates, "formula")) {
        mf_x <- model.frame(covariates, data = data, na.action = na.pass)
        X_mat <- model.matrix(covariates, mf_x)
        if ("(Intercept)" %in% colnames(X_mat)) {
          X_mat <- X_mat[, colnames(X_mat) != "(Intercept)", drop = FALSE]
        }
     } else if (is.character(covariates) && !is.null(data)) {
        X_mat <- as.matrix(data[, covariates, drop = FALSE])
     }
  }

  if (!is.null(X_mat)) {
    # Ensure numeric
    X_mat <- X_mat[, sapply(as.data.frame(X_mat), is.numeric), drop = FALSE]
  }

  # Parse ID (NSE support for ID = group)
  id_val <- NULL
  if (!is.null(ID)) {
    if (!is.null(data)) {
      # Try to evaluate ID in data
      id_val <- try(eval(id_expr, data, parent.frame()), silent = TRUE)
      if (inherits(id_val, "try-error") || is.null(id_val)) {
        # Fallback: check if ID name is in data
        id_name <- as.character(id_expr)
        if (id_name %in% names(data)) {
          id_val <- data[[id_name]]
        } else {
          id_val <- ID # Literal value
        }
      }
    } else if (is.data.frame(Y_mat)) {
      # If data is NULL but x is a dataframe, check if ID is a column name
      id_name <- as.character(id_expr)
      if (id_name %in% names(Y_mat)) {
        id_val <- Y_mat[[id_name]]
        # We will remove this column from Y_mat later
      } else {
        id_val <- ID
      }
    } else {
      id_val <- ID
    }

    # If the response matrix still contains the ID column, remove it
    id_name_str <- as.character(id_expr)
    if (is.data.frame(Y_mat) && id_name_str %in% colnames(Y_mat)) {
       Y_mat[[id_name_str]] <- NULL
    }
  }

  # Ensure Y_mat is a numeric matrix
  Y_mat <- as.data.frame(Y_mat)
  num_cols <- sapply(Y_mat, is.numeric)
  Y_mat <- Y_mat[, num_cols, drop = FALSE]
  Y_mat <- as.matrix(Y_mat)

  P_y <- ncol(Y_mat)
  target_names <- colnames(Y_mat)
  if (is.null(target_names)) target_names <- paste0("Y", 1:P_y)

  # Combine with covariates for Joint MVN
  if (!is.null(X_mat)) {
     P_x <- ncol(X_mat)
     control_names <- colnames(X_mat)
     if (is.null(control_names)) control_names <- paste0("X", 1:P_x)
     Y_mat <- cbind(Y_mat, X_mat)
     colnames(Y_mat) <- c(target_names, control_names)
  } else {
     P_x <- 0
     control_names <- NULL
  }

  N <- nrow(Y_mat)
  P <- ncol(Y_mat)
  if (P < 1) stop("No numeric columns found for correlation analysis.")

  var_names <- colnames(Y_mat)
  if (is.null(var_names)) var_names <- paste0("V", 1:P)
  colnames(Y_mat) <- var_names



  # --- RTMB Implementation ---
  if (is.null(prior)) {
    prior <- prior_uniform()
  }

  if (!is.null(id_val)) {
     # Multilevel correlation mode
     group_factor <- as.factor(id_val)
     group_id <- as.integer(group_factor)
     group_names <- levels(group_factor)
     J <- length(group_names)

     # Automatically switch to prior_weak() if y_range is provided and prior is default uniform
     if (!is.null(y_range) && inherits(prior, "rtmb_prior") && prior$type == "uniform" &&
         is.null(prior$Intercept_sd) && is.null(prior$b_sd) &&
         is.null(prior$sigma_rate) && is.null(prior$tau_rate)) {
       prior <- prior_weak()
     }

     if (!inherits(prior, "rtmb_prior")) {
       stop("prior must be an object of class 'rtmb_prior'. Use prior_weak() or prior_uniform().")
     }

     default_prior_settings <- list(Intercept_sd = 10, mu_sd = 10, b_sd = 10, sigma_rate = 5, sd_rate = 5, lkj_eta = 1.0, sd_ratio = 0.5)
     prior <- .merge_prior(default_prior_settings, prior)

     if (!is.null(prior$Intercept_sd) && prior$Intercept_sd != 10) prior$mu_sd <- prior$Intercept_sd
     if (!is.null(prior$mu_sd) && prior$mu_sd != 10) prior$Intercept_sd <- prior$mu_sd

     prior_type <- prior$type
     use_weak_info <- prior_type %in% c("weak")
     multivariate <- P > 1

     setup_exprs <- list()
     if (use_weak_info) {
       if (is.null(y_range)) {
         stop("Specifying 'y_range' is required when using prior_weak().")
       }
       if (is.list(y_range)) {
         Y_range_mat <- do.call(rbind, lapply(y_range, as.numeric))
       } else if (is.vector(y_range) && length(y_range) == 2) {
         Y_range_mat <- matrix(rep(as.numeric(y_range), each = P), P, 2)
       } else {
         Y_range_mat <- as.matrix(y_range)
       }

       half_d_y <- (Y_range_mat[, 2] - Y_range_mat[, 1]) / 2
       mid_y_val <- (Y_range_mat[, 2] + Y_range_mat[, 1]) / 2
       base_scale <- half_d_y * prior$sd_ratio

       setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mid_y)
       setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- alpha_prior_sd)
       setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate_vec <- sigma_rate_vec)
     }
     setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

     param_exprs <- list(as.name("{"))
     param_exprs[[length(param_exprs) + 1]] <- bquote(mu <- Dim(.(P), random = TRUE))
     param_exprs[[length(param_exprs) + 1]] <- bquote(u <- Dim(c(J, .(P)), random = TRUE))

     param_exprs[[length(param_exprs) + 1]] <- bquote(sigma_between <- Dim(.(P), lower = 0))
     param_exprs[[length(param_exprs) + 1]] <- bquote(sigma_within <- Dim(.(P), lower = 0))

     if (multivariate) {
       param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr_between <- Dim(c(.(P), .(P)), type = "CF_corr"))
       param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr_within <- Dim(c(.(P), .(P)), type = "CF_corr"))
     }
     param_ast <- as.call(param_exprs)

     model_exprs <- list(as.name("{"))
     if (multivariate) {
       model_exprs[[length(model_exprs) + 1]] <- bquote(L_corr_between ~ lkj_CF_corr(.(prior$lkj_eta)))
       model_exprs[[length(model_exprs) + 1]] <- bquote(L_corr_within ~ lkj_CF_corr(.(prior$lkj_eta)))
       model_exprs[[length(model_exprs) + 1]] <- quote(u ~ multi_normal_CF(mean = rep(0, P), sd = sigma_between, CF_Omega = L_corr_between))
     } else {
       model_exprs[[length(model_exprs) + 1]] <- quote(u ~ normal(0, sigma_between))
     }

     model_exprs[[length(model_exprs) + 1]] <- quote(Y_pred <- u[group_id, ] + rep(mu, each = N))
     if (multivariate) {
       model_exprs[[length(model_exprs) + 1]] <- quote(Y ~ multi_normal_CF(mean = Y_pred, sd = sigma_within, CF_Omega = L_corr_within))
     } else {
       model_exprs[[length(model_exprs) + 1]] <- quote(Y ~ normal(Y_pred, sigma_within))
     }

     if (prior_type == "weak") {
       model_exprs[[length(model_exprs) + 1]] <- quote(sigma_between ~ exponential(sigma_rate_vec))
       model_exprs[[length(model_exprs) + 1]] <- quote(sigma_within ~ exponential(sigma_rate_vec))
       model_exprs[[length(model_exprs) + 1]] <- quote(mu ~ normal(mid_y, alpha_prior_sd))
     } else if (prior_type == "uniform") {
       if (!is.null(prior$sigma_rate)) {
         model_exprs[[length(model_exprs) + 1]] <- bquote(sigma_between ~ exponential(.(prior$sigma_rate)))
         model_exprs[[length(model_exprs) + 1]] <- bquote(sigma_within ~ exponential(.(prior$sigma_rate)))
       }
       if (!is.null(prior$Intercept_sd)) {
         model_exprs[[length(model_exprs) + 1]] <- bquote(mu ~ normal(0, .(prior$Intercept_sd)))
       }
     }
     model_ast <- as.call(model_exprs)

      generate_exprs <- list(as.name("{"))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(ICC <- (sigma_between^2) / (sigma_between^2 + sigma_within^2))
      if (multivariate) {
        generate_exprs[[length(generate_exprs) + 1]] <- quote(B_corr <- L_corr_between %*% t(L_corr_between))
        generate_exprs[[length(generate_exprs) + 1]] <- quote(W_corr <- L_corr_within %*% t(L_corr_within))

        if (P_x > 0) {
           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_yy <- W_corr[1:P_y, 1:P_y])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_yx <- W_corr[1:P_y, (P_y+1):P])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_xx <- W_corr[(P_y+1):P, (P_y+1):P])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(P_cov_w <- R_yy - R_yx %*% solve(R_xx) %*% t(R_yx))
           generate_exprs[[length(generate_exprs) + 1]] <- quote(D_w <- diag(1 / sqrt(diag(P_cov_w))))
           generate_exprs[[length(generate_exprs) + 1]] <- quote(W_pcorr <- D_w %*% P_cov_w %*% D_w)

           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_yy_b <- B_corr[1:P_y, 1:P_y])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_yx_b <- B_corr[1:P_y, (P_y+1):P])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_xx_b <- B_corr[(P_y+1):P, (P_y+1):P])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(P_cov_b <- R_yy_b - R_yx_b %*% solve(R_xx_b) %*% t(R_yx_b))
           generate_exprs[[length(generate_exprs) + 1]] <- quote(D_b <- diag(1 / sqrt(diag(P_cov_b))))
           generate_exprs[[length(generate_exprs) + 1]] <- quote(B_pcorr <- D_b %*% P_cov_b %*% D_b)
           generate_exprs[[length(generate_exprs) + 1]] <- quote(list(ICC = ICC, B_corr = B_corr, W_corr = W_corr, W_pcorr = W_pcorr, B_pcorr = B_pcorr))
        } else {
           generate_exprs[[length(generate_exprs) + 1]] <- quote(list(ICC = ICC, B_corr = B_corr, W_corr = W_corr))
        }
      } else {
        generate_exprs[[length(generate_exprs) + 1]] <- quote(list(ICC = ICC))
      }
      generate_ast <- as.call(generate_exprs)

     mdl_code <- list(setup = setup_ast, parameters = param_ast, model = model_ast, generate = generate_ast, env = parent.frame())
     class(mdl_code) <- "rtmb_code"

     v_names <- list(mu = var_names, sigma_between = var_names, sigma_within = var_names)
     v_names$u <- list(group_names, var_names)
     v_names$ICC <- var_names

     if (multivariate) {
       v_names$B_corr <- list(var_names, var_names)
       v_names$W_corr <- list(var_names, var_names)
       if (P_x > 0) {
          v_names$B_pcorr <- list(target_names, target_names)
          v_names$W_pcorr <- list(target_names, target_names)
       }
     }

     data_list <- list(Y = Y_mat, group_id = group_id, N = N, P = P, J = J, P_y = P_y, P_x = P_x)
     if (use_weak_info) {
       data_list$mid_y <- mid_y_val
       data_list$alpha_prior_sd <- half_d_y
       data_list$sigma_rate_vec <- 1.0 / base_scale
     }

     init_list <- list(mu = colMeans(Y_mat))
     init_list$sigma_between <- apply(Y_mat, 2, sd) * 0.5
     init_list$sigma_within <- apply(Y_mat, 2, sd) * 0.5
     init_list$u <- matrix(0, J, P)

     view_order <- c("ICC")
     if (multivariate) {
       if (P_x > 0) view_order <- c(view_order, "W_pcorr", "B_pcorr")
       view_order <- c(view_order, "B_corr", "W_corr")
     }
     view_order <- c(view_order, "mu", "sigma_between", "sigma_within")

     obj <- rtmb_model(data_list, mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order)
     obj$raw_data <- data

     obj$type <- "corr"
     obj$extra$df_pars <- "mu"

     return(obj)
  } else {
     # Simple correlation mode
     if (!is.null(y_range) && inherits(prior, "rtmb_prior") && prior$type == "uniform" &&
         is.null(prior$Intercept_sd) && is.null(prior$mu_sd) && is.null(prior$b_sd) &&
         is.null(prior$sigma_rate)) {
       prior <- prior_weak()
     }

     if (!inherits(prior, "rtmb_prior")) {
       stop("prior must be an object of class 'rtmb_prior'. Use prior_weak() or prior_uniform().")
     }

     prior_type <- prior$type
     default_prior <- list(Intercept_sd = 10, mu_sd = 10, sigma_rate = 1.0, lkj_eta = 1.0)
     prior <- .merge_prior(default_prior, prior)
     if (!is.null(prior$Intercept_sd) && prior$Intercept_sd != 10) prior$mu_sd <- prior$Intercept_sd
     if (!is.null(prior$mu_sd) && prior$mu_sd != 10) prior$Intercept_sd <- prior$mu_sd

     setup_exprs <- list(as.name("{"))
     setup_exprs[[length(setup_exprs) + 1]] <- quote(N <- N)
     setup_exprs[[length(setup_exprs) + 1]] <- quote(P <- P)

     if (prior_type == "weak") {
       if (is.null(y_range)) {
         stop("y_range is required when using prior_weak().")
       }
       if (is.list(y_range)) {
         y_range_mat <- do.call(rbind, y_range)
       } else if (is.vector(y_range) && length(y_range) == 2) {
         y_range_mat <- matrix(y_range, P, 2, byrow = TRUE)
       } else {
         y_range_mat <- as.matrix(y_range)
       }

       mid_y_val <- (y_range_mat[, 1] + y_range_mat[, 2]) / 2
       alpha_prior_sd_val <- (y_range_mat[, 2] - y_range_mat[, 1]) / 2
       sigma_rate_val <- 1 / (alpha_prior_sd_val * 0.5)

       setup_exprs[[length(setup_exprs) + 1]] <- quote(mid_y <- mid_y)
       setup_exprs[[length(setup_exprs) + 1]] <- quote(alpha_prior_sd <- alpha_prior_sd)
       setup_exprs[[length(setup_exprs) + 1]] <- quote(sigma_rate_vec <- sigma_rate_vec)
     }
     setup_ast <- as.call(setup_exprs)

     param_exprs <- list(as.name("{"))
     param_exprs[[length(param_exprs) + 1]] <- bquote(mean <- Dim(.(P), random = TRUE))
     param_exprs[[length(param_exprs) + 1]] <- bquote(sd   <- Dim(.(P), lower = 0))

     if (P == 2) {
       param_exprs[[length(param_exprs) + 1]] <- bquote(corr <- Dim(lower = -1, upper = 1))
     } else {
       param_exprs[[length(param_exprs) + 1]] <- bquote(CF_corr <- Dim(c(.(P), .(P)), type = "CF_corr"))
     }
     param_ast <- as.call(param_exprs)

     model_exprs <- list(as.name("{"))
     if (prior_type == "weak") {
       model_exprs[[length(model_exprs) + 1]] <- quote(mean ~ normal(mid_y, alpha_prior_sd))
       model_exprs[[length(model_exprs) + 1]] <- quote(sd ~ exponential(sigma_rate_vec))
     } else {
       if (!is.null(prior$mu_sd)) {
         model_exprs[[length(model_exprs) + 1]] <- bquote(mean ~ normal(0, .(prior$mu_sd)))
       }
       if (!is.null(prior$sigma_rate)) {
         model_exprs[[length(model_exprs) + 1]] <- bquote(sd ~ exponential(.(prior$sigma_rate)))
       }
     }

     if (P == 2) {
       model_exprs[[length(model_exprs) + 1]] <- bquote(corr ~ lkj_corr(.(prior$lkj_eta)))
       model_exprs[[length(model_exprs) + 1]] <- quote(CF_corr <- matrix(corr * 0, nrow = 2, ncol = 2))
       model_exprs[[length(model_exprs) + 1]] <- quote(CF_corr[1, 1] <- 1)
       model_exprs[[length(model_exprs) + 1]] <- quote(CF_corr[2, 1] <- corr)
       model_exprs[[length(model_exprs) + 1]] <- quote(CF_corr[2, 2] <- sqrt(1 - corr^2))
     } else {
       model_exprs[[length(model_exprs) + 1]] <- bquote(CF_corr ~ lkj_CF_corr(.(prior$lkj_eta)))
     }

     model_exprs[[length(model_exprs) + 1]] <- quote(S_Y ~ sufficient_multi_normal_CF(N, Y_bar, mean, sd, CF_corr))
     model_ast <- as.call(model_exprs)

     generate_exprs <- list(as.name("{"))
      if (P > 2) {
        generate_exprs[[length(generate_exprs) + 1]] <- quote(M_CF <- as.matrix(CF_corr))
        generate_exprs[[length(generate_exprs) + 1]] <- quote(corr <- M_CF %*% t(M_CF))

        if (P_x > 0) {
           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_yy <- corr[1:P_y, 1:P_y])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_yx <- corr[1:P_y, (P_y+1):P])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(R_xx <- corr[(P_y+1):P, (P_y+1):P])
           generate_exprs[[length(generate_exprs) + 1]] <- quote(P_cov <- R_yy - R_yx %*% solve(R_xx) %*% t(R_yx))
           generate_exprs[[length(generate_exprs) + 1]] <- quote(D <- diag(1 / sqrt(diag(P_cov))))
           generate_exprs[[length(generate_exprs) + 1]] <- quote(pcorr <- D %*% P_cov %*% D)
           generate_exprs[[length(generate_exprs) + 1]] <- quote(list(corr = corr, pcorr = pcorr))
        } else {
           generate_exprs[[length(generate_exprs) + 1]] <- quote(list(corr = corr))
        }
      }
      generate_ast <- as.call(generate_exprs)

     mdl_code <- list(setup = setup_ast, parameters = param_ast, model = model_ast, generate = generate_ast, env = parent.frame())
     class(mdl_code) <- "rtmb_code"

     dat_list <- list(
       N = N, P = P, Y_bar = colMeans(Y_mat), S_Y = cov(Y_mat) * (N - 1), P_y = P_y, P_x = P_x
     )
     if (prior_type == "weak") {
       dat_list$mid_y <- mid_y_val
       dat_list$alpha_prior_sd <- alpha_prior_sd_val
       dat_list$sigma_rate_vec <- sigma_rate_val
     }

     v_names <- list(mean = var_names, sd = var_names)
     if (P == 2) {
       v_names$corr <- "rho"
     } else {
       v_names$corr <- list(var_names, var_names)
       if (P_x > 0) v_names$pcorr <- list(target_names, target_names)
     }

     init_list <- if (is.null(init)) {
       list(mean = colMeans(Y_mat), sd = apply(Y_mat, 2, sd))
     } else init

     view_vars <- if (P_x > 0) "pcorr" else "corr"
     obj <- rtmb_model(data = dat_list, code = mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_vars)

     if (!is.null(null)) {
       obj <- obj$null_model(target = null)
     }

    obj$raw_data <- data
    obj$type <- "corr"
    obj$extra$df_pars <- "mean"
    obj$extra$lkj_eta <- prior$lkj_eta

    return(obj)
  }
}

#' RTMB-based Bayesian two-sample t-test wrapper function
#'
#' @description
#' Performs a Bayesian or Frequentist two-sample t-test using RTMB.
#'
#' @param x Numeric vector of responses for group 1, a formula (e.g., `y ~ group`), or a column name (unquoted) if `data` is provided.
#' @param y Numeric vector of responses for group 2, or a column name (unquoted) if `data` is provided. Required if `x` is not a formula.
#' @param data Data frame containing the variables.
#' @param r Numeric; Cauchy prior scale for the effect size (delta). Default is 0.707.
#' @param paired Logical; whether to perform a paired t-test.
#' @param ID Character; name of the ID variable for paired t-tests (required for formula input with paired = TRUE).
#' @param y_range Theoretical minimum and maximum values of the response variable as a vector c(min, max). Required when using weakly informative priors.
#' @param prior An object of class "rtmb_prior" (e.g., `prior_uniform()`, `prior_jzs()`, or `prior_weak()`).
#' @param init List of initial values.
#' @param null Character string specifying the target parameter for the null model (e.g., "delta").
#' @param var.equal Logical; whether to assume equal variances. Default is TRUE.
#' @param ... Additional arguments.
#' @return An `RTMB_Model` object.
#' @export
rtmb_ttest <- function(x, y = NULL, data = NULL, r = 0.707,
                       paired = FALSE, ID = NULL,
                       y_range = NULL,
                       prior = prior_uniform(),
                       init = NULL, fixed = NULL, null = NULL,
                       var.equal = TRUE, ...) {

  x_expr <- substitute(x)
  y_expr <- substitute(y)

  # --- 1. Data Extraction ---
  is_formula <- is.call(x_expr) && (x_expr[[1]] == quote(`~`))
  if (is_formula) {
    x <- eval(x_expr, parent.frame())
    mf <- if (is.null(data)) model.frame(x, parent.frame()) else model.frame(x, data)
    response <- mf[[1]]; group <- as.factor(mf[[2]])
    levs <- levels(group); x_label <- levs[1]; y_label <- levs[2]
    if (paired) {
      id_val <- if (!is.null(data)) data[[ID]] else eval(substitute(ID), parent.frame())
      d1 <- mf[group == levs[1], ]; d1$id <- id_val[group == levs[1]]
      d2 <- mf[group == levs[2], ]; d2$id <- id_val[group == levs[2]]
      common_ids <- intersect(d1$id, d2$id)
      Y1 <- d1[match(common_ids, d1$id), 1]; Y2 <- d2[match(common_ids, d2$id), 1]
    } else {
      Y1 <- as.numeric(na.omit(response[group == levs[1]]))
      Y2 <- as.numeric(na.omit(response[group == levs[2]]))
    }
  } else {
    x_label <- deparse(x_expr); y_label <- deparse(y_expr)
    if (!is.null(data)) {
      Y1 <- if (x_label %in% names(data)) data[[x_label]] else eval(x_expr, parent.frame())
      Y2 <- if (y_label %in% names(data)) data[[y_label]] else eval(y_expr, parent.frame())
    } else {
      Y1 <- eval(x_expr, parent.frame()); Y2 <- eval(y_expr, parent.frame())
    }
    Y1 <- as.numeric(na.omit(Y1)); Y2 <- as.numeric(na.omit(Y2))
  }
  levs <- c(x_label, y_label)

  # --- 2. Prior Setup ---
  if (is.null(prior)) prior <- prior_uniform()

  # Automatically switch to prior_weak() if y_range is provided and prior is default uniform
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && prior$type == "uniform") {
    prior <- prior_weak()
  }
  prior_type <- prior$type
  is_jzs <- prior_type == "jzs"
  is_weak <- prior_type == "weak"
  use_delta_param <- is_jzs
  if (is_jzs) r <- prior$r

  # --- 3. AST Construction ---
  setup_ast <- quote({})
  if (is_weak) {
    if (is.null(y_range)) {
      stop("Specifying 'y_range' (theoretical minimum and maximum values of the response) is required when using prior_weak().", call. = FALSE)
    }
    
    # Common setup values
    sd_ratio_val <- prior$sd_ratio
    max_beta_val <- prior$max_beta
    
    setup_ast <- bquote({
      half_d_y <- .(diff(y_range) / 2)
      base_scale <- half_d_y * .(sd_ratio_val)
      alpha_prior_sd <- half_d_y
      
      # Determine x_sd based on design
      .(if (paired) {
        quote(x_sd <- 0.5) # Balanced Condition effect (Condition 0/1)
      } else {
        quote({
          n1 <- length(Y1)
          n2 <- length(Y2)
          p1 <- n1 / (n1 + n2)
          x_sd <- sqrt(p1 * (1 - p1))
        })
      })
      
      diff_prior_sd <- .(max_beta_val) * base_scale / x_sd
      mid_y <- .(if(paired) 0 else mean(y_range))
      sigma_rate_weak <- 1.0 / base_scale
    })
  }

  if (paired) {
    dat <- list(diffs = Y1 - Y2, r = r)
    param_ast <- if (use_delta_param) {
      quote({ sd_diff = Dim(1, lower = 0); delta = Dim(1) })
    } else {
      quote({ sd_diff = Dim(1, lower = 0); diff = Dim(1) })
    }
    tran_ast <- if (use_delta_param) {
      quote({ diff <- delta * sd_diff; report(diff) })
    } else {
      quote({ delta <- diff / sd_diff; report(delta) })
    }
    model_body <- list(quote(diffs ~ normal(diff, sd_diff)))
    if (is_jzs) model_body[[length(model_body)+1]] <- quote(delta ~ cauchy(0, r))
    if (is_weak) {
      model_body[[length(model_body)+1]] <- quote(diff ~ normal(0, diff_prior_sd))
      model_body[[length(model_body)+1]] <- quote(sd_diff ~ exponential(sigma_rate_weak))
    }
    view_vars <- c("diff", "delta", "sd_diff")
    df_pars <- if (use_delta_param) "delta" else "diff"
  } else {
    dat <- list(Y1 = Y1, Y2 = Y2, r = r)
    param_ast <- if (use_delta_param) {
      quote({ mean = Dim(1); sd = Dim(1, lower = 0); delta = Dim(1) })
    } else {
      quote({ mean = Dim(1); sd = Dim(1, lower = 0); diff = Dim(1) })
    }
    tran_ast <- if (use_delta_param) {
      quote({ diff <- delta * sd; mean0 <- mean + diff/2; mean1 <- mean - diff/2; report(diff) })
    } else {
      quote({ delta <- diff / sd; mean0 <- mean + diff/2; mean1 <- mean - diff/2; report(delta) })
    }
    model_body <- list(quote(Y1 ~ normal(mean0, sd)), quote(Y2 ~ normal(mean1, sd)))
    if (is_jzs) model_body[[length(model_body)+1]] <- quote(delta ~ cauchy(0, r))
    if (is_weak) {
      model_body[[length(model_body)+1]] <- quote(mean ~ normal(mid_y, alpha_prior_sd))
      model_body[[length(model_body)+1]] <- quote(diff ~ normal(0, diff_prior_sd))
      model_body[[length(model_body)+1]] <- quote(sd ~ exponential(sigma_rate_weak))
    }
    view_vars <- c("diff", "delta", "mean", "sd")
    df_pars <- if (use_delta_param) c("mean", "delta") else c("mean", "diff")
  }

  model_ast <- as.call(c(list(as.name("{")), model_body))
  code_obj <- eval(substitute(
    rtmb_code(setup = S, parameters = P, transform = T, model = M),
    list(S = setup_ast, P = param_ast, T = tran_ast, M = model_ast)
  ))

  tmp_env <- list2env(dat)
  ordered_data <- env_to_ordered_list(tmp_env, dat, setup_ast)
  obj <- rtmb_model(data = ordered_data, code = code_obj, fixed = fixed, view = view_vars)

  obj$type <- "ttest"
  obj$extra$df_pars <- df_pars
  obj$extra$levs <- levs
  if (!is.null(null)) obj <- obj$null_model(target = null)

  return(obj)
}

#' Mixture Model Wrapper for RTMB
#'
#' @description
#' Provides a user-friendly interface for fitting Gaussian mixture models
#' with optional covariates on class membership probabilities and various
#' covariance structures.
#'
#' @param formula A formula specifying the response variable(s). For multivariate, use \code{cbind(y1, y2) ~ 1}.
#' @param k Number of mixture components.
#' @param data A data frame containing the variables in the model.
#' @param covariance Covariance structure: "diagonal" (default), "diagonal_equal", "full", "full_equal", or "full_equal_corr".
#' @param prior Prior configuration object: \code{prior_uniform()} (default), \code{prior_weak()}, \code{prior_rhs()}, or \code{prior_ssp()}.
#' @param ... Additional arguments passed to \code{rtmb_model}.
#' @return A \code{RTMB_Model} object.
#' @export
rtmb_mixture <- function(formula, k = 2, data = NULL,
                         covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
                         prior = prior_uniform(), fixed = NULL, ...) {

  if (is.null(prior)) {
    prior <- prior_uniform()
  }

  if (!inherits(prior, "rtmb_prior")) {
    stop("prior must be an object of class 'rtmb_prior'. Use prior_weak(), prior_rhs(), or prior_ssp().")
  }

  # NSE for formula: handle case where formula is just a variable name in data
  formula_expr <- substitute(formula)
  formula_val <- try(formula, silent = TRUE)

  if (inherits(formula_val, "try-error") ||
      (!inherits(formula_val, "formula") && !is.matrix(formula_val) && !is.vector(formula_val))) {
    # Check if variables in formula_expr are in data
    if (!is.null(data) && all(all.vars(formula_expr) %in% names(data))) {
      formula <- as.formula(call("~", formula_expr, 1))
    } else if (!inherits(formula_val, "try-error")) {
      formula <- formula_val
    } else {
      stop(paste0("Object '", deparse(formula_expr), "' not found."))
    }
  } else {
    formula <- formula_val
  }

  covariance <- match.arg(covariance)
  K_mix <- k
  prior_type <- prior$type
  regularization <- if (prior_type %in% c("rhs", "ssp")) prior_type else "none"
  use_weak_info <- prior_type %in% c("weak", "rhs", "ssp")

  default_prior <- list(Intercept_sd = 10, mu_sd = 10, b_sd = 10, sigma_rate = 1, lkj_eta = 1.0, dirichlet_alpha = 1.0)
  prior <- .merge_prior(default_prior, prior)
  # Sync mu_sd and Intercept_sd
  if (!is.null(prior$Intercept_sd) && prior$Intercept_sd != 10) prior$mu_sd <- prior$Intercept_sd
  if (!is.null(prior$mu_sd) && prior$mu_sd != 10) prior$Intercept_sd <- prior$mu_sd

  if (use_weak_info) {
    if (is.null(prior$max_beta)) prior$max_beta <- 1.0
    if (is.null(prior$expected_vars)) prior$expected_vars <- 3
    if (is.null(prior$slab_scale)) prior$slab_scale <- 2.0
    if (is.null(prior$slab_df)) prior$slab_df <- 4.0
    if (is.null(prior$ssp_ratio)) prior$ssp_ratio <- 0.25
  }

  # Parse formulas and prepare data
  if (is.matrix(formula) || is.vector(formula)) {
    Y_mat <- as.matrix(formula)
    X_prob <- NULL
  } else if (!inherits(formula, "formula")) {
    # Case where a single variable name is passed
    var_name <- as.character(formula)
    if (!is.null(data) && var_name %in% names(data)) {
      Y_mat <- as.matrix(data[[var_name]])
    } else {
      Y_mat <- as.matrix(formula)
    }
    X_prob <- NULL
  } else {
    mf_y <- model.frame(formula, data = data)
    Y_mat <- model.response(mf_y)
    if (is.null(Y_mat)) Y_mat <- as.matrix(mf_y)

    # Check if formula has RHS covariates for prob
    # RHS is formula[[3]] for two-sided formula
    rhs <- formula[[3]]
    is_empty_rhs <- (is.numeric(rhs) && rhs == 1) || (is.symbol(rhs) && rhs == "1")

    if (!is_empty_rhs) {
      X_prob <- model.matrix(formula, data = data)
    } else {
      X_prob <- NULL
    }
  }

  Y_mat <- as.matrix(Y_mat)
  N_obs <- nrow(Y_mat)
  P_dim <- ncol(Y_mat)

  if (is.null(colnames(Y_mat))) {
    colnames(Y_mat) <- paste0("Y", 1:P_dim)
  }

  multivariate <- P_dim > 1
  has_cov_prob <- !is.null(X_prob)
  K_prob <- if (has_cov_prob) ncol(X_prob) else 0
  is_sigma_equal <- covariance %in% c("diagonal_equal", "full_equal", "full_equal_corr")

  # --- 1. Dynamic AST Construction: setup ---
  setup_exprs <- list(
    as.name("{"),
    quote(N <- N),
    quote(K <- K),
    quote(P <- P)
  )
  if (has_cov_prob) {
    setup_exprs[[length(setup_exprs) + 1]] <- quote(K_prob <- K_prob)
    setup_exprs[[length(setup_exprs) + 1]] <- quote(X_means <- matrix(X_means, 1, K_prob))
    if (regularization == "rhs") {
      # RHS setup
      setup_exprs[[length(setup_exprs) + 1]] <- quote(tau0 <- tau0)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(half_slab_df <- half_slab_df)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(half_slab_scale2 <- half_slab_scale2)
    } else if (regularization == "ssp") {
      # SSP setup
      setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_scale <- tau_scale)
    }
  }
  setup_ast <- as.call(setup_exprs)

  # --- 2. Dynamic AST Construction: parameters ---
  param_exprs <- list(as.name("{"))
  if (has_cov_prob) {
      if (regularization == "rhs") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(z_b <- Dim(c(.(K_prob), .(K_mix - 1))))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_lambda_b <- Dim(c(.(K_prob), .(K_mix - 1)), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_tau_b <- Dim(.(K_mix - 1), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(c2_b <- Dim(.(K_mix - 1), lower = 0))
      } else if (regularization == "ssp") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b_raw <- Dim(c(.(K_prob), .(K_mix - 1))))
        param_exprs[[length(param_exprs) + 1]] <- bquote(r_b <- Dim(c(.(K_prob), .(K_mix - 1)), lower = 0, upper = 1))
        param_exprs[[length(param_exprs) + 1]] <- bquote(tau_b <- Dim(.(K_mix - 1), lower = 0))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b <- Dim(c(.(K_prob), .(K_mix - 1))))
      }
  } else {
    param_exprs[[length(param_exprs) + 1]] <- bquote(theta <- Dim(.(K_mix), type = "simplex"))
  }

  if (multivariate) {
    param_exprs[[length(param_exprs) + 1]] <- bquote(mu <- Dim(c(K, P)))
    s_dim <- if (is_sigma_equal && covariance != "diagonal_equal") P_dim else c(K_mix, P_dim)
    if (covariance == "diagonal_equal") s_dim <- P_dim
    param_exprs[[length(param_exprs) + 1]] <- bquote(sigma <- Dim(.(s_dim), lower = 0))
  } else {
    param_exprs[[length(param_exprs) + 1]] <- bquote(mu <- Dim(K))
    s_dim <- if (is_sigma_equal) NULL else K_mix
    param_exprs[[length(param_exprs) + 1]] <- if (is.null(s_dim)) quote(sigma <- Dim(lower = 0)) else bquote(sigma <- Dim(.(s_dim), lower = 0))
  }

  if (multivariate && covariance %in% c("full", "full_equal", "full_equal_corr")) {
    if (covariance == "full") {
      param_exprs[[length(param_exprs) + 1]] <- quote(L_corr <- Dim(c(K, P, P), type = "CF_corr"))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- quote(L_corr <- Dim(c(P, P), type = "CF_corr"))
    }
  }
  param_ast <- as.call(param_exprs)

  # --- 3. Dynamic AST Construction: model ---
  model_exprs <- list(as.name("{"))

  # Regularization transformations
  if (has_cov_prob) {
    if (regularization == "rhs") {
      model_exprs[[length(model_exprs) + 1]] <- quote(lambda_b <- 1 / sqrt(w_lambda_b))
      model_exprs[[length(model_exprs) + 1]] <- quote(tau_hs_b <- tau0 / sqrt(w_tau_b))
      model_exprs[[length(model_exprs) + 1]] <- quote(b <- matrix(0, K_prob, K-1)) # Placeholder for correct dim if needed, but b is local
      model_exprs[[length(model_exprs) + 1]] <- quote(for (k in 1:(K-1)) {
        lambda_tilde_k <- sqrt(c2_b[k] * lambda_b[, k]^2 / (c2_b[k] + tau_hs_b[k]^2 * lambda_b[, k]^2))
        b[, k] <- z_b[, k] * tau_hs_b[k] * lambda_tilde_k
      })
    } else if (regularization == "ssp") {
      model_exprs[[length(model_exprs) + 1]] <- quote(b <- b_raw * r_b * tau_b)
    }
  }

  model_exprs[[length(model_exprs) + 1]] <- quote(lp <- mu[1] * 0)
  model_exprs[[length(model_exprs) + 1]] <- quote(log_dens_mat <- matrix(mu[1] * 0, N, K))

  is_diag <- covariance %in% c("diagonal", "diagonal_equal")
  s_k_expr <- if (is_sigma_equal) quote(sigma) else quote(sigma[k, ])
  if (!multivariate && !is_sigma_equal) s_k_expr <- quote(sigma[k])

  dist_body <- if (multivariate) {
    if (is_diag) {
      bquote({
        ld <- Y[1] * 0
        s_vec <- .(s_k_expr)
        m_vec <- mu[k, ]
        for (p in 1:P) {
          ld <- ld + normal_lpdf(Y[, p], m_vec[p], s_vec[p], sum = FALSE)
        }
        log_dens_mat[, k] <- ld
      })
    } else {
      L_corr_k_expr <- if (covariance == "full") quote(matrix(L_corr[k, , ], P, P)) else quote(L_corr)
      bquote({
        log_dens_mat[, k] <- multi_normal_CF_lpdf(Y, mean = mu[k, ], sd = .(s_k_expr), CF_Omega = .(L_corr_k_expr), sum = FALSE)
      })
    }
  } else {
    bquote({
      log_dens_mat[, k] <- normal_lpdf(Y, mu[k], .(s_k_expr), sum = FALSE)
    })
  }

  model_exprs[[length(model_exprs) + 1]] <- as.call(list(as.name("for"), as.name("k"), quote(1:K), dist_body))

  if (has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- quote(eta_prob <- X_prob %*% b)
    model_exprs[[length(model_exprs) + 1]] <- quote(log_pi_mat <- matrix(eta_prob[1] * 0, N, K))
    model_exprs[[length(model_exprs) + 1]] <- quote(for (i in 1:N) {
      log_pi_mat[i, ] <- log_softmax(c(0, eta_prob[i, ]))
    })
    model_exprs[[length(model_exprs) + 1]] <- quote(lp <- lp + sum(log_sum_exp(log_pi_mat + log_dens_mat)))
  } else {
    model_exprs[[length(model_exprs) + 1]] <- quote(log_theta <- log(theta))
    model_exprs[[length(model_exprs) + 1]] <- quote(lp <- lp + sum(log_sum_exp(t(t(log_dens_mat) + log_theta))))
  }

  # Priors
  if (use_weak_info) {
    # Component priors
    model_exprs[[length(model_exprs) + 1]] <- bquote(mu ~ normal(0, .(prior$Intercept_sd)))
    model_exprs[[length(model_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))

    if (multivariate && covariance %in% c("full", "full_equal", "full_equal_corr")) {
      if (covariance == "full") {
        model_exprs[[length(model_exprs) + 1]] <- bquote(for (k in 1:K) matrix(L_corr[k, , ], P, P) ~ lkj_CF_corr(.(prior$lkj_eta)))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(L_corr ~ lkj_CF_corr(.(prior$lkj_eta)))
      }
    }

    if (has_cov_prob) {
      if (regularization == "rhs") {
        model_exprs[[length(model_exprs) + 1]] <- quote(z_b ~ normal(0, 1))
        model_exprs[[length(model_exprs) + 1]] <- quote(w_lambda_b ~ gamma(0.5, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- quote(w_tau_b ~ gamma(0.5, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- quote(c2_b ~ inverse_gamma(half_slab_df, half_slab_scale2))
      } else if (regularization == "ssp") {
        model_exprs[[length(model_exprs) + 1]] <- quote(b_raw ~ laplace(0, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- bquote(mu_r <- log(.(prior$ssp_ratio) / (1 - .(prior$ssp_ratio))))
        model_exprs[[length(model_exprs) + 1]] <- quote(r_b ~ logit_normal(mu_r, 3))
        model_exprs[[length(model_exprs) + 1]] <- quote(tau_b ~ exponential(1 / tau_scale))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
      }
    } else {
      model_exprs[[length(model_exprs) + 1]] <- bquote(theta ~ dirichlet(rep(.(prior$dirichlet_alpha), K)))
    }
  } else if (prior_type == "uniform") {
    # Default is uniform (flat), but if users specified specific sd/rates in prior object, apply them
    if (!is.null(prior$Intercept_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(mu ~ normal(0, .(prior$Intercept_sd)))
    if (!is.null(prior$sigma_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    if (has_cov_prob && !is.null(prior$b_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
  }

  model_ast <- as.call(model_exprs)

  # --- 4. Dynamic AST Construction: generate ---
  generate_exprs <- list(as.name("{"))
  if (has_cov_prob) {
    if (regularization == "rhs") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(lambda_b <- 1 / sqrt(w_lambda_b))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(tau_hs_b <- tau0 / sqrt(w_tau_b))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- matrix(0, P, K-1))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:(K-1)) {
        lambda_tilde_k <- sqrt(c2_b[k] * lambda_b[, k]^2 / (c2_b[k] + tau_hs_b[k]^2 * lambda_b[, k]^2))
        b[, k] <- z_b[, k] * tau_hs_b[k] * lambda_tilde_k
      })
    } else if (regularization == "ssp") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- b_raw * r_b * tau_b)
    }
    generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- softmax(c(0, X_means %*% b)))
  } else {
    generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- theta)
  }

  if (multivariate && !is_diag) {
    if (covariance == "full") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(corr <- array(mu[1] * 0, dim = c(K, P, P)))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:K) {
        corr[k, , ] <- matrix(L_corr[k, , ], P, P) %*% t(matrix(L_corr[k, , ], P, P))
      })
    } else {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(corr <- L_corr %*% t(L_corr))
    }
  }
  generate_ast <- as.call(generate_exprs)

  mdl_code <- list(
    setup = setup_ast,
    parameters = param_ast,
    model = model_ast,
    generate = generate_ast,
    env = parent.frame()
  )
  class(mdl_code) <- "rtmb_code"

  v_names <- list()
  if (multivariate) {
    v_names$mu <- list(paste0("C", 1:K_mix), colnames(Y_mat))
    v_names$sigma <- if (is_sigma_equal && covariance != "diagonal_equal") colnames(Y_mat) else list(paste0("C", 1:K_mix), colnames(Y_mat))
    if (covariance == "diagonal_equal") v_names$sigma <- colnames(Y_mat)
  } else {
    v_names$mu <- paste0("C", 1:K_mix)
    v_names$sigma <- if (is_sigma_equal) "sigma" else paste0("C", 1:K_mix)
  }

  if (has_cov_prob) {
    colnames(X_prob)[colnames(X_prob) == "(Intercept)"] <- "Intercept"
    v_names$b <- list(colnames(X_prob), paste0("C", 2:K_mix))
  }
  if (multivariate && !is_diag) {
    if (covariance == "full") {
      v_names$L_corr <- list(paste0("C", 1:K_mix), colnames(Y_mat), colnames(Y_mat))
      v_names$corr <- list(paste0("C", 1:K_mix), colnames(Y_mat), colnames(Y_mat))
    } else {
      v_names$L_corr <- list(colnames(Y_mat), colnames(Y_mat))
      v_names$corr <- list(colnames(Y_mat), colnames(Y_mat))
    }
  }

  data_list <- list(Y = Y_mat, N = N_obs, K = K_mix, P = P_dim)
  if (has_cov_prob) {
    data_list$X_prob <- X_prob
    data_list$K_prob <- K_prob
    data_list$X_means <- as.vector(colMeans(X_prob))
    if (regularization == "rhs") {
      p0 <- prior$expected_vars
      data_list$tau0 <- p0 / (K_prob - p0) / sqrt(N_obs)
      data_list$half_slab_df <- prior$slab_df / 2
      data_list$half_slab_scale2 <- prior$slab_scale^2 / 2
    } else if (regularization == "ssp") {
      data_list$tau_scale <- prior$max_beta / 1.96
    }
  }

  init_list <- list()
  if (nrow(Y_mat) >= K_mix) {
    # Use kmeans for better initial values
    km <- try(kmeans(Y_mat, centers = K_mix, nstart = 5), silent = TRUE)
    if (!inherits(km, "try-error")) {
      if (multivariate) {
        init_list$mu <- km$centers
      } else {
        init_list$mu <- as.vector(km$centers)
      }
      # Initial theta / probabilities
      props <- as.vector(table(factor(km$cluster, levels = 1:K_mix))) / nrow(Y_mat)
      init_list$theta <- props

      # Initial b if covariates present
      if (has_cov_prob) {
        b_init <- matrix(0, ncol(X_prob), K_mix - 1)
        for (k in 2:K_mix) {
          z_k <- as.numeric(km$cluster == k)
          fit_k <- try(suppressWarnings(glm(z_k ~ X_prob - 1 + offset(rep(qlogis(max(0.01, props[k])), length(z_k))), family = binomial)), silent = TRUE)
          if (!inherits(fit_k, "try-error")) {
            b_init[, k-1] <- coef(fit_k)
          }
        }
        if (regularization == "rhs") {
          init_list$z_b <- b_init
          init_list$w_lambda_b <- matrix(1, nrow(b_init), ncol(b_init))
          init_list$w_tau_b <- rep(1, ncol(b_init))
          init_list$c2_b <- rep(1, ncol(b_init))
        } else if (regularization == "ssp") {
          init_list$b_raw <- b_init
          init_list$r_b <- matrix(0.5, nrow(b_init), ncol(b_init))
          init_list$tau_b <- rep(1, ncol(b_init))
        } else {
          init_list$b <- b_init
        }
      }
    }
  }

  if (is.null(init_list$mu)) {
    # Fallback to even spread
    if (multivariate) {
      mu_init <- matrix(0, K_mix, P_dim)
      for (p in 1:P_dim) {
        r_y <- range(Y_mat[, p])
        mu_init[, p] <- seq(r_y[1], r_y[2], length.out = K_mix + 2)[2:(K_mix + 1)]
      }
      init_list$mu <- mu_init
    } else {
      r_y <- range(Y_mat)
      init_list$mu <- seq(r_y[1], r_y[2], length.out = K_mix + 2)[2:(K_mix + 1)]
    }
  }

  if (!is.null(Y_mat)) {
     init_list$sigma <- if (is_sigma_equal) apply(Y_mat, 2, sd) else matrix(apply(Y_mat, 2, sd), K_mix, P_dim, byrow = TRUE)
  }

  view_order <- if (has_cov_prob) c("b", "prob_mean", "mu", "sigma") else c("prob_mean", "mu", "sigma")
  if (multivariate && !is_diag) view_order <- c(view_order, "corr")

  mdl <- rtmb_model(data_list, mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order)
  return(mdl)
}

#' Fit a Latent Rank Theory (LRT) Model
#'
#' @description
#' Fits a Latent Rank Theory model, which is a mixture model with ordered ranks
#' and Gaussian Process smoothing on the mean profiles.
#'
#' @param formula A formula specifying the response variable(s).
#' @param k Number of ranks (mixture components).
#' @param data A data frame containing the variables.
#' @param magnitude Signal standard deviation for the GP prior. If NULL, it is estimated.
#' @param smoothing Length-scale for the GP prior. If NULL, it is estimated.
#' @param noise Measurement noise for the GP prior (default is 0.01).
#' @param prior Prior configuration object: \code{prior_uniform()} (default), \code{prior_weak()}, \code{prior_rhs()}, or \code{prior_ssp()}.
#' @param ... Additional arguments passed to \code{rtmb_model}.
#' @return A \code{RTMB_Model} object.
#' @export
rtmb_lrt <- function(formula, k = 3, data = NULL,
                     rank_coords = NULL,
                     covariance = c("diagonal", "diagonal_equal", "full", "full_equal", "full_equal_corr"),
                     magnitude = NULL, smoothing = NULL, noise = 0.01,
                     prob_smoothing = FALSE,
                     link = c("ordered", "sequential"),
                     prior = prior_uniform(), fixed = NULL, ...) {

  if (is.null(rank_coords)) rank_coords <- 1:k

  if (is.null(prior)) {
    prior <- prior_uniform()
  }

  if (!inherits(prior, "rtmb_prior")) {
    stop("prior must be an object of class 'rtmb_prior'. Use prior_weak(), prior_rhs(), or prior_ssp().")
  }

  # NSE for formula: handle case where formula is just a variable name in data
  formula_expr <- substitute(formula)
  formula_val <- try(formula, silent = TRUE)

  if (inherits(formula_val, "try-error") ||
      (!inherits(formula_val, "formula") && !is.matrix(formula_val) && !is.vector(formula_val))) {
    if (!is.null(data) && all(all.vars(formula_expr) %in% names(data))) {
      formula <- as.formula(call("~", formula_expr, 1))
    } else if (!inherits(formula_val, "try-error")) {
      formula <- formula_val
    } else {
      stop(paste0("Object '", deparse(formula_expr), "' not found."))
    }
  } else {
    formula <- formula_val
  }

  covariance <- match.arg(covariance)
  link <- match.arg(link)
  K_mix <- k
  prior_type <- prior$type
  regularization <- if (prior_type %in% c("rhs", "ssp")) prior_type else "none"
  use_weak_info <- prior_type %in% c("weak", "rhs", "ssp")

  default_prior <- list(Intercept_sd = 10, mu_sd = 10, b_sd = 10, sigma_rate = 1, lkj_eta = 1.0, mag_rate = 1.0, smooth_rate = 1.0, cutpoint_sd = 2.5)
  prior <- .merge_prior(default_prior, prior)
  # Sync mu_sd and Intercept_sd
  if (!is.null(prior$Intercept_sd) && prior$Intercept_sd != 10) prior$mu_sd <- prior$Intercept_sd
  if (!is.null(prior$mu_sd) && prior$mu_sd != 10) prior$Intercept_sd <- prior$mu_sd

  if (use_weak_info) {
    if (is.null(prior$max_beta)) prior$max_beta <- 1.0
    if (is.null(prior$expected_vars)) prior$expected_vars <- 3
    if (is.null(prior$slab_scale)) prior$slab_scale <- 2.0
    if (is.null(prior$slab_df)) prior$slab_df <- 4.0
    if (is.null(prior$ssp_ratio)) prior$ssp_ratio <- 0.25
  }

  # Data parsing
  if (is.matrix(formula) || is.vector(formula)) {
    Y_mat <- as.matrix(formula)
    X_prob <- NULL
  } else if (!inherits(formula, "formula")) {
    var_name <- as.character(formula)
    if (!is.null(data) && var_name %in% names(data)) {
      Y_mat <- as.matrix(data[[var_name]])
    } else {
      Y_mat <- as.matrix(formula)
    }
    X_prob <- NULL
  } else {
    mf_y <- model.frame(formula, data = data)
    Y_mat <- model.response(mf_y)
    if (is.null(Y_mat)) Y_mat <- as.matrix(mf_y)

    # Check if formula has RHS covariates for prob
    rhs <- formula[[3]]
    is_empty_rhs <- (is.numeric(rhs) && rhs == 1) || (is.symbol(rhs) && rhs == "1")

    if (!is_empty_rhs) {
      X_prob <- model.matrix(formula, data = data)
      if ("(Intercept)" %in% colnames(X_prob)) {
        X_prob <- X_prob[, colnames(X_prob) != "(Intercept)", drop = FALSE]
      }
      if (ncol(X_prob) == 0) X_prob <- NULL
    } else {
      X_prob <- NULL
    }
  }

  Y_mat <- as.matrix(Y_mat)
  N_obs <- nrow(Y_mat)
  P_dim <- ncol(Y_mat)
  if (is.null(colnames(Y_mat))) colnames(Y_mat) <- paste0("Y", 1:P_dim)
  has_cov_prob <- !is.null(X_prob)
  K_prob <- if (has_cov_prob) ncol(X_prob) else 0

  multivariate <- P_dim > 1
  is_diag <- covariance %in% c("diagonal", "diagonal_equal")
  is_sigma_equal <- covariance %in% c("diagonal_equal", "full_equal", "full_equal_corr")

  # --- 1. Setup ---
  setup_exprs <- list(
    as.name("{"),
    quote(N <- N),
    quote(K <- K),
    quote(P <- P),
    quote(rank_coords <- rank_coords)
  )
  if (has_cov_prob) {
    setup_exprs[[length(setup_exprs) + 1]] <- quote(K_prob <- K_prob)
    setup_exprs[[length(setup_exprs) + 1]] <- quote(X_means <- matrix(X_means, 1, K_prob))
    if (regularization == "rhs") {
      # RHS setup
      setup_exprs[[length(setup_exprs) + 1]] <- quote(tau0 <- tau0)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(half_slab_df <- half_slab_df)
      setup_exprs[[length(setup_exprs) + 1]] <- quote(half_slab_scale2 <- half_slab_scale2)
    } else if (regularization == "ssp") {
      # SSP setup
      setup_exprs[[length(setup_exprs) + 1]] <- quote(tau_scale <- tau_scale)
    }
  }
  setup_ast <- as.call(setup_exprs)

  # --- 2. Parameters ---
  param_exprs <- list(as.name("{"))
  param_exprs[[length(param_exprs) + 1]] <- bquote(mu_p <- Dim(c(.(P_dim), .(K_mix)), type = "ordered"))

  # Sigma parameterization based on covariance
  if (is_sigma_equal) {
    param_exprs[[length(param_exprs) + 1]] <- bquote(sigma <- Dim(.(P_dim), lower = 0))
  } else {
    param_exprs[[length(param_exprs) + 1]] <- bquote(sigma <- Dim(c(.(K_mix), .(P_dim)), lower = 0))
  }

  if (multivariate && !is_diag) {
    if (covariance == "full") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr <- Dim(c(.(K_mix), .(P_dim), .(P_dim)), type = "CF_corr"))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- bquote(L_corr <- Dim(c(.(P_dim), .(P_dim)), type = "CF_corr"))
    }
  }

  if (has_cov_prob) {
    if (link == "ordered") {
      if (regularization == "rhs") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(z_b <- Dim(.(K_prob)))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_lambda_b <- Dim(.(K_prob), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_tau_b <- Dim(lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(c2_b <- Dim(lower = 0))
      } else if (regularization == "ssp") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b_raw <- Dim(.(K_prob)))
        param_exprs[[length(param_exprs) + 1]] <- bquote(r_b <- Dim(.(K_prob), lower = 0, upper = 1))
        param_exprs[[length(param_exprs) + 1]] <- bquote(tau_b <- Dim(lower = 0))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b <- Dim(.(K_prob)))
      }
      param_exprs[[length(param_exprs) + 1]] <- bquote(cutpoints <- Dim(.(K_mix - 1), type = "ordered"))
    } else {
      # Sequential: Separate b for each step
      if (regularization == "rhs") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(z_b <- Dim(c(.(K_prob), .(K_mix - 1))))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_lambda_b <- Dim(c(.(K_prob), .(K_mix - 1)), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(w_tau_b <- Dim(.(K_mix - 1), lower = 0))
        param_exprs[[length(param_exprs) + 1]] <- bquote(c2_b <- Dim(.(K_mix - 1), lower = 0))
      } else if (regularization == "ssp") {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b_raw <- Dim(c(.(K_prob), .(K_mix - 1))))
        param_exprs[[length(param_exprs) + 1]] <- bquote(r_b <- Dim(c(.(K_prob), .(K_mix - 1)), lower = 0, upper = 1))
        param_exprs[[length(param_exprs) + 1]] <- bquote(tau_b <- Dim(.(K_mix - 1), lower = 0))
      } else {
        param_exprs[[length(param_exprs) + 1]] <- bquote(b <- Dim(c(.(K_prob), .(K_mix - 1))))
      }
      param_exprs[[length(param_exprs) + 1]] <- bquote(alpha <- Dim(.(K_mix - 1)))
    }
  } else {
    if (prob_smoothing) {
      param_exprs[[length(param_exprs) + 1]] <- bquote(logit_theta <- Dim(.(K_mix)))
      param_exprs[[length(param_exprs) + 1]] <- quote(magnitude_theta <- Dim(lower = 0))
    } else {
      param_exprs[[length(param_exprs) + 1]] <- bquote(theta <- Dim(.(K_mix), type = "simplex"))
    }
  }

  if (is.null(magnitude)) param_exprs[[length(param_exprs) + 1]] <- quote(magnitude <- Dim(lower = 0))
  if (is.null(smoothing)) param_exprs[[length(param_exprs) + 1]] <- quote(smoothing <- Dim(lower = 0))
  param_ast <- as.call(param_exprs)

  # --- 3. Model ---
  model_exprs <- list(as.name("{"))

  # Regularization transformations
  if (has_cov_prob) {
    if (regularization == "rhs") {
      model_exprs[[length(model_exprs) + 1]] <- quote(lambda_b <- 1 / sqrt(w_lambda_b))
      model_exprs[[length(model_exprs) + 1]] <- quote(tau_hs_b <- tau0 / sqrt(w_tau_b))
      if (link == "ordered") {
        model_exprs[[length(model_exprs) + 1]] <- quote(lambda_tilde_b <- sqrt(c2_b * lambda_b^2 / (c2_b + tau_hs_b^2 * lambda_b^2)))
        model_exprs[[length(model_exprs) + 1]] <- quote(b <- z_b * tau_hs_b * lambda_tilde_b)
      } else {
        model_exprs[[length(model_exprs) + 1]] <- quote(b <- matrix(0, K_prob, K-1))
        model_exprs[[length(model_exprs) + 1]] <- quote(for (k in 1:(K-1)) {
          lambda_tilde_k <- sqrt(c2_b[k] * lambda_b[, k]^2 / (c2_b[k] + tau_hs_b[k]^2 * lambda_b[, k]^2))
          b[, k] <- z_b[, k] * tau_hs_b[k] * lambda_tilde_k
        })
      }
    } else if (regularization == "ssp") {
      model_exprs[[length(model_exprs) + 1]] <- quote(b <- b_raw * r_b * tau_b)
    }
  }

  model_exprs[[length(model_exprs) + 1]] <- quote(mu <- t(mu_p))

  if (has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- quote(pi_mat <- matrix(mu[1] * 0, N, K))
    if (link == "ordered") {
      model_exprs[[length(model_exprs) + 1]] <- quote(eta <- X_prob %*% b)
      model_exprs[[length(model_exprs) + 1]] <- quote(for (i in 1:N) {
        F_prev <- 0
        for (k in 1:(K-1)) {
          F_k <- inv_logit(cutpoints[k] - eta[i])
          pi_mat[i, k] <- F_k - F_prev
          F_prev <- F_k
        }
        pi_mat[i, K] <- 1 - F_prev
      })
    } else {
      # Sequential with Step-specific b
      model_exprs[[length(model_exprs) + 1]] <- quote(eta_mat <- X_prob %*% b)
      model_exprs[[length(model_exprs) + 1]] <- quote(for (i in 1:N) {
        q_cum <- 1
        for (k in 1:(K-1)) {
          q_k <- inv_logit(alpha[k] + eta_mat[i, k])
          pi_mat[i, k] <- q_cum * (1 - q_k)
          q_cum <- q_cum * q_k
        }
        pi_mat[i, K] <- q_cum
      })
    }
    model_exprs[[length(model_exprs) + 1]] <- quote(log_pi_mat <- log(pi_mat + 1e-15))
  } else if (prob_smoothing) {
    model_exprs[[length(model_exprs) + 1]] <- quote(theta <- softmax(logit_theta))
  }

  model_exprs[[length(model_exprs) + 1]] <- quote(lp <- mu[1] * 0)
  model_exprs[[length(model_exprs) + 1]] <- quote(log_dens_mat <- matrix(mu[1] * 0, N, K))

  # Likelihood Calculation
  if (is_diag) {
    dist_body <- bquote({
      ld <- Y[1] * 0
      m_vec <- mu[k, ]
      s_vec <- if (.(is_sigma_equal)) sigma else sigma[k, ]
      for (p in 1:P) {
        ld <- ld + normal_lpdf(Y[, p], m_vec[p], s_vec[p], sum = FALSE)
      }
      log_dens_mat[, k] <- ld
    })
  } else {
    # Full covariance
    dist_body <- bquote({
      m_vec <- mu[k, ]
      s_vec <- if (.(is_sigma_equal)) sigma else sigma[k, ]
      L_mat <- if (.(covariance == "full")) matrix(L_corr[k, , ], P, P) else matrix(L_corr, P, P)
      log_dens_mat[, k] <- multi_normal_CF_lpdf(Y, mean = as.vector(m_vec), sd = as.vector(s_vec), CF_Omega = L_mat, sum = FALSE)
    })
  }
  model_exprs[[length(model_exprs) + 1]] <- as.call(list(as.name("for"), as.name("k"), quote(1:K), dist_body))

  if (has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- quote(lp <- lp + sum(log_sum_exp(log_pi_mat + log_dens_mat)))
  } else {
    model_exprs[[length(model_exprs) + 1]] <- quote(log_theta <- log(theta))
    model_exprs[[length(model_exprs) + 1]] <- quote(lp <- lp + sum(log_sum_exp(t(t(log_dens_mat) + log_theta))))
  }

  # GP Prior (Items)
  gp_mag <- if (is.null(magnitude)) quote(magnitude) else magnitude
  gp_sm  <- if (is.null(smoothing)) quote(smoothing) else smoothing
  model_exprs[[length(model_exprs) + 1]] <- bquote(for (p in 1:P) {
    lp <- lp + gaussian_process_lpdf(mu[, p], x = rank_coords, magnitude = .(gp_mag), smoothing = .(gp_sm), noise = .(noise))
  })

  if (prob_smoothing && !has_cov_prob) {
    model_exprs[[length(model_exprs) + 1]] <- bquote(lp <- lp + gaussian_process_lpdf(logit_theta, x = rank_coords, magnitude = magnitude_theta, smoothing = .(gp_sm), noise = .(noise)))
  }

  # Priors
  if (use_weak_info) {
    model_exprs[[length(model_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    if (is.null(magnitude)) model_exprs[[length(model_exprs) + 1]] <- bquote(magnitude ~ exponential(.(prior$mag_rate)))
    if (is.null(smoothing)) model_exprs[[length(model_exprs) + 1]] <- bquote(smoothing ~ exponential(.(prior$smooth_rate)))
    if (prob_smoothing && !has_cov_prob) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(magnitude_theta ~ exponential(.(prior$mag_rate)))
    }

    if (multivariate && !is_diag) {
      if (covariance == "full") {
        model_exprs[[length(model_exprs) + 1]] <- bquote(for (k in 1:K) matrix(L_corr[k, , ], P, P) ~ lkj_CF_corr(.(prior$lkj_eta)))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(L_corr ~ lkj_CF_corr(.(prior$lkj_eta)))
      }
    }

    if (has_cov_prob) {
      if (regularization == "rhs") {
        model_exprs[[length(model_exprs) + 1]] <- quote(z_b ~ normal(0, 1))
        model_exprs[[length(model_exprs) + 1]] <- quote(w_lambda_b ~ gamma(0.5, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- quote(w_tau_b ~ gamma(0.5, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- quote(c2_b ~ inverse_gamma(half_slab_df, half_slab_scale2))
      } else if (regularization == "ssp") {
        model_exprs[[length(model_exprs) + 1]] <- quote(b_raw ~ laplace(0, 0.5))
        model_exprs[[length(model_exprs) + 1]] <- bquote(mu_r <- log(.(prior$ssp_ratio) / (1 - .(prior$ssp_ratio))))
        model_exprs[[length(model_exprs) + 1]] <- quote(r_b ~ logit_normal(mu_r, 3))
        model_exprs[[length(model_exprs) + 1]] <- quote(tau_b ~ exponential(1 / tau_scale))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
      }

      if (link == "ordered") {
        model_exprs[[length(model_exprs) + 1]] <- bquote(cutpoints ~ normal(0, .(prior$cutpoint_sd)))
      } else {
        model_exprs[[length(model_exprs) + 1]] <- bquote(alpha ~ normal(0, .(prior$cutpoint_sd)))
      }
    }
    if (!has_cov_prob && !prob_smoothing) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(theta ~ dirichlet(rep(1, K)))
    }
  } else if (prior_type == "uniform") {
    if (!is.null(prior$sigma_rate)) model_exprs[[length(model_exprs) + 1]] <- bquote(sigma ~ exponential(.(prior$sigma_rate)))
    if (has_cov_prob && !is.null(prior$b_sd)) model_exprs[[length(model_exprs) + 1]] <- bquote(b ~ normal(0, .(prior$b_sd)))
  }
  model_ast <- as.call(model_exprs)

  # --- 4. Generate ---
  generate_exprs <- list(as.name("{"))
  generate_exprs[[length(generate_exprs) + 1]] <- quote(mu <- t(mu_p))
  if (has_cov_prob) {
    if (regularization == "rhs") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(lambda_b <- 1 / sqrt(w_lambda_b))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(tau_hs_b <- tau0 / sqrt(w_tau_b))
      if (link == "ordered") {
        generate_exprs[[length(generate_exprs) + 1]] <- quote(lambda_tilde_b <- sqrt(c2_b * lambda_b^2 / (c2_b + tau_hs_b^2 * lambda_b^2)))
        generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- z_b * tau_hs_b * lambda_tilde_b)
      } else {
        generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- matrix(0, K_prob, K-1))
        generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:(K-1)) {
          lambda_tilde_k <- sqrt(c2_b[k] * lambda_b[, k]^2 / (c2_b[k] + tau_hs_b[k]^2 * lambda_b[, k]^2))
          b[, k] <- z_b[, k] * tau_hs_b[k] * lambda_tilde_k
        })
      }
    } else if (regularization == "ssp") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(b <- b_raw * r_b * tau_b)
    }
    if (link == "ordered") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(eta_mean <- X_means %*% b)
      generate_exprs[[length(generate_exprs) + 1]] <- quote(F_prev <- 0)
      generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- numeric(K))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:(K - 1)) {
        F_k <- inv_logit(cutpoints[k] - eta_mean)
        prob_mean[k] <- F_k - F_prev
        F_prev <- F_k
      })
      generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean[K] <- 1 - F_prev)
    } else {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(eta_mean_mat <- X_means %*% b)
      generate_exprs[[length(generate_exprs) + 1]] <- quote(q_cum <- 1)
      generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- numeric(K))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:(K - 1)) {
        q_k <- inv_logit(alpha[k] + eta_mean_mat[k])
        prob_mean[k] <- q_cum * (1 - q_k)
        q_cum <- q_cum * q_k
      })
      generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean[K] <- q_cum)
    }
  } else {
    if (prob_smoothing) generate_exprs[[length(generate_exprs) + 1]] <- quote(theta <- softmax(logit_theta))
    generate_exprs[[length(generate_exprs) + 1]] <- quote(prob_mean <- theta)
  }

  if (!is_diag) {
    generate_exprs[[length(generate_exprs) + 1]] <- bquote(corr <- if (.(covariance == "full")) array(mu[1] * 0, dim = c(K, P, P)) else matrix(mu[1] * 0, P, P))
    if (covariance == "full") {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(for (k in 1:K) {
        L_mat <- matrix(L_corr[k, , ], P, P)
        corr[k, , ] <- L_mat %*% t(L_mat)
      })
    } else {
      generate_exprs[[length(generate_exprs) + 1]] <- quote(L_mat <- matrix(L_corr, P, P))
      generate_exprs[[length(generate_exprs) + 1]] <- quote(corr <- L_mat %*% t(L_mat))
    }
  }
  generate_ast <- as.call(generate_exprs)

  mdl_code <- list(setup = setup_ast, parameters = param_ast, model = model_ast, generate = generate_ast, env = parent.frame())
  class(mdl_code) <- "rtmb_code"

  v_names <- list(mu = list(paste0("Rank", 1:K_mix), colnames(Y_mat)))
  v_names$sigma <- if (is_sigma_equal) colnames(Y_mat) else list(paste0("Rank", 1:K_mix), colnames(Y_mat))

  if (!is_diag) {
    if (covariance == "full") {
      v_names$L_corr <- list(paste0("Rank", 1:K_mix), colnames(Y_mat), colnames(Y_mat))
      v_names$corr <- list(paste0("Rank", 1:K_mix), colnames(Y_mat), colnames(Y_mat))
    } else {
      v_names$L_corr <- list(colnames(Y_mat), colnames(Y_mat))
      v_names$corr <- list(colnames(Y_mat), colnames(Y_mat))
    }
  }

  if (has_cov_prob) {
    if (link == "ordered") {
      v_names$b <- colnames(X_prob)
      v_names$cutpoints <- paste0("C", 1:(K_mix - 1))
    } else {
      step_labels <- paste0(1:(K_mix - 1), "->", 2:K_mix)
      v_names$b <- list(colnames(X_prob), step_labels)
      v_names$alpha <- step_labels
    }
  }

  data_list <- list(Y = Y_mat, N = N_obs, K = K_mix, P = P_dim, rank_coords = rank_coords)
  if (has_cov_prob) {
    data_list$X_prob <- X_prob
    data_list$K_prob <- K_prob
    data_list$X_means <- as.vector(colMeans(X_prob))
    if (regularization == "rhs") {
      p0 <- prior$expected_vars
      data_list$tau0 <- p0 / (K_prob - p0) / sqrt(N_obs)
      data_list$half_slab_df <- prior$slab_df / 2
      data_list$half_slab_scale2 <- prior$slab_scale^2 / 2
    } else if (regularization == "ssp") {
      data_list$tau_scale <- prior$max_beta / 1.96
    }
  }

  init_list <- list()
  if (nrow(Y_mat) >= K_mix) {
    # PCA-based ordering for LRT initial values
    # We use PC1 to capture the dominant direction of variation
    pc1 <- try(prcomp(Y_mat, rank. = 1)$x[, 1], silent = TRUE)
    if (!inherits(pc1, "try-error")) {
      # Split by quantiles of PC1
      groups <- cut(pc1, breaks = quantile(pc1, probs = seq(0, 1, length.out = K_mix + 1)), include.lowest = TRUE, labels = FALSE)
      mu_init_p <- matrix(0, P_dim, K_mix)
      for (k in 1:K_mix) {
        idx <- which(groups == k)
        if (length(idx) > 0) {
          mu_init_p[, k] <- colMeans(Y_mat[idx, , drop = FALSE])
        } else {
          mu_init_p[, k] <- colMeans(Y_mat)
        }
      }
      # Initial mu_p: handle potential NAs from colMeans
      mu_init_p[is.na(mu_init_p)] <- 0
      # Ensure it's ordered along the average direction
      if (mean(mu_init_p[, 1]) > mean(mu_init_p[, K_mix])) {
        mu_init_p <- mu_init_p[, K_mix:1]
      }
      # Add small jitter to ensure strict ordering for AD setup
      # AND handle any individual items that might be decreasing
      for (k in 2:K_mix) {
        for (p in 1:P_dim) {
          if (mu_init_p[p, k] <= mu_init_p[p, k-1]) {
            mu_init_p[p, k] <- mu_init_p[p, k-1] + 1e-3
          }
        }
      }
      init_list$mu_p <- mu_init_p

      # Initial b and cutpoints if covariates present
      if (has_cov_prob) {
        # Proxy regression of ranks on X
        fit_b <- try(lm(groups ~ X_prob), silent = TRUE)
        if (!inherits(fit_b, "try-error")) {
          b_est <- coef(fit_b)
          # Exclude intercept if it was included (lm usually does)
          if ("(Intercept)" %in% names(b_est)) b_est <- b_est[names(b_est) != "(Intercept)"]
          b_est[is.na(b_est)] <- 0

          if (regularization == "rhs") {
            init_list$z_b <- if (link == "ordered") as.vector(b_est) else matrix(b_est, ncol(X_prob), K_mix - 1)
            init_list$w_lambda_b <- if (link == "ordered") rep(1, length(init_list$z_b)) else matrix(1, ncol(X_prob), K_mix - 1)
            init_list$w_tau_b <- if (link == "ordered") 1 else rep(1, K_mix - 1)
            init_list$c2_b <- if (link == "ordered") 1 else rep(1, K_mix - 1)
          } else if (regularization == "ssp") {
            init_list$b_raw <- if (link == "ordered") as.vector(b_est) else matrix(b_est, ncol(X_prob), K_mix - 1)
            init_list$r_b <- if (link == "ordered") rep(0.5, length(init_list$b_raw)) else matrix(0.5, ncol(X_prob), K_mix - 1)
            init_list$tau_b <- if (link == "ordered") 1 else rep(1, K_mix - 1)
          } else {
            if (link == "ordered") {
              init_list$b <- as.vector(b_est)
            } else {
              init_list$b <- matrix(b_est, ncol(X_prob), K_mix - 1)
            }
          }
        }

        # Thresholds from group proportions
        props <- table(factor(groups, levels = 1:K_mix)) / length(groups)
        F_k <- cumsum(props)
        jitter <- (1:(K_mix - 1)) * 1e-4
        if (link == "ordered") {
          init_list$cutpoints <- qlogis(pmin(pmax(F_k[1:(K_mix-1)], 0.01), 0.99)) + jitter
        } else {
          # Sequential alphas: q_k = P(Z>k|Z>=k)
          q_vals <- (1 - F_k[1:(K_mix-1)]) / (1 - c(0, F_k)[1:(K_mix-1)])
          init_list$alpha <- qlogis(pmin(pmax(q_vals, 0.01), 0.99)) + jitter
        }
      }
    }
  }

  if (is.null(init_list$mu_p)) {
    mu_init <- matrix(0, P_dim, K_mix)
    for (p in 1:P_dim) {
      r_y <- range(Y_mat[, p])
      mu_init[p, ] <- seq(r_y[1], r_y[2], length.out = K_mix)
    }
    init_list$mu_p <- mu_init
  }

  if (!is.null(Y_mat)) {
    s_init <- apply(Y_mat, 2, sd)
    if (is_sigma_equal) {
      init_list$sigma <- s_init
    } else {
      init_list$sigma <- matrix(s_init, K_mix, P_dim, byrow = TRUE)
    }
  }

  view_order <- c("prob_mean", "mu", "sigma")
  if (has_cov_prob) {
    view_order <- if (link == "ordered") c("b", "cutpoints", view_order) else c("b", "alpha", view_order)
  }
  if (!is_diag) view_order <- c(view_order, "corr")

  mdl <- rtmb_model(data_list, mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order)
  return(mdl)
}

#' Multilevel Correlation Analysis (Kenny's Model)
#'
#' @description
#' This function fits a multilevel correlation model to estimate between-group and within-group
#' covariance/correlation structures, along with Intraclass Correlation Coefficients (ICC).


#' RTMB-based Mediation Analysis Wrapper
#'
#' @description
#' `rtmb_mediation` performs mediation analysis by simultaneously estimating multiple
#' GLM regression equations. It automatically identifies mediation paths and calculates
#' indirect, direct, and total effects.
#'
#' @param formula A list of formulas defining the regression paths (e.g., `list(M ~ X, Y ~ X + M)`).
#' @param data A data frame containing the variables.
#' @param family A single character string or a list of character strings specifying the error distribution
#'   for each equation (e.g., `family = list("gaussian", "binomial")`). Default is "gaussian".
#' @param prior An object of class "rtmb_prior" specifying the prior distribution.
#' @param y_range Theoretical minimum and maximum values of the response variable.
#' @param fixed A named list of parameter values to fix (optional).
#' @param view Character vector of parameter names to prioritize in summary.

#' @param ... Additional arguments passed to the model construction.
#'
#' @details
#' The function identifies mediation paths by looking for variables that are
#' responses in one equation and predictors in another. Indirect effects are
#' calculated as the product of coefficients along these paths ($a \times b$).
#'
#' \strong{Uncertainty Estimation}:
#' When using `$optimize(ci_method = "sampling")`, the function provides asymmetric
#' confidence intervals for indirect effects based on the distribution of products,
#' which is more accurate than the standard Sobel test (Delta Method).
#'
#' @return An `RTMB_Model` object.
#' @export
rtmb_mediation <- function(formula, data, family = "gaussian", prior = prior_uniform(), y_range = NULL, fixed = NULL, view = NULL, ...) {


  if (!is.list(formula)) stop("formula must be a list of formulas (e.g., list(M ~ X, Y ~ X + M)).")
  n_eq <- length(formula)

  # Validate: No random effects allowed in this version
  for (i in 1:n_eq) {
    if (any(grepl("\\|", as.character(formula[[i]])))) {
      stop("Multilevel mediation (random effects) is not yet supported in 'rtmb_mediation'. Please use fixed-effects formulas.")
    }
  }

  if (is.null(prior)) prior <- prior_uniform()

  # Automatically switch to prior_weak() if y_range is provided and prior is default uniform
  if (!is.null(y_range) && inherits(prior, "rtmb_prior") && prior$type == "uniform" &&
      is.null(prior$Intercept_sd) && is.null(prior$b_sd) && is.null(prior$sigma_rate)) {
    prior <- prior_weak()
  }

  # Prepare family list
  if (!is.list(family)) {
    family_list <- rep(list(family), n_eq)
  } else {
    if (length(family) != n_eq) stop("Length of family list must match length of formula list.")
    family_list <- family
  }

  N <- nrow(data)
  data_list <- list(N = N)
  resp_names <- character(n_eq)
  X_list <- list()
  X_colnames <- list()

  # 1. Parse Formulas and Prepare Data
  for (i in 1:n_eq) {
    f <- formula[[i]]
    mf <- model.frame(f, data = data)
    y_name <- as.character(f[[2]])
    resp_names[i] <- y_name

    Y_val <- as.numeric(model.response(mf))
    data_list[[paste0("Y_", i)]] <- Y_val
    X_mat <- model.matrix(f, data = data)
    data_list[[paste0("X_", i)]] <- X_mat
    cols <- colnames(X_mat)
    cols[cols == "(Intercept)"] <- "Intercept"
    X_colnames[[i]] <- cols
    X_list[[i]] <- X_mat

    # Add range variables for weak priors dynamically
    if (inherits(prior, "rtmb_prior") && prior$type %in% c("weak", "rhs", "ssp")) {
       f_type <- family_list[[i]]
       if (!(f_type %in% c("bernoulli", "binomial", "poisson"))) {
         range_i <- if (is.list(y_range)) y_range[[y_name]] else y_range
         if (is.null(range_i)) {
           stop(paste0("y_range is required for response variable '", y_name, "' when using weakly informative priors. ",
                       "Please provide y_range as a vector or a named list (e.g., y_range = list(", y_name, " = c(1, 5)))."))
         }
         data_list[[paste0("half_d_y_", i)]] <- diff(range_i) / 2
         data_list[[paste0("mid_y_", i)]] <- mean(range_i)
       }
    }
  }

  prior_type <- if (inherits(prior, "rtmb_prior")) prior$type else "uniform"
  if (prior_type == "weak") {
    if (is.null(prior$max_beta)) prior$max_beta <- 1.0
    if (is.null(prior$sd_ratio)) prior$sd_ratio <- 0.5
  }

  # 2. Setup AST Block
  setup_exprs <- list()
  setup_exprs[[1]] <- quote(N <- N)

  for (i in 1:n_eq) {
    f_type <- family_list[[i]]
    p_name <- paste0("b", i)
    X_name <- as.name(paste0("X_", i))

    has_b_prior <- prior_type == "weak" || !is.null(prior$b_sd) || !is.null(prior$Intercept_sd)
    has_sigma_prior <- f_type == "gaussian" && (prior_type == "weak" || !is.null(prior$sigma_rate))

    if (has_b_prior) {
      X_sd_name <- as.name(paste0("X_sd_", i))
      X_mean_name <- as.name(paste0("X_mean_", i))
      b_prior_sd_name <- as.name(paste0(p_name, "_prior_sd"))
      b_prior_mean_name <- as.name(paste0(p_name, "_prior_mean"))
      intercept_prior_sd_name <- as.name(paste0("intercept_prior_sd_", i))

      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_sd_name) <- apply(.(X_name), 2, sd))
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_sd_name)[.(X_sd_name) == 0] <- 1)
      setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_mean_name) <- rep(0, ncol(.(X_name))))

      if (prior_type == "weak") {
        X_c_name <- as.name(paste0("X_c_", i))
        mid_y_name <- as.name(paste0("mid_y_", i))

        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_mean_name) <- apply(.(X_name), 2, mean))
        has_intercept <- "Intercept" %in% X_colnames[[i]]
        if (has_intercept) {
           idx <- which(X_colnames[[i]] == "Intercept")
           setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_mean_name)[.(idx)] <- 0)
        }
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(X_c_name) <- .(X_name) - rep(1, N) %*% t(.(X_mean_name)))

        base_scale_name <- as.name(paste0("base_scale_", i))
        if (f_type %in% c("bernoulli", "binomial")) {
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(base_scale_name) <- pi / sqrt(3))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(intercept_prior_sd_name) <- .(base_scale_name) * .(prior$max_beta))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name) <- (.(prior$max_beta) * .(base_scale_name)) / .(X_sd_name))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(mid_y_name) <- 0)
        } else if (f_type == "poisson") {
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(base_scale_name) <- 1.0)
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(intercept_prior_sd_name) <- 1.0)
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name) <- 1.0 / .(X_sd_name))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(mid_y_name) <- 0)
        } else {
          half_d_y_name <- as.name(paste0("half_d_y_", i))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(base_scale_name) <- .(half_d_y_name) * .(prior$sd_ratio))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(intercept_prior_sd_name) <- .(half_d_y_name))
          setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name) <- (.(prior$max_beta) * .(base_scale_name)) / .(X_sd_name))
          if (has_sigma_prior) {
            sigma_rate_name <- as.name(paste0("sigma", i, "_rate"))
            setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(sigma_rate_name) <- 1.0 / .(base_scale_name))
          }
        }
        if (has_intercept) {
           idx <- which(X_colnames[[i]] == "Intercept")
           setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_mean_name)[.(idx)] <- .(mid_y_name))
        }
      } else {
        # Uniform with explicit SDs
        b_sd_val <- if (!is.null(prior$b_sd)) prior$b_sd else 10
        int_sd_val <- if (!is.null(prior$Intercept_sd)) prior$Intercept_sd else 10
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name) <- rep(.(b_sd_val), ncol(.(X_name))))
        setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(intercept_prior_sd_name) <- .(int_sd_val))
      }
      # Overwrite intercept SD if needed
      has_intercept <- "Intercept" %in% X_colnames[[i]]
      if (has_intercept) {
         idx <- which(X_colnames[[i]] == "Intercept")
         setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(b_prior_sd_name)[.(idx)] <- .(intercept_prior_sd_name))
      }
    }

    if (has_sigma_prior && prior_type != "weak") {
       sigma_rate_val <- if (!is.null(prior$sigma_rate)) prior$sigma_rate else 1
       sigma_rate_name <- as.name(paste0("sigma", i, "_rate"))
       setup_exprs[[length(setup_exprs) + 1]] <- bquote(.(sigma_rate_name) <- .(sigma_rate_val))
    }
  }

  setup_ast <- as.call(c(list(as.name("{")), setup_exprs))

  tmp_env <- list2env(data_list)
  eval(setup_ast, tmp_env)

  # 3. Parameters, Transform and Model Block AST
  param_exprs <- list()
  tran_exprs <- list()
  model_exprs <- list()
  generate_exprs <- list()

  v_names <- list()
  init_list <- list()
  b_vars <- c()
  s_vars <- c()

  for (i in 1:n_eq) {
    y_name <- resp_names[i]
    p_name <- paste0("b", i)
    p_c_name <- paste0("b_c", i)
    s_name <- paste0("sigma", i)
    f_type <- family_list[[i]]
    P_dim <- ncol(X_list[[i]])

    has_intercept <- "Intercept" %in% X_colnames[[i]]
    is_centered <- (prior_type == "weak" && has_intercept)
    target_p_name <- if (is_centered) p_c_name else p_name

    param_exprs[[length(param_exprs) + 1]] <- bquote(.(as.name(target_p_name)) <- Dim(.(P_dim)))

    if (is_centered) {
       v_names[[p_c_name]] <- X_colnames[[i]]
       v_names[[p_name]] <- X_colnames[[i]]
       init_list[[p_c_name]] <- rep(0, P_dim)

       idx <- which(X_colnames[[i]] == "Intercept")
       X_mean_name <- as.name(paste0("X_mean_", i))

       tran_exprs[[length(tran_exprs) + 1]] <- bquote(.(as.name(p_name)) <- .(as.name(p_c_name)))
       tran_exprs[[length(tran_exprs) + 1]] <- bquote(.(as.name(p_name))[.(idx)] <- .(as.name(p_c_name))[.(idx)] - sum(.(X_mean_name) * .(as.name(p_c_name))))

       lin_pred_expr <- bquote(.(as.name(paste0("X_c_", i))) %*% .(as.name(p_c_name)))
    } else {
       v_names[[p_name]] <- X_colnames[[i]]
       init_list[[p_name]] <- rep(0, P_dim)
       lin_pred_expr <- bquote(.(as.name(paste0("X_", i))) %*% .(as.name(p_name)))
    }

    if (f_type == "gaussian") {
      param_exprs[[length(param_exprs) + 1]] <- bquote(.(as.name(s_name)) <- Dim(lower = 0))
      init_list[[s_name]] <- 1.0
      model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(paste0("Y_", i))) ~ normal(.(lin_pred_expr), .(as.name(s_name))))
    } else if (f_type %in% c("binomial", "bernoulli")) {
      model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(paste0("Y_", i))) ~ bernoulli_logit(.(lin_pred_expr)))
    } else if (f_type == "poisson") {
      model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(paste0("Y_", i))) ~ poisson_log(.(lin_pred_expr)))
    }

    has_b_prior <- prior_type == "weak" || !is.null(prior$b_sd) || !is.null(prior$Intercept_sd)
    if (has_b_prior) {
      b_prior_sd_name <- as.name(paste0(p_name, "_prior_sd"))
      b_prior_mean_name <- as.name(paste0(p_name, "_prior_mean"))
      model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(target_p_name)) ~ normal(.(b_prior_mean_name), .(b_prior_sd_name)))
    }

    if (f_type == "gaussian") {
      has_sigma_prior <- prior_type == "weak" || !is.null(prior$sigma_rate)
      if (has_sigma_prior) {
        rate_name <- as.name(paste0(s_name, "_rate"))
        model_exprs[[length(model_exprs) + 1]] <- bquote(.(as.name(s_name)) ~ exponential(.(rate_name)))
      }
    }

    b_vars <- c(b_vars, p_name)
    if (f_type == "gaussian") s_vars <- c(s_vars, s_name)
  }

  # 4. Path Identification and Indirect Effects
  all_preds <- unique(unlist(X_colnames))
  all_resps <- unique(resp_names)
  mediators <- intersect(all_resps, all_preds)
  indeps <- setdiff(all_preds, c(all_resps, "Intercept"))
  if (length(indeps) == 0 && length(mediators) > 0) indeps <- setdiff(X_colnames[[1]], "Intercept")

  for (iv in indeps) {
    for (m in mediators) {
      idx_m_resp <- which(resp_names == m)
      idx_m_pred <- which(sapply(X_colnames, function(x) m %in% x))

      if (length(idx_m_resp) > 0 && length(idx_m_pred) > 0) {
        pos_iv <- which(X_colnames[[idx_m_resp]] == iv)
        if (length(pos_iv) > 0) {
          a_val <- bquote(.(as.name(paste0("b", idx_m_resp)))[.(pos_iv)])
          for (dv_idx in idx_m_pred) {
            if (dv_idx == idx_m_resp) next
            dv_name <- resp_names[dv_idx]
            pos_m <- which(X_colnames[[dv_idx]] == m)
            b_val <- bquote(.(as.name(paste0("b", dv_idx)))[.(pos_m)])
            ie_name <- paste0("IE_", iv, "_", m, "_", dv_name)
            generate_exprs[[ie_name]] <- bquote(.(as.name(ie_name)) <- .(a_val) * .(b_val))
            pos_iv_direct <- which(X_colnames[[dv_idx]] == iv)
            if (length(pos_iv_direct) > 0) {
              de_name <- paste0("DE_", iv, "_", dv_name)
              de_val <- bquote(.(as.name(paste0("b", dv_idx)))[.(pos_iv_direct)])
              generate_exprs[[de_name]] <- bquote(.(as.name(de_name)) <- .(de_val))
              te_name <- paste0("TE_", iv, "_", m, "_", dv_name)
              generate_exprs[[te_name]] <- bquote(.(as.name(te_name)) <- .(as.name(ie_name)) + .(de_val))
            }
          }
        }
      }
    }
  }

  # Calculate specific degrees of freedom for each equation
  df_map <- list()
  for (i in 1:n_eq) {
    df_val <- N - ncol(X_list[[i]])
    df_map[[paste0("b", i)]] <- df_val
    df_map[[paste0("sigma", i)]] <- df_val
  }
  # For IE/DE/TE, use the DF of the outcome equation
  for (iv in indeps) {
    for (m in mediators) {
      idx_m_resp <- which(resp_names == m)
      idx_m_pred <- which(sapply(X_colnames, function(x) m %in% x))
      if (length(idx_m_resp) > 0 && length(idx_m_pred) > 0) {
        for (dv_idx in idx_m_pred) {
          if (dv_idx == idx_m_resp) next
          dv_name <- resp_names[dv_idx]
          df_dv <- N - ncol(X_list[[dv_idx]])
          df_map[[paste0("IE_", iv, "_", m, "_", dv_name)]] <- df_dv
          df_map[[paste0("DE_", iv, "_", dv_name)]] <- df_dv
          df_map[[paste0("TE_", iv, "_", m, "_", dv_name)]] <- df_dv
        }
      }
    }
  }

  mdl_code <- list(
    setup = setup_ast,
    parameters = as.call(c(list(as.name("{")), param_exprs))
  )
  if (length(tran_exprs) > 0) {
    mdl_code$transform <- as.call(c(list(as.name("{")), tran_exprs))
  }
  mdl_code$model <- as.call(c(list(as.name("{")), model_exprs))
  mdl_code$generate <- as.call(c(list(as.name("{")), generate_exprs))
  mdl_code$env <- tmp_env

  view_order <- c(b_vars, names(generate_exprs), s_vars)

  # Ensure env is correctly formatted for rtmb_model
  # Using the evaluated objects directly from tmp_env
  mdl <- rtmb_model(data = as.list(tmp_env), code = mdl_code, par_names = v_names, init = init_list, fixed = fixed, view = view_order, silent = FALSE)
  mdl$formula <- formula
  mdl$raw_data <- data

  mdl$type <- "mediation"
  mdl$extra$df_pars <- if (prior_type == "weak") paste0("b_c", 1:n_eq) else paste0("b", 1:n_eq)

  return(mdl)
}

