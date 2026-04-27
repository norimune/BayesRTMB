# Package index

## Core Modeling

Functions for defining and fitting custom models.

- [`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
  : Define an RTMB Model with Stan-like Syntax
- [`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md)
  : Create an RTMB_Model Object
- [`safe_rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/safe_rtmb_model.md)
  : Safe RTMB model construction (with error message translation)
- [`Dim()`](https://norimune.github.io/BayesRTMB/reference/Dim.md) :
  Define parameter dimensions and types

## Wrapper Functions

Functions for quick, standard analyses.

- [`rtmb_lm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lm.md)
  : RTMB-based Linear Regression wrapper function
- [`rtmb_glm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glm.md)
  : RTMB-based GLM wrapper function (no random effects)
- [`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md)
  : RTMB-based GLMM wrapper function
- [`rtmb_ttest()`](https://norimune.github.io/BayesRTMB/reference/rtmb_ttest.md)
  : RTMB-based Bayesian two-sample t-test wrapper function
- [`rtmb_corr()`](https://norimune.github.io/BayesRTMB/reference/rtmb_corr.md)
  : Wrapper for estimating correlation matrix (multivariate normal
  distribution)
- [`rtmb_fa()`](https://norimune.github.io/BayesRTMB/reference/rtmb_fa.md)
  : RTMB-based Factor Analysis Wrapper
- [`rtmb_irt()`](https://norimune.github.io/BayesRTMB/reference/rtmb_irt.md)
  : RTMB-based IRT (Item Response Theory) Wrapper

## Plotting and Diagnostics

Functions for visualizing MCMC samples and model results.

- [`plot_acf()`](https://norimune.github.io/BayesRTMB/reference/plot_acf.md)
  : Plot autocorrelation for one variable across chains
- [`plot_dens()`](https://norimune.github.io/BayesRTMB/reference/plot_dens.md)
  : Plot posterior densities for MCMC samples
- [`plot_forest()`](https://norimune.github.io/BayesRTMB/reference/plot_forest.md)
  : Plot parameter estimates and credible intervals (Forest Plot)
- [`plot_pairs()`](https://norimune.github.io/BayesRTMB/reference/plot_pairs.md)
  : Plot pairs for posterior samples
- [`plot_trace()`](https://norimune.github.io/BayesRTMB/reference/plot_trace.md)
  : Plot MCMC trace plots

## Post-Estimation and Utilities

Tools for model evaluation and post-processing.

- [`bayes_factor()`](https://norimune.github.io/BayesRTMB/reference/bayes_factor.md)
  : Calculate Bayes factor from log marginal likelihoods
- [`conditional_effects()`](https://norimune.github.io/BayesRTMB/reference/conditional_effects.md)
  : Calculate Conditional Effects
- [`item_curve()`](https://norimune.github.io/BayesRTMB/reference/item_curve.md)
  : Calculate Item Response Curve / Category Response Curve
- [`item_info()`](https://norimune.github.io/BayesRTMB/reference/item_info.md)
  : Calculate Item Information Function
- [`test_info()`](https://norimune.github.io/BayesRTMB/reference/test_info.md)
  : Calculate Test Information Function
- [`sort_loadings()`](https://norimune.github.io/BayesRTMB/reference/sort_loadings.md)
  : Sort and display factor loadings neatly
- [`read_mcmc_csv()`](https://norimune.github.io/BayesRTMB/reference/read_mcmc_csv.md)
  : Restore MCMC Fit from CSV

## Math and Transformation Functions

Functions used within rtmb_code for stability and transformations.

- [`logit()`](https://norimune.github.io/BayesRTMB/reference/logit.md) :
  Logit function
- [`inv_logit()`](https://norimune.github.io/BayesRTMB/reference/inv_logit.md)
  : Inverse logit function
- [`log1m()`](https://norimune.github.io/BayesRTMB/reference/log1m.md) :
  Log of one minus x
- [`log1m_exp()`](https://norimune.github.io/BayesRTMB/reference/log1m_exp.md)
  : Log of one minus exponential of x
- [`log1p_exp()`](https://norimune.github.io/BayesRTMB/reference/log1p_exp.md)
  : Log of one plus exponential of x
- [`log_sum_exp()`](https://norimune.github.io/BayesRTMB/reference/log_sum_exp.md)
  : Log-sum-exp function
- [`log_softmax()`](https://norimune.github.io/BayesRTMB/reference/log_softmax.md)
  : Log-softmax function
- [`softmax()`](https://norimune.github.io/BayesRTMB/reference/softmax.md)
  : Softmax function
- [`log_mix()`](https://norimune.github.io/BayesRTMB/reference/log_mix.md)
  : Log mixture of two probabilities
- [`log_det_chol()`](https://norimune.github.io/BayesRTMB/reference/log_det_chol.md)
  : Log determinant of a Cholesky factor
- [`quad_form_chol()`](https://norimune.github.io/BayesRTMB/reference/quad_form_chol.md)
  : Quadratic form using a Cholesky factor
- [`quad_form_diag()`](https://norimune.github.io/BayesRTMB/reference/quad_form_diag.md)
  : Quadratic form with a diagonal matrix
- [`distance()`](https://norimune.github.io/BayesRTMB/reference/distance.md)
  : Euclidean distance
- [`squared_distance()`](https://norimune.github.io/BayesRTMB/reference/squared_distance.md)
  : Squared Euclidean distance
- [`fabs()`](https://norimune.github.io/BayesRTMB/reference/fabs.md) :
  Smooth absolute value function
- [`to_centered_matrix()`](https://norimune.github.io/BayesRTMB/reference/to_centered_matrix.md)
  : Vector to centered matrix (RTMB compatible)
- [`to_centered_tri()`](https://norimune.github.io/BayesRTMB/reference/to_centered_tri.md)
  : Vector to centered triangular matrix (RTMB compatible)
- [`sum_to_zero()`](https://norimune.github.io/BayesRTMB/reference/sum_to_zero.md)
  : Sum-to-zero transformation
- [`stz_basis()`](https://norimune.github.io/BayesRTMB/reference/stz_basis.md)
  : stz basis function

## Classes

- [`RTMB_Fit_Base`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.md)
  : Base class for RTMB Fit objects
- [`MAP_Fit`](https://norimune.github.io/BayesRTMB/reference/MAP_Fit.md)
  : MAP fit object
- [`MCMC_Fit`](https://norimune.github.io/BayesRTMB/reference/MCMC_Fit.md)
  : MCMC fit object
- [`VB_Fit`](https://norimune.github.io/BayesRTMB/reference/VB_Fit.md) :
  VB fit object

## Datasets

- [`beverage`](https://norimune.github.io/BayesRTMB/reference/beverage.md)
  : Beverage Preference Data
- [`BigFive`](https://norimune.github.io/BayesRTMB/reference/BigFive.md)
  : Big Five Personality Traits Data
- [`discussion`](https://norimune.github.io/BayesRTMB/reference/discussion.md)
  : Group Discussion Simulation Data

## Documentation Topics

- [`distributions`](https://norimune.github.io/BayesRTMB/reference/distributions.md)
  : Probability Distributions for RTMB Models
- [`math_functions`](https://norimune.github.io/BayesRTMB/reference/math_functions.md)
  : Mathematical and Matrix Utility Functions for RTMB Models
- [`model_code()`](https://norimune.github.io/BayesRTMB/reference/model_code.md)
  : Model Code Wrapper for RTMB
- [`parameters_code()`](https://norimune.github.io/BayesRTMB/reference/parameters_code.md)
  : Code block for parameter definitions
- [`parameter_types`](https://norimune.github.io/BayesRTMB/reference/parameter_types.md)
  : Parameter Types and Constraints in RTMB Models

## Other Objects and Internals

- [`ADVI_method()`](https://norimune.github.io/BayesRTMB/reference/ADVI_method.md)
  : Automatic Differentiation Variational Inference (ADVI)
- [`BigFive`](https://norimune.github.io/BayesRTMB/reference/BigFive.md)
  : Big Five Personality Traits Data
- [`Dim()`](https://norimune.github.io/BayesRTMB/reference/Dim.md) :
  Define parameter dimensions and types
- [`MAP_Fit`](https://norimune.github.io/BayesRTMB/reference/MAP_Fit.md)
  : MAP fit object
- [`MCMC_Fit`](https://norimune.github.io/BayesRTMB/reference/MCMC_Fit.md)
  : MCMC fit object
- [`RTMB_Fit_Base`](https://norimune.github.io/BayesRTMB/reference/RTMB_Fit_Base.md)
  : Base class for RTMB Fit objects
- [`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md)
  : Create an RTMB_Model Object
- [`VB_Fit`](https://norimune.github.io/BayesRTMB/reference/VB_Fit.md) :
  VB fit object
- [`bayes_factor()`](https://norimune.github.io/BayesRTMB/reference/bayes_factor.md)
  : Calculate Bayes factor from log marginal likelihoods
- [`beverage`](https://norimune.github.io/BayesRTMB/reference/beverage.md)
  : Beverage Preference Data
- [`conditional_effects()`](https://norimune.github.io/BayesRTMB/reference/conditional_effects.md)
  : Calculate Conditional Effects
- [`conditional_effects(`*`<mcmc_fit>`*`)`](https://norimune.github.io/BayesRTMB/reference/conditional_effects.mcmc_fit.md)
  : Calculate conditional effects for MCMC fit objects
- [`discussion`](https://norimune.github.io/BayesRTMB/reference/discussion.md)
  : Group Discussion Simulation Data
- [`distance()`](https://norimune.github.io/BayesRTMB/reference/distance.md)
  : Euclidean distance
- [`distributions`](https://norimune.github.io/BayesRTMB/reference/distributions.md)
  : Probability Distributions for RTMB Models
- [`fabs()`](https://norimune.github.io/BayesRTMB/reference/fabs.md) :
  Smooth absolute value function
- [`generate_random_init()`](https://norimune.github.io/BayesRTMB/reference/generate_random_init.md)
  : Generate Random Initial Values
- [`inv_logit()`](https://norimune.github.io/BayesRTMB/reference/inv_logit.md)
  : Inverse logit function
- [`item_curve(`*`<RTMB_Fit_Base>`*`)`](https://norimune.github.io/BayesRTMB/reference/item_curve.RTMB_Fit_Base.md)
  : Item Response Curve for RTMB_Fit_Base
- [`item_curve()`](https://norimune.github.io/BayesRTMB/reference/item_curve.md)
  : Calculate Item Response Curve / Category Response Curve
- [`item_info(`*`<RTMB_Fit_Base>`*`)`](https://norimune.github.io/BayesRTMB/reference/item_info.RTMB_Fit_Base.md)
  : Item Information Function for RTMB_Fit_Base
- [`item_info()`](https://norimune.github.io/BayesRTMB/reference/item_info.md)
  : Calculate Item Information Function
- [`log1m()`](https://norimune.github.io/BayesRTMB/reference/log1m.md) :
  Log of one minus x
- [`log1m_exp()`](https://norimune.github.io/BayesRTMB/reference/log1m_exp.md)
  : Log of one minus exponential of x
- [`log1p_exp()`](https://norimune.github.io/BayesRTMB/reference/log1p_exp.md)
  : Log of one plus exponential of x
- [`log_det_chol()`](https://norimune.github.io/BayesRTMB/reference/log_det_chol.md)
  : Log determinant of a Cholesky factor
- [`log_mix()`](https://norimune.github.io/BayesRTMB/reference/log_mix.md)
  : Log mixture of two probabilities
- [`log_softmax()`](https://norimune.github.io/BayesRTMB/reference/log_softmax.md)
  : Log-softmax function
- [`log_sum_exp()`](https://norimune.github.io/BayesRTMB/reference/log_sum_exp.md)
  : Log-sum-exp function
- [`logit()`](https://norimune.github.io/BayesRTMB/reference/logit.md) :
  Logit function
- [`math_functions`](https://norimune.github.io/BayesRTMB/reference/math_functions.md)
  : Mathematical and Matrix Utility Functions for RTMB Models
- [`model_code()`](https://norimune.github.io/BayesRTMB/reference/model_code.md)
  : Model Code Wrapper for RTMB
- [`parameter_types`](https://norimune.github.io/BayesRTMB/reference/parameter_types.md)
  : Parameter Types and Constraints in RTMB Models
- [`parameters_code()`](https://norimune.github.io/BayesRTMB/reference/parameters_code.md)
  : Code block for parameter definitions
- [`plot(`*`<ce_rtmb>`*`)`](https://norimune.github.io/BayesRTMB/reference/plot.ce_rtmb.md)
  : Plot method for ce_rtmb class (Base R)
- [`plot_acf()`](https://norimune.github.io/BayesRTMB/reference/plot_acf.md)
  : Plot autocorrelation for one variable across chains
- [`plot_dens()`](https://norimune.github.io/BayesRTMB/reference/plot_dens.md)
  : Plot posterior densities for MCMC samples
- [`plot_forest()`](https://norimune.github.io/BayesRTMB/reference/plot_forest.md)
  : Plot parameter estimates and credible intervals (Forest Plot)
- [`plot_pairs()`](https://norimune.github.io/BayesRTMB/reference/plot_pairs.md)
  : Plot pairs for posterior samples
- [`plot_trace()`](https://norimune.github.io/BayesRTMB/reference/plot_trace.md)
  : Plot MCMC trace plots
- [`print(`*`<bayes_factor>`*`)`](https://norimune.github.io/BayesRTMB/reference/print.bayes_factor.md)
  : Print method for bayes_factor objects
- [`print(`*`<ce_rtmb>`*`)`](https://norimune.github.io/BayesRTMB/reference/print.ce_rtmb.md)
  : Print method for ce_rtmb class (automatically calls plot)
- [`print(`*`<summary_BayesRTMB>`*`)`](https://norimune.github.io/BayesRTMB/reference/print.summary_BayesRTMB.md)
  : print for summary_BayesRTMB class
- [`quad_form_chol()`](https://norimune.github.io/BayesRTMB/reference/quad_form_chol.md)
  : Quadratic form using a Cholesky factor
- [`quad_form_diag()`](https://norimune.github.io/BayesRTMB/reference/quad_form_diag.md)
  : Quadratic form with a diagonal matrix
- [`read_mcmc_csv()`](https://norimune.github.io/BayesRTMB/reference/read_mcmc_csv.md)
  : Restore MCMC Fit from CSV
- [`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
  : Define an RTMB Model with Stan-like Syntax
- [`rtmb_corr()`](https://norimune.github.io/BayesRTMB/reference/rtmb_corr.md)
  : Wrapper for estimating correlation matrix (multivariate normal
  distribution)
- [`rtmb_fa()`](https://norimune.github.io/BayesRTMB/reference/rtmb_fa.md)
  : RTMB-based Factor Analysis Wrapper
- [`rtmb_glm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glm.md)
  : RTMB-based GLM wrapper function (no random effects)
- [`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md)
  : RTMB-based GLMM wrapper function
- [`rtmb_irt()`](https://norimune.github.io/BayesRTMB/reference/rtmb_irt.md)
  : RTMB-based IRT (Item Response Theory) Wrapper
- [`rtmb_lm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lm.md)
  : RTMB-based Linear Regression wrapper function
- [`rtmb_syntax`](https://norimune.github.io/BayesRTMB/reference/rtmb_syntax.md)
  : Guidelines for Writing RTMB-Compatible Code
- [`rtmb_ttest()`](https://norimune.github.io/BayesRTMB/reference/rtmb_ttest.md)
  : RTMB-based Bayesian two-sample t-test wrapper function
- [`rtmb_wrappers`](https://norimune.github.io/BayesRTMB/reference/rtmb_wrappers.md)
  : Common Features and Arguments of RTMB Wrapper Functions
- [`safe_rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/safe_rtmb_model.md)
  : Safe RTMB model construction (with error message translation)
- [`softmax()`](https://norimune.github.io/BayesRTMB/reference/softmax.md)
  : Softmax function
- [`sort_loadings()`](https://norimune.github.io/BayesRTMB/reference/sort_loadings.md)
  : Sort and display factor loadings neatly
- [`squared_distance()`](https://norimune.github.io/BayesRTMB/reference/squared_distance.md)
  : Squared Euclidean distance
- [`stz_basis()`](https://norimune.github.io/BayesRTMB/reference/stz_basis.md)
  : stz basis function
- [`sum_to_zero()`](https://norimune.github.io/BayesRTMB/reference/sum_to_zero.md)
  : Sum-to-zero transformation
- [`test_info()`](https://norimune.github.io/BayesRTMB/reference/test_info.md)
  : Calculate Test Information Function
- [`to_centered_matrix()`](https://norimune.github.io/BayesRTMB/reference/to_centered_matrix.md)
  : Vector to centered matrix (RTMB compatible)
- [`to_centered_tri()`](https://norimune.github.io/BayesRTMB/reference/to_centered_tri.md)
  : Vector to centered triangular matrix (RTMB compatible)
- [`to_lower_tri()`](https://norimune.github.io/BayesRTMB/reference/to_lower_tri.md)
  : Vector to lower triangular matrix (RTMB compatible)
- [`transform_code()`](https://norimune.github.io/BayesRTMB/reference/transform_code.md)
  : Transformed Code Wrapper for RTMB
- [`validate_data()`](https://norimune.github.io/BayesRTMB/reference/validate_data.md)
  : Pre-validation of data and parameters
