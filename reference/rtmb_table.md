# RTMB-based Contingency Table Analysis (Chi-squared Test)

\`rtmb_table\` performs a chi-squared test of independence between two
categorical variables. It provides both classic (frequentist) Pearson
chi-squared tests and Bayesian multinomial-style models.

## Usage

``` r
rtmb_table(
  x,
  y = NULL,
  data = NULL,
  correct = TRUE,
  prior = prior_uniform(),
  fixed = NULL,
  ...
)
```

## Arguments

- x:

  Variable name or formula.

- y:

  Variable name (optional if x is a formula).

- data:

  A data frame.

- correct:

  Logical; if TRUE, apply Yates' continuity correction (for 2x2 classic
  only).

- prior:

  Prior specification (Bayesian mode). Default is \`prior_uniform()\`.

- ...:

  Additional arguments.

- classic:

  Logical; if TRUE, perform frequentist chi-squared and Fisher's exact
  tests.

## Value

A \`Classic_Fit\` or \`MCMC_Fit\` object.

## Examples

``` r
# \donttest{
# Classic chi-squared test
rtmb_table(skill, cond, data = debate, classic = TRUE)
#> Pre-checking model code...
#> Checking RTMB setup...
#> <rtmb_table>
#>   Inherits from: <RTMB_Model>
#>   Public:
#>     build_ad_obj: function (init = NULL, laplace = FALSE, jacobian_target = "all", 
#>     calculate_reml_satterthwaite_df: function (ad_obj, opt_par, beta_idx) 
#>     calculate_satterthwaite_df: function (ad_obj, idx_fix_active = NULL, L_u_total = NULL, opt_par = NULL, 
#>     classic: function (df = "auto", df_pars = "none", REML = TRUE, view = NULL, 
#>     clone: function (deep = FALSE) 
#>     code: list
#>     constrained_vector_to_list: function (vec) 
#>     contrasts: NULL
#>     data: list
#>     extra: list
#>     family: NULL
#>     fixed_model: function (fixed_list) 
#>     formula: NULL
#>     generate: function (dat, para) 
#>     get_ad_obj: function (...) 
#>     get_n_obs: function () 
#>     get_par_list: function (init = NULL) 
#>     init: 0.852412934102303 0.118440864476572 0.0275952171199781 0 ...
#>     initialize: function (data, par_list, log_prob, transform = NULL, generate = NULL, 
#>     log_prob: function (dat, para) 
#>     map: NULL
#>     null_model: function (target, value = 0) 
#>     optimize: function (laplace = TRUE, init = NULL, num_estimate = 1, control = list(), 
#>     par_list: list
#>     par_names: list
#>     pl_full: list
#>     prepare_init: function (init_arg) 
#>     print_code: function () 
#>     print_generate: function () 
#>     print_log_prob: function () 
#>     print_transform: function () 
#>     prior_correction: 0
#>     raw_data: NULL
#>     requested_contrasts: NULL
#>     sample: function (sampling = 1000, warmup = 1000, chains = 4, thin = 1, 
#>     to_constrained: function (par_list_unconstrained) 
#>     to_unconstrained: function (par_list_constrained) 
#>     transform: function (dat, para) 
#>     type: table
#>     unconstrained_vector_to_list: function (vec) 
#>     variational: function (iter = 3000, tol_rel_obj = 0.005, window_size = 100, 
#>     view: NULL
# }
```
