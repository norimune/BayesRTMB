# Define an RTMB Model with Stan-like Syntax

`rtmb_code` is the primary interface for defining Bayesian or
Frequentist models within the RTMB framework. It allows you to organize
model logic into functional blocks while utilizing a Stan-like `~`
operator for specifying priors and likelihoods.

## Usage

``` r
rtmb_code(...)
```

## Arguments

- ...:

  Model definition blocks.

## Value

A list of unevaluated code blocks.

## Details

**Key Blocks and Logic:**

- `setup`: Pre-compilation data processing.

- `parameters`: Declaration of estimated parameters. See
  [`parameter_types`](https://norimune.github.io/BayesRTMB/reference/parameter_types.md)
  for available constraints.

- `transform`: Definition of intermediate variables. Use
  [`math_functions`](https://norimune.github.io/BayesRTMB/reference/math_functions.md)
  for stable AD calculations.

- `model`: Likelihood and Priors. Refer to
  [`distributions`](https://norimune.github.io/BayesRTMB/reference/distributions.md)
  for available sampling statements.

- `generate`: Post-hoc calculation of predictions or diagnostics.

**Automatic Differentiation (AD) Note:** All code within `parameters`,
`transform`, and `model` is automatically differentiated by RTMB. To
ensure numerical stability, use provided utility functions like
[`log_sum_exp`](https://norimune.github.io/BayesRTMB/reference/log_sum_exp.md)
or
[`inv_logit`](https://norimune.github.io/BayesRTMB/reference/inv_logit.md)
instead of raw algebraic implementations.

## See also

[`rtmb_model`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md)
for building the model object.
[`parameter_types`](https://norimune.github.io/BayesRTMB/reference/parameter_types.md)
for constraining parameters.
[`distributions`](https://norimune.github.io/BayesRTMB/reference/distributions.md)
for likelihood functions.
[`math_functions`](https://norimune.github.io/BayesRTMB/reference/math_functions.md)
for numerical utilities.

## Examples

``` r
# \donttest{
# Simulate data for a linear regression
set.seed(123)
N <- 100
x <- rnorm(N)
y <- 2.0 + 1.5 * x + rnorm(N, mean = 0, sd = 1)
data_list <- list(N = N, x = x, y = y)

# Define the model using rtmb_code
code <- rtmb_code(
  setup = {
    # Center the predictor variable (executed once)
    x_centered <- x - mean(x)
  },
  parameters = {
    # Define parameters and their constraints
    alpha = Dim(1)
    beta  = Dim(1)
    sigma = Dim(1, lower = 0)
  },
  transform = {
    # Calculate the linear predictor
    mu <- alpha + beta * x_centered
  },
  model = {
    # Priors
    alpha ~ normal(0, 10)
    beta  ~ normal(0, 10)
    sigma ~ exponential(1)

    # Likelihood (Vectorized)
    y ~ normal(mu, sigma)
  },
  generate = {
    # Calculate generated quantities
    y_pred <- mu

    # Must return a named list
    list(y_pred = y_pred)
  }
)

# Create the model object
mod <- rtmb_model(data = data_list, code = code)
#> Pre-checking model code...
#> Checking RTMB setup...

# Fit the model using MAP estimation
map_res <- mod$optimize()
#> Starting optimization...
#> 
#> Optimization converged. Final objective: 145.34
map_res$summary(pars = c("alpha", "beta", "sigma"))
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 145.34
#> Approx. Log Marginal Likelihood (Laplace): -149.89
#> 
#> Point Estimates and 95% Wald CI:
#> variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> alpha      2.02788     0.09564    1.84043    2.21532 
#> beta       1.44737     0.10530    1.24099    1.65374 
#> sigma      0.95544     0.06708    0.83261    1.09639 
#> 

# The generated quantity 'y_pred' can also be summarized
map_res$summary("y_pred", max_rows = 5)
#> 
#> Call:
#> MAP Estimation via RTMB
#> 
#> Negative Log-Posterior: 145.34
#> Approx. Log Marginal Likelihood (Laplace): -149.89
#> 
#> Point Estimates and 95% Wald CI:
#>  variable  Estimate  Std. Error  Lower 95%  Upper 95% 
#> y_pred[1]   1.08581          NA         NA         NA 
#> y_pred[2]   1.56387          NA         NA         NA 
#> y_pred[3]   4.15305          NA         NA         NA 
#> y_pred[4]   1.99908          NA         NA         NA 
#> y_pred[5]   2.08415          NA         NA         NA 
#> 
# }
```
