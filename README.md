# BayesRTMB

**Language / 言語**: [English introduction](https://norimune.github.io/BayesRTMB/articles/introduction.html)
| [パッケージ紹介](https://norimune.github.io/BayesRTMB/articles/ja-introduction.html)

**BayesRTMB** is an R package for fast, flexible, and compile-free
Bayesian modeling.

Powered by **RTMB** (R template Model Builder) as its automatic
differentiation engine, BayesRTMB allows you to write models natively in
R with a syntax similar to Stan, while eliminating the need for lengthy
C++ compilation times.

## Key Features

- **No C++ Compilation**: Write models in pure R and execute them
  instantly.
- **Intuitive Block Syntax**: Define your model cleanly using
  `rtmb_code()` with `parameters`, `transform`, `model`, and `generate`
  blocks.
- **Multiple Inference Engines**: Seamlessly switch between:
  - **NUTS** (No-U-Turn Sampler) for robust MCMC sampling.
  - **ADVI** (Automatic Differentiation Variational Inference) for fast
    approximate inference.
  - **MAP** (Maximum a Posteriori) estimation.
- **Laplace Approximation**: Built-in support for marginalizing random
  effects to speed up complex hierarchical models.
- **Model Comparison**: Compute marginal likelihoods via Bridge Sampling
  and calculate Bayes Factors.
- **Rich Diagnostics & Plotting**: Integrated functions for trace plots,
  density plots, autocorrelation, and forest plots.
- **Convenient Wrappers**: Ready-to-use functions like `rtmb_ttest` and
  `rtmb_glmer` for standard analyses with automatic weakly informative
  priors.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("norimune/BayesRTMB")
```

### Important Note for Windows Users

Since BayesRTMB depends on C++ libraries through the `RTMB` and `TMB`
ecosystems, Windows users **must install Rtools** to compile the
underlying dependencies.

1.  Download and install the appropriate version of
    [Rtools](https://cran.r-project.org/bin/windows/Rtools/) that
    matches your R version (e.g., Rtools44 for R 4.4.x).
2.  After installation, restart RStudio.
3.  You can verify if Rtools is correctly installed and recognized by
    your system by running:

``` r
# If TRUE is returned, your system is ready to build packages.
pkgbuild::check_build_tools(debug = TRUE)
```

## Basic Usage

A model in BayesRTMB is defined by combining your data with an
`rtmb_code()` block. Here is a simple example of estimating a mean and
standard deviation.

``` r
library(BayesRTMB)

# 1. Prepare data
Y <- c(5.2, 4.8, 5.5, 6.1, 4.9, 5.3)
dat <- list(Y = Y, N = length(Y))

# 2. Define the model using rtmb_code
code <- rtmb_code(
  parameters = {
    mu = Dim(1)
    sigma = Dim(1, lower = 0)
  },
  model = {
    # Priors
    mu ~ normal(0, 10)
    sigma ~ exponential(0.1)
    
    # Likelihood
    Y ~ normal(mu, sigma)
  }
)

# 3. Build the model object
mdl <- rtmb_model(data = dat, code = code)

# 4. Run NUTS MCMC sampling
fit_mcmc <- mdl$sample(iter = 2000, warmup = 1000, chains = 4)

# Print summary
fit_mcmc$summary()
```

## Advanced Inference

BayesRTMB provides a unified interface for various inference methods
using R6 objects. You can easily apply different algorithms to the exact
same model object.

``` r
# Maximum a Posteriori (MAP) estimation
fit_map <- mdl$optimize()

# Automatic Differentiation Variational Inference (ADVI)
# Supports "meanfield", "fullrank", and "hybrid" approaches
fit_advi <- mdl$advi(method = "meanfield", iter = 5000)
```

### Dealing with Random Effects

For hierarchical models, you can declare parameters as random effects.
BayesRTMB can utilize RTMB’s **Laplace approximation** to integrate them
out, significantly speeding up both MAP and NUTS.

``` r
code_hierarchical <- rtmb_code(
  parameters = {
    mu_global = Dim(1)
    sigma_global = Dim(1, lower = 0)
    # Declare a random effect vector
    alpha = Dim(J, random = TRUE) 
  },
  model = {
    alpha ~ normal(mu_global, sigma_global)
    # ... rest of the model ...
  }
)

# Enable Laplace approximation during sampling
fit_mcmc_laplace <- mdl$sample(laplace = TRUE)
```

## Model Comparison & Bayes Factors

BayesRTMB integrates bridge sampling to easily compute marginal
likelihoods for model comparison. You can evaluate the strength of
evidence for one model over another using Bayes factors.

``` r
# Assume fit_model1 and fit_model2 are MCMC_Fit objects
# Compute marginal likelihoods via Bridge Sampling
bs1 <- fit_model1$bridgesampling()
bs2 <- fit_model2$bridgesampling()

# Calculate Bayes Factor
# This function outputs the BF and an interpretation based on Kass & Raftery (1995)
bayes_factor(bs1, bs2)
```

## Visualization

The package includes built-in functions to visualize posterior
distributions and diagnose MCMC performance quickly.

``` r
# Extract the underlying array of samples
samples <- fit_mcmc$fit

# Diagnostic plots
plot_trace(samples)
plot_dens(samples)
plot_acf(samples, var_idx = 1)

# Compare parameter estimates
plot_forest(samples, prob = 0.95)
plot_pairs(samples)
```

## Built-in Wrappers

If you don’t want to write model code from scratch, BayesRTMB offers
wrapper functions for common tasks, complete with automated weakly
informative priors.

``` r
# Bayesian Two-Sample T-Test
# Automatically scales priors based on y_range
mdl_test <- rtmb_ttest(
  y1 = group_A, 
  y2 = group_B, 
  y_range = c(0,10)
)

fit_test <- mdl_test$sample()
```
