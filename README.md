# BayesRTMB

**Language / 言語**: [English
introduction](https://norimune.github.io/BayesRTMB/articles/introduction.html)
| [パッケージ紹介](https://norimune.github.io/BayesRTMB/articles/ja-introduction.html)

**BayesRTMB** is an R package for writing and fitting statistical models
with RTMB as the automatic differentiation engine.

You can start from wrapper functions such as `rtmb_lm()`,
`rtmb_glmer()`, `rtmb_corr()`, and `rtmb_ttest()`, or write your own
model with `rtmb_code()`. The same model object can then be used for
MCMC, MAP estimation, variational inference, and frequency-oriented
classical analyses where supported.

## Key Features

- **Model code in R**: Write models with `setup`, `parameters`,
  `transform`, `model`, and `generate` blocks.
- **Multiple estimation methods**: Use `$sample()`, `$optimize()`,
  `$variational()`, and `$classic()` from a common model object.
- **Wrapper functions**: Fit regression models, GLM/GLMMs, t tests,
  correlations, contingency tables, factor analysis, IRT, latent rank
  models, and multidimensional unfolding models.
- **Classical analyses from wrappers**: Use `$classic()` for
  frequency-oriented outputs, including `AIC()`, `BIC()`, and `anova()`
  for supported fits.
- **Random effects and Laplace approximation**: Declare parameters with
  `random = TRUE` and use RTMB’s Laplace machinery for latent variables
  and mixed models.
- **Diagnostics and visualization**: Plot posterior draws, convergence
  diagnostics, forest plots, conditional effects, and related summaries.

## Installation

You can install the development version from GitHub with either `pak` or
`remotes`.

``` r
pak::pak("norimune/BayesRTMB")
```

``` r
remotes::install_github("norimune/BayesRTMB")
```

### Windows Users

Windows users need Rtools because BayesRTMB depends on the RTMB/TMB
ecosystem.

``` r
pkgbuild::check_build_tools(debug = TRUE)
```

If this check fails, install the Rtools version that matches your R
version from the [Rtools
page](https://cran.r-project.org/bin/windows/Rtools/), restart R, and
try again.

## Quick Example

For standard analyses, start with a wrapper function.

``` r
library(BayesRTMB)
data(debate)

mdl <- rtmb_lm(sat ~ talk * perf, data = debate)

fit_mcmc <- mdl$sample()
fit_map  <- mdl$optimize()
fit_lm   <- mdl$classic()
```

The generated model code can be inspected.

``` r
mdl$print_code()
```

You can also write a model directly.

``` r
data_list <- list(Y = debate$sat, N = nrow(debate))

code <- rtmb_code(
  parameters = {
    mu = Dim()
    sigma = Dim(lower = 0)
  },
  model = {
    Y ~ normal(mu, sigma)
    mu ~ normal(0, 10)
    sigma ~ exponential(1)
  }
)

mdl_custom <- rtmb_model(data_list, code)
fit_custom <- mdl_custom$sample()
```

## Articles

- [Quick Start](https://norimune.github.io/BayesRTMB/articles/quick_start.html):
  installation, a minimal model, MCMC diagnostics, visualization, and t
  tests.
- [Wrapper Functions](https://norimune.github.io/BayesRTMB/articles/wrapper_functions.html):
  using wrappers for Bayesian and classical analyses.
- [Writing Model Codes](https://norimune.github.io/BayesRTMB/articles/writing_models.html):
  writing custom models with `rtmb_code()`.
- [日本語: パッケージ紹介](https://norimune.github.io/BayesRTMB/articles/ja-introduction.html)
- [日本語: クイックスタート](https://norimune.github.io/BayesRTMB/articles/ja-quick_start.html)
- [日本語: ラッパー関数の使い方](https://norimune.github.io/BayesRTMB/articles/ja-wrapper_functions.html)
- [日本語: モデルコードの書き方](https://norimune.github.io/BayesRTMB/articles/ja-writing_models.html)
