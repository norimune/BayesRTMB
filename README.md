
# BayesRTMB

**Language / 言語**: [English
introduction](https://norimune.github.io/BayesRTMB/articles/introduction.html)
\| [パッケージ紹介](https://norimune.github.io/BayesRTMB/articles/ja-introduction.html)

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

You can install BayesRTMB from CRAN.

``` r
install.packages("BayesRTMB")
```

The development version can be installed from GitHub with either `pak` or
`remotes`.

``` r
pak::pak("norimune/BayesRTMB")
```

``` r
remotes::install_github("norimune/BayesRTMB")
```

### Windows Users

For ordinary use, Windows users can install the CRAN binary package
without Rtools. Rtools is only needed for source installation,
development, or compiling custom TMB C++ templates.

``` r
pkgbuild::check_build_tools(debug = TRUE)
```

If you install BayesRTMB from source and this check fails, install the
Rtools version that matches your R version from the [Rtools
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

You can also write a model directly.

``` r
Y <- debate$sat
X <- debate[c("talk","perf")] |> as.matrix()

data_list <- list(Y = Y, X = X)

code <- rtmb_code(
  setup = {
    N <- length(Y)
    K <- ncol(X)
  },
  parameters = {
    Intercept <- Dim(1)
    b <- Dim(K)
    sigma = Dim(lower = 0)
  },
  model = {
    mu <- Intercept + X %*% b
    Y ~ normal(mu, sigma)
    Intercept ~ normal(0, 10)
    b ~ normal(0, 10)
    sigma ~ exponential(1)
  }
)

mdl_custom <- rtmb_model(data_list, code)
fit_custom <- mdl_custom$sample()
```

## Articles

- [Introduction](https://norimune.github.io/BayesRTMB/articles/introduction.html):
  overall concepts, entry points, and inference workflow.
- [Quick Start](https://norimune.github.io/BayesRTMB/articles/quick_start.html):
  installation, a minimal model, MCMC diagnostics, visualization, and t
  tests.
- [Wrapper Functions](https://norimune.github.io/BayesRTMB/articles/wrapper_functions.html):
  using wrappers for Bayesian and classical analyses.
- [Hierarchical Models and
  GLMMs](https://norimune.github.io/BayesRTMB/articles/rtmb_glmer.html):
  detailed use of `rtmb_glmer()` for mixed models, GLMMs, priors,
  residual correlation, and visualization.
- [Writing Model
  Codes](https://norimune.github.io/BayesRTMB/articles/writing_models.html):
  writing custom models with `rtmb_code()`.
- [RTMB Internals and Inference
  Algorithms](https://norimune.github.io/BayesRTMB/articles/rtmb_internals.html):
  internal model representation, Laplace approximation, and inference
  pipelines.
- [日本語:
  パッケージ紹介](https://norimune.github.io/BayesRTMB/articles/ja-introduction.html)
- [日本語:
  クイックスタート](https://norimune.github.io/BayesRTMB/articles/ja-quick_start.html)
- [日本語:
  ラッパー関数の使い方](https://norimune.github.io/BayesRTMB/articles/ja-wrapper_functions.html)
- [日本語:
  モデルコードの書き方](https://norimune.github.io/BayesRTMB/articles/ja-writing_models.html)
