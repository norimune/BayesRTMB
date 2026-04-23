# BayesRTMB 日本語紹介

## BayesRTMB 日本語紹介

**BayesRTMB** は、RTMB を自動微分エンジンとして用いる R パッケージです。
Stan に近い感覚でモデルを書きつつ、R
の中でそのままベイズ推定を進められます。

このページでは、BayesRTMB
の位置づけと、最初に押さえておくとよい使い方を日本語でまとめます。

### BayesRTMB でできること

BayesRTMB には、次のような特徴があります。

- **C++ のコンパイルなしでモデルを書ける**
- **[`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
  によるブロック構文**でモデルを整理して書ける
- **NUTS / ADVI / MAP** を同じモデルオブジェクトから使い分けられる
- **Laplace 近似**により、階層モデルの random effect を効率よく扱える
- **Bridge Sampling** を使った周辺尤度や Bayes factor
  の計算に対応している
- **[`rtmb_lm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lm.md)
  や
  [`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md)
  などのラッパー関数**で、標準的な分析をすぐに始められる

### インストール

GitHub 版をインストールする場合は、次のようにします。

``` r
# install.packages("remotes")
remotes::install_github("norimune/BayesRTMB")
```

``` r
library(BayesRTMB)
```

### 基本の流れ

BayesRTMB の基本的な流れは次の 3 段階です。

1.  **データ**を用意する
2.  **[`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)**
    でモデルを書く
3.  **[`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md)**
    でモデルオブジェクトを作り、推定する

最小例として、平均 `mu` と標準偏差 `sigma`
を推定する単純な正規モデルを示します。

``` r
library(BayesRTMB)

# 1. データ
Y <- c(5.2, 4.8, 5.5, 6.1, 4.9, 5.3)
dat <- list(Y = Y, N = length(Y))

# 2. モデル
code <- rtmb_code(
  parameters = {
    mu    = Dim(1)
    sigma = Dim(1, lower = 0)
  },
  model = {
    mu    ~ normal(0, 10)
    sigma ~ exponential(0.1)
    Y     ~ normal(mu, sigma)
  }
)

# 3. モデルオブジェクトの作成
mdl <- rtmb_model(data = dat, code = code)
```

### 推定法の使い分け

同じ `mdl` に対して、複数の推定法を使えます。

#### MAP 推定

まずモデルが正しく動くかを素早く確認したいときに便利です。

``` r
fit_map <- mdl$optimize()
fit_map$summary()
```

#### NUTS による MCMC

最も標準的なベイズ推定です。事後分布の要約、区間推定、収束診断を行いたいときに使います。

``` r
fit_mcmc <- mdl$sample(sampling = 1000, warmup = 1000, chains = 4)
fit_mcmc$summary()
```

#### ADVI による変分推論

近似推定でよいので高速に結果を得たいときに向いています。

``` r
fit_vb <- mdl$variational(
  method = "meanfield",
  iter = 3000,
  num_estimate = 4
)
fit_vb$summary()
```

### random effect を含むモデル

階層モデルでは、パラメータを `random = TRUE` として宣言できます。

``` r
code_hier <- rtmb_code(
  parameters = {
    mu_global    = Dim(1)
    sigma_global = Dim(1, lower = 0)
    alpha        = Dim(J, random = TRUE)
  },
  model = {
    alpha        ~ normal(mu_global, sigma_global)
    # ... likelihood ...
  }
)
```

このようなモデルでは、Laplace 近似を使って random effect
を周辺化できます。

``` r
fit_map_laplace  <- mdl$optimize(laplace = TRUE)
fit_mcmc_laplace <- mdl$sample(sampling = 1000, warmup = 1000, chains = 4, laplace = TRUE)
```

### ラッパー関数から始める

BayesRTMB には、標準的なモデルを簡単に使うためのラッパー関数があります。

たとえば線形回帰なら次のように書けます。

``` r
fit_lm <- rtmb_lm(mpg ~ wt + cyl, data = mtcars)

map_lm <- fit_lm$optimize()
map_lm$summary()

mcmc_lm <- fit_lm$sample(sampling = 1000, warmup = 1000, chains = 4)
mcmc_lm$summary()
```

一般化線形混合モデルは
[`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md)、t
検定は
[`rtmb_ttest()`](https://norimune.github.io/BayesRTMB/reference/rtmb_ttest.md)、因子分析は
[`rtmb_fa()`](https://norimune.github.io/BayesRTMB/reference/rtmb_fa.md)
のように、分析目的に応じた関数が用意されています。

### どのページから読むとよいか

初めて使う場合は、次の順番がわかりやすいです。

1.  この日本語紹介ページで全体像をつかむ
2.  `Reference` で
    [`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md),
    [`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md),
    [`Dim()`](https://norimune.github.io/BayesRTMB/reference/Dim.md)
    を確認する
3.  ラッパー関数や example を見て、自分の分析に近い形から試す

### 次の候補

日本語ページを今後増やすなら、次の 3 つが特に有用です。

- **日本語クイックスタート**
- **モデル記法の日本語解説**
- **NUTS / ADVI / MAP の使い分け**

これらを追加すると、英語の Reference
を読む前に全体像を理解しやすくなります。
