# BayesRTMB クイックスタート

## このページの目的

このページでは、BayesRTMB を使い始めるための最小限の流れを確認します。

ここで扱うのは、次の内容です。

1.  [`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
    で小さなモデルを書く
2.  [`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md)
    でモデルオブジェクトを作る
3.  [`optimize()`](https://rdrr.io/r/stats/optimize.html) と
    [`sample()`](https://rdrr.io/r/base/sample.html) で推定する
4.  MCMC の事後分布を可視化する
5.  ラッパー関数で重回帰と交互作用を扱う
6.  `classic()` による頻度主義的な t 検定を行う
7.  JZS prior と MCMC によって Bayes factor を計算する

詳しいモデルコードの書き方や、混合モデル、GLMM、モデル比較の詳細は、別の
vignette で扱います。

## 0. インストールと環境確認

BayesRTMB は GitHub からインストールできます。 `remotes`
を使う場合は、次のようにします。

``` r

install.packages("remotes")
remotes::install_github("norimune/BayesRTMB")
```

`pak` を使う場合は、次のようにインストールできます。

``` r

install.packages("pak")
pak::pak("norimune/BayesRTMB")
```

### Windows ユーザー向け: Rtools の確認

BayesRTMB は内部で RTMB / TMB を利用します。 Windows
では、これらの依存パッケージをコンパイルするために Rtools が必要です。

Rtools が利用できるかどうかは、次のコードで確認できます。

``` r

pkgbuild::check_build_tools(debug = TRUE)
```

`TRUE`
が返る、またはビルドツールが利用可能であることが表示されれば、準備はできています。
Rtools が見つからない場合は、利用している R のバージョンに対応する
Rtools をインストールし、R または RStudio を再起動してください。

## 1. 最小モデルを書く

まず、もっとも小さな例として二項モデルを書きます。 ここでは、10
回の試行のうち成功が 6 回観測された状況を考えます。

``` r

library(BayesRTMB)

Trial <- 10
Y <- 6

dat <- list(Trial = Trial, Y = Y)

code <- rtmb_code(
  parameters = {
    theta <- Dim(lower = 0, upper = 1)
  },
  model = {
    Y ~ binomial(Trial, theta)
    theta ~ beta(1, 1)
  }
)
```

`parameters` ブロックでは推定するパラメータを宣言します。
ここでは成功確率 `theta` を、0 から 1
の範囲を持つパラメータとして定義しています。

`model` ブロックでは、観測データの分布と事前分布を書きます。
`Y ~ binomial(Trial, theta)` は、成功数 `Y`
が二項分布に従うことを表しています。

## 2. モデルオブジェクトを作る

[`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md)
にデータとモデルコードを渡すと、推定用のモデルオブジェクトが作られます。

``` r

mdl <- rtmb_model(dat, code)
```

``` text
## Pre-checking model code...
## Checking RTMB setup...
```

この段階では、まだ推定は行われていません。 `mdl`
は、モデル定義とデータを保持した `RTMB_Model` オブジェクトです。

## 3. MAP 推定を行う

点推定をすばやく得たいときは、[`optimize()`](https://rdrr.io/r/stats/optimize.html)
を使います。 BayesRTMB では、事前分布がある場合は MAP 推定、flat prior
の場合は最尤推定に近い推定として扱えます。

``` r

fit_map <- mdl$optimize()
fit_map
```

``` text
## Starting RTMB optimization...
## 
## 
## Call:
## MAP Estimation via RTMB
## 
## Negative Log-Posterior: 1.38
## Approx. Log Marginal Likelihood (Laplace): -2.33
## 
## Point Estimates and 95% Wald CI:
## variable  Estimate  Std. Error  Lower 95%  Upper 95% 
## theta      0.60000     0.15492    0.29740    0.84166 
```

この例では、成功確率 `theta` の推定値はちょうど 0.60 です。

## 4. MCMC で事後分布を見る

事後分布全体を見たいときは、[`sample()`](https://rdrr.io/r/base/sample.html)
を使います。
クイックスタートでは短い設定にしていますが、実際の分析ではより多くの反復数を使ってください。

``` r

set.seed(1)

fit_mcmc <- mdl$sample(
  sampling = 200,
  warmup = 200,
  chains = 2
)

fit_mcmc$summary()
```

``` text
## Starting sequential sampling (chains = 2)...
## chain 1 started... 
## chain 1: iter 100 warmup 
## chain 1: iter 200 warmup 
## chain 1: iter 300 sampling 
## chain 1: iter 400 sampling 
## chain 2 started... 
## chain 2: iter 100 warmup 
## chain 2: iter 200 warmup 
## chain 2: iter 300 sampling 
## chain 2: iter 400 sampling 
## variable   mean    sd    map   q2.5  q97.5  ess_bulk  ess_tail  rhat 
## lp        -3.30  0.69  -2.88  -5.25  -2.80       132       188  1.00 
## theta      0.58  0.13   0.60   0.32   0.83       165       177  1.01 
```

MCMC の結果では、平均、標準偏差、事後分位点、ESS、R-hat
などを確認できます。 R-hat が 1 に近く、ESS が十分に大きいほど、MCMC
の診断としては安心しやすくなります。

### MCMC を並列化する場合

通常の MCMC は追加設定なしで実行できます。
ただし、`sample(parallel = TRUE)` のように MCMC
を並列実行する場合は、追加で `future`, `future.apply`, `progressr`
が必要です。 `progressr` は、内部で
[`progressr::progressor()`](https://progressr.futureverse.org/reference/progressor.html)
による進捗表示にも使われます。

``` r

install.packages(c("future", "future.apply", "progressr"))
```

これらは BayesRTMB の Suggests
パッケージなので、並列化を使わない場合は必須ではありません。

``` r

fit_mcmc <- mdl$sample(
  sampling = 1000,
  warmup = 1000,
  chains = 4,
  parallel = TRUE
)
```

## 5. MCMC の結果を可視化する

MCMC の結果は、数値だけでなく図でも確認します。 `draws()`
で対象パラメータのサンプルを取り出し、密度、トレース、自己相関、区間推定を確認できます。

``` r

theta_draws <- fit_mcmc$draws("theta")

plot_dens(theta_draws)
plot_trace(theta_draws)
plot_acf(theta_draws)
plot_forest(theta_draws)
```

それぞれの図は、次の目的で使います。

- [`plot_dens()`](https://norimune.github.io/BayesRTMB/reference/plot_dens.md)
  は、事後分布の形を確認します。
- [`plot_trace()`](https://norimune.github.io/BayesRTMB/reference/plot_trace.md)
  は、チェインがよく混ざっているかを確認します。
- [`plot_acf()`](https://norimune.github.io/BayesRTMB/reference/plot_acf.md)
  は、自己相関が強すぎないかを確認します。
- [`plot_forest()`](https://norimune.github.io/BayesRTMB/reference/plot_forest.md)
  は、点推定値と信用区間を一覧します。

実際には、たとえば次のような図として確認できます。

![Posterior density plot](plot_dens.png)

Posterior density plot

![Trace plot](plot_trace.png)

Trace plot

![Autocorrelation plot](plot_acf.png)

Autocorrelation plot

![Forest plot](plot_forest.png)

Forest plot

## 6. ラッパー関数で重回帰を行う

標準的な分析では、[`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
を自分で書かずにラッパー関数から始められます。

ここでは、`debate` データを使い、満足度 `sat` を
`talk`、`perf`、およびその交互作用で説明する重回帰モデルを推定します。

``` r

data(debate)

mdl_lm <- rtmb_lm(
  sat ~ talk * perf,
  data = debate,
  prior = prior_normal()
)

fit_lm <- mdl_lm$optimize(
  se_method = "sampling",
  num_samples = 1000,
  seed = 1
)
fit_lm
```

``` text
## Starting RTMB optimization...
## 
## Using simulation-based error propagation (1000 samples)...
## 
## 
## Call:
## MAP Estimation via RTMB
## 
## Negative Log-Posterior: 402.24
## Approx. Log Marginal Likelihood (Laplace): -414.02
## 
## Point Estimates and 95% Sampling-based CI:
##     variable  Estimate  Std. Error  Lower 95%  Upper 95% 
## Intercept_c    3.43325     0.05232    3.33254    3.53601 
## b[talk]       -0.34490     0.15982   -0.65456   -0.01513 
## b[perf]       -0.25266     0.10109   -0.45514   -0.04584 
## b[talk:perf]   0.13020     0.03059    0.06865    0.19039 
## sigma          0.87135     0.03637    0.80450    0.94693 
## Intercept      3.79883     0.51483    2.75177    4.83008 
```

係数を一覧したいときは、[`plot_forest()`](https://norimune.github.io/BayesRTMB/reference/plot_forest.md)
が便利です。

``` r

fit_lm$draws(c("b[talk]", "b[perf]", "b[talk:perf]")) |>
  plot_forest(point_estimate = "MAP")
```

![Regression coefficient forest plot](plot_forest.png)

Regression coefficient forest plot

## 7. 交互作用を図で確認する

交互作用は、係数表だけでは解釈しにくいことがあります。
その場合は、[`conditional_effects()`](https://norimune.github.io/BayesRTMB/reference/conditional_effects.md)
を使って予測値を図示します。

``` r

ce <- conditional_effects(fit_lm, effect = "talk:perf")
plot(ce)
```

![Conditional effect plot](conditional_effect.png)

Conditional effect plot

`effect = "talk:perf"` と書くと、`talk` の効果が `perf`
の値によってどのように変わるかを確認できます。

より詳しく調べたい場合は、単純傾斜を計算する
[`simple_effects()`](https://norimune.github.io/BayesRTMB/reference/simple_effects.md)
も使えます。

``` r

simple_effects(fit_lm, effect = "talk:perf")
```

## 8. t 検定を頻度主義的に行う

BayesRTMB のラッパー関数は、頻度主義的な分析にも使えます。
[`prior_flat()`](https://norimune.github.io/BayesRTMB/reference/prior_flat.md)
を使った t 検定に対して `classic()` を呼ぶと、通常の t
検定に近い形で結果を表示できます。

``` r

mdl_t <- rtmb_ttest(
  sat ~ cond,
  data = debate,
  prior = prior_flat()
)

fit_t_classic <- mdl_t$classic()
fit_t_classic
```

``` text
## Pre-checking model code...
## Checking RTMB setup...
## Starting RTMB optimization...
## 
## 
## Call:
## Classical estimation via ttest 
## 
## Log-Likelihood: -421.320, AIC: 848.640, BIC: 842.640
## 
## Point Estimates and Confidence Intervals:
##       Estimate Std. Error Lower 95% Upper 95%  df  t value     Pr    
## diff  -0.37333    0.11297  -0.59564  -0.15102 298 -3.30484 .00107  **
## delta -0.38161    0.11652  -0.61092  -0.15230 298 -3.27497 .00118  **
## mean   3.43333    0.05648   3.32218   3.54449 298 60.78547 <.0001 ***
## sd     0.97831    0.04007   0.90254   1.06044 298 -0.53533 .59282    
## mean0  3.24667    0.07988   3.08947   3.40386 298 40.64495 <.0001 ***
## mean1  3.62000    0.07988   3.46280   3.77720 298 45.31870 <.0001 ***
```

ここでは、`diff` が2群の平均差、`delta` が標準化効果量を表します。
`classic()` は、BayesRTMB
のモデルを頻度主義的な推定として確認したいときに使えます。

## 9. JZS prior で Bayes factor を計算する

同じ t 検定でも、JZS prior を使うと、効果量 `delta` に Cauchy
事前分布を置いた Bayes factor を計算できます。

``` r

mdl_t_jzs <- rtmb_ttest(
  sat ~ cond,
  data = debate,
  prior = prior_jzs()
)

set.seed(2)

fit_t_jzs <- mdl_t_jzs$sample(
  sampling = 200,
  warmup = 200,
  chains = 2
)

bf <- fit_t_jzs$bayes_factor(fixed = list(delta = 0))
bf
```

``` text
## Calculating marginal likelihood for the full model...
## Bridge Sampling Converged: LogML = -424.669 (Error = 0.0180, ESS = 56.0)
## 
## --- Sampling from the comparison model ---
## Starting sequential sampling (chains = 2)...
## 
## --- Calculating marginal likelihood for the comparison model ---
## Bridge Sampling Converged: LogML = -427.732 (Error = 0.0111, ESS = 175.7)
## --- Bayes Factor Analysis (Bridge Sampling) ---
## Bayes Factor (BF12) : 21.3858 
## Log Bayes Factor    : 3.0627 (Approx. Error = 0.0212)
## Evidence            : Strong evidence for Model 1 
## Comparison model    : Parameters fixed at list(delta = 0) 
```

`fixed = list(delta = 0)` は、効果量を 0
に固定した帰無モデルと比較する指定です。
この例では、`Model 1`、つまり効果量を推定するモデルのほうが支持されています。

実際の分析では、Bayes factor
の安定性を確認するために、ここで示した例よりも多い MCMC
サンプル数を使うことをおすすめします。

## 次に読むページ

このページでは、BayesRTMB の入口だけを扱いました。
目的に応じて、次のページに進んでください。

1.  **[ラッパー関数の使い方](https://norimune.github.io/BayesRTMB/articles/ja-wrapper_functions.md)**  
    回帰、GLM、混合モデル、t 検定、相関、因子分析、IRT
    など、標準的な分析をラッパー関数で行う方法を確認できます。

2.  **[モデルコードの書き方](https://norimune.github.io/BayesRTMB/articles/ja-writing_models.md)**  
    [`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
    の `setup`, `parameters`, `transform`, `model`, `generate`
    を使って、独自モデルを書く方法を学べます。

3.  **[RTMB
    の仕組みと推定アルゴリズム](https://norimune.github.io/BayesRTMB/articles/ja-rtmb_internals.md)**  
    MAP 推定、Laplace 近似、MCMC、変分推論などの内部処理を確認できます。
