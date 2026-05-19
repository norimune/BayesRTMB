# ラッパー関数の使い方

BayesRTMB
では、[`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
を使ってモデルを直接書くこともできますが、よく使う分析についてはラッパー関数が用意されています。

ラッパー関数は、[`lm()`](https://rdrr.io/r/stats/lm.html) や
[`glm()`](https://rdrr.io/r/stats/glm.html)
のような短い書き方から、内部的に一貫した
[`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
を生成します。そのため、同じモデルに対して
MCMC、MAP推定、変分推論、そして頻度主義的分析を切り替えて使うことができます。

このページでは、各ラッパー関数の細かい理論やオプションではなく、「どのような雰囲気で使うか」を中心に紹介します。

## ラッパー関数とは

ラッパー関数は、典型的な統計モデルを簡潔に書くための入り口です。たとえば、重回帰なら次のように書けます。

``` r

library(BayesRTMB)
data(debate)

mdl <- rtmb_lm(sat ~ talk * perf, data = debate)
```

この時点では、まだ推定は実行されていません。モデルを作ったあとに、目的に応じて推定方法を選びます。

``` r

fit_mcmc <- mdl$sample()
fit_map  <- mdl$optimize()
fit_cl   <- mdl$classic()
```

同じ `mdl`
から、ベイズ推定と頻度主義的分析の両方を扱えるところが、BayesRTMB
のラッパー関数の大きな特徴です。

## ラッパー関数一覧

主なラッパー関数は次のとおりです。ここでは、細かなオプションよりも典型的な使い方を示します。

| 関数 | 分析 | 典型的なコード |
|:---|:---|:---|
| [`rtmb_ttest()`](https://norimune.github.io/BayesRTMB/reference/rtmb_ttest.md) | t検定 | `rtmb_ttest(sat ~ cond, data = debate)` |
| [`rtmb_lm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lm.md) | 線形回帰 | `rtmb_lm(sat ~ talk * perf, data = debate)` |
| [`rtmb_glm()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glm.md) | 一般化線形モデル | `rtmb_glm(y ~ x, data = df, family = "bernoulli")` |
| [`rtmb_lmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lmer.md) | 線形混合モデル | `rtmb_lmer(y ~ x + (1 | group), data = df)` |
| [`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md) | 一般化線形混合モデル | `rtmb_glmer(y ~ x + (1 | group), data = df, family = "poisson")` |
| [`rtmb_corr()`](https://norimune.github.io/BayesRTMB/reference/rtmb_corr.md) | 相関・偏相関 | `rtmb_corr(cbind(sat, perf), data = debate)` |
| [`rtmb_table()`](https://norimune.github.io/BayesRTMB/reference/rtmb_table.md) | クロス表分析 | `rtmb_table(skill, cond, data = debate)` |
| [`rtmb_loglinear()`](https://norimune.github.io/BayesRTMB/reference/rtmb_loglinear.md) | 対数線形モデル | `rtmb_loglinear(~ A + B + A:B, data = df)` |
| [`rtmb_fa()`](https://norimune.github.io/BayesRTMB/reference/rtmb_fa.md) | 因子分析 | `rtmb_fa(items, nfactors = 3)` |
| [`rtmb_irt()`](https://norimune.github.io/BayesRTMB/reference/rtmb_irt.md) | 項目反応理論 | `rtmb_irt(items, model = "2PL")` |
| [`rtmb_lrt()`](https://norimune.github.io/BayesRTMB/reference/rtmb_lrt.md) | 潜在ランク理論 | `rtmb_lrt(y, k = 3)` |
| [`rtmb_mdu()`](https://norimune.github.io/BayesRTMB/reference/rtmb_mdu.md) | 多次元展開法 | `rtmb_mdu(Y, ndim = 2)` |
| [`rtmb_mediation()`](https://norimune.github.io/BayesRTMB/reference/rtmb_mediation.md) | 媒介分析 | `rtmb_mediation(list(m ~ x, y ~ x + m), data = df)` |
| [`rtmb_mixture()`](https://norimune.github.io/BayesRTMB/reference/rtmb_mixture.md) | 混合分布モデル | `rtmb_mixture(y ~ x, k = 3, data = df)` |

## 推定のイメージ

ラッパー関数でモデルを作ったあとは、メソッドを選んで推定します。

MCMC を使う場合は [`sample()`](https://rdrr.io/r/base/sample.html)
です。

``` r

mdl <- rtmb_lm(sat ~ talk * perf, data = debate)
fit <- mdl$sample()
fit$summary()
```

MAP推定を使う場合は
[`optimize()`](https://rdrr.io/r/stats/optimize.html) です。

``` r

fit <- mdl$optimize()
fit$summary()
```

頻度主義的分析を使う場合は `classic()` です。

``` r

fit <- mdl$classic()
fit$summary()
```

このように、モデルの書き方は同じで、推定方法だけを切り替えられます。

## 頻度主義的分析

`classic()` は、BayesRTMB
のラッパー関数を頻度主義的な分析として使うためのメソッドです。

``` r

rtmb_corr(cbind(sat, perf), data = debate)$classic()
```

[`rtmb_corr()`](https://norimune.github.io/BayesRTMB/reference/rtmb_corr.md)
のデフォルトは `method = "pearson"` なので、`classic()` は
[`cor.test()`](https://rdrr.io/r/stats/cor.test.html) と同じ Pearson
相関の結果を返します。

``` text
## Call:
## Classical estimation via corr 
## 
## Point Estimates and Confidence Intervals:
##           Estimate Std. Error Lower 95% Upper 95%  df t value     Pr    
## corr[rho]  0.29742         NA   0.19059   0.39728 298 5.37755 <.0001 ***
```

統制変数を入れた場合は、共変量だけを統制した偏相関になります。

``` r

rtmb_corr(cbind(sat, perf),
          data = debate,
          covariates = ~ skill)$classic()
```

``` text
##           Estimate Std. Error Lower 95% Upper 95%  df t value     Pr    
## pcorr[rho]  0.29604         NA   0.18896   0.39617 297 5.34134 <.0001 ***
```

`method = "spearman"` を指定すると Spearman
相関を使います。`method = "reml"` を指定すると、RTMB
のREML推定に基づく相関モデルとして推定されます。

``` r

rtmb_corr(cbind(sat, perf), data = debate, method = "spearman")$classic()
rtmb_corr(cbind(sat, perf), data = debate, method = "reml")$classic()
```

### S3メソッド

`classic()` の結果は `Classic_Fit`
オブジェクトなので、[`anova()`](https://rdrr.io/r/stats/anova.html)、[`AIC()`](https://rdrr.io/r/stats/AIC.html)、[`BIC()`](https://rdrr.io/r/stats/AIC.html)
などのS3メソッドを使えます。

クロス表分析では、1つの対象に対して
[`anova()`](https://rdrr.io/r/stats/anova.html)
を使うと、独立性の検定を表示できます。

``` r

fit_tab <- rtmb_table(skill, cond, data = debate)$classic()
anova(fit_tab)
```

``` text
## Contingency Table Analysis Tests
## - 
##                                        X-squared df p-value
## Pearson's Chi-squared (from est. prob)   1.20478  2  0.5475
## 
## Fisher's Exact Test for Count Data
## p-value =  0.54310 
```

複数のモデルを比較する場合は、`anova(fit1, fit2)`
と書けます。たとえば、因子分析で因子数の異なるモデルを比較できます。

``` r

data(BigFive)

items <- BigFive[, 1:10]
fit_fa1 <- rtmb_fa(items, nfactors = 1)$classic()
fit_fa2 <- rtmb_fa(items, nfactors = 2)$classic()

anova(fit_fa1, fit_fa2)
```

``` text
## Likelihood Ratio Tests
##    Model Df  logLik    AIC    BIC  Chisq Chi Df Pr(>Chisq)
## fit_fa1 30 -2485.5 5030.9 4970.9     NA     NA         NA
## fit_fa2 40 -2440.4 4960.9 4880.9 90.089     10 5.1421e-15
```

[`AIC()`](https://rdrr.io/r/stats/AIC.html) や
[`BIC()`](https://rdrr.io/r/stats/AIC.html) も同じように使えます。

``` r

c(AIC_1 = AIC(fit_fa1),
  AIC_2 = AIC(fit_fa2),
  BIC_1 = BIC(fit_fa1),
  BIC_2 = BIC(fit_fa2))
```

``` text
##    AIC_1    AIC_2    BIC_1    BIC_2 
## 5030.942 4960.853 4970.942 4880.853
```

## print_code()

ラッパー関数は簡単に使えますが、内部ではすべて
[`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
に展開されています。どのようなモデルが生成されたかは `print_code()`
で確認できます。

``` r

mdl <- rtmb_ttest(sat ~ cond, data = debate)
mdl$print_code()
```

``` text
## === RTMB Model Code ===
## 
## rtmb_code(
##   setup = {
##     Y1 <- as.numeric(na.omit(Y[G == 1]))
##     Y2 <- as.numeric(na.omit(Y[G == 2]))
##   }, 
##   parameters = {
##     mean0 = Dim(1)
##     mean1 = Dim(1)
##     sd = Dim(2, lower = 0)
##   }, 
##   transform = {
##     diff <- mean0 - mean1
##     mean <- (mean0 + mean1)/2
##     sd_pooled <- sqrt((sd[1]^2 + sd[2]^2)/2)
##     delta <- diff/sd_pooled
##   }, 
##   model = {
##     Y1 ~ normal(mean0, sd[1])
##     Y2 ~ normal(mean1, sd[2])
##   }
## )
```

この例は Welch の
t検定に対応するモデルです。通常の検定としては短く呼び出せますが、内部では2群それぞれの平均と標準偏差を推定するモデルとして表現されています。

このように、ラッパー関数は「簡単に使う」ための入り口であると同時に、「どのモデルを推定しているか」を確認しながら、必要に応じて独自の
[`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
へ進むための足場にもなります。

## 欠損値の扱い (Missing Data Handling)

BayesRTMB のラッパー関数は、欠損値（`NA`）を適切かつ統一的に扱うための
`missing`
引数を備えています。分析手法やデータの形式に応じて、最適な手法を選択できます。

多くのラッパー関数ではデフォルトが
`"listwise"`（リストワイズ削除）となっていますが、多変量解析を扱う一部のラッパー関数では、**完全情報最尤法（FIML）**
や **ペアワイズ削除** も選択可能です。

| ラッパー関数 | 指定可能な `missing` オプション | デフォルト | 内容説明 |
|:---|:---|:---|:---|
| `rtmb_lm`, `rtmb_glm`, `rtmb_lmer`, `rtmb_glmer`, `rtmb_ttest` | `"listwise"` | `"listwise"` | モデル式（formula）に含まれる変数に `NA` が1つでもある行（観測値）をすべて除外します。 |
| `rtmb_fa`, `rtmb_corr` | `"listwise"`, `"fiml"` | `"listwise"` | `"listwise"` は完全データのみから共分散行列（十分統計量）を計算し高速に推定します。`"fiml"` は行ごとの欠損パターンに基づき、観測されている変数のみの周辺対数尤度を動的に足し合わせて推定します（完全情報最尤法）。 |
| `rtmb_corr` （追加オプション） | `"pairwise"` |  | `method = "pearson"` または `"spearman"` のときのみ使用可能。R標準の `cor(..., use = "pairwise.complete.obs")` と完全に一致するペアワイズ削除結果を `classic()` メソッドなどで出力します。 |
| `rtmb_mdu` | `"listwise"`, `"fiml"` | `"listwise"` | `"fiml"` は評定データの欠損値（`NA`）を尤度計算ループ内で安全にスキップし、有効なデータのみで推定します。`"listwise"` は欠損を含む対象者を丸ごと除外します。 |
| `rtmb_irt` | `"listwise"`, `"fiml"` | `"fiml"` | `"fiml"` は項目ごとの欠損値を自動的にスキップし、手元にある回答のみを用いて能力や項目パラメータを推定します（デフォルト）。`"listwise"` を指定すると、1項目でも欠損がある回答者をすべて除外します。 |

#### 分析例：因子分析における完全情報最尤法（FIML）

データに一部の欠損が含まれている場合、`missing = "fiml"`
を指定することで、一部しか回答していない参加者も含めた全データを用いて因子負荷量を頑健に推定できます。

``` r

data(BigFive)
# 意図的にランダムな欠損を生成
set.seed(42)
BigFive_miss <- BigFive
for (col in 1:ncol(BigFive_miss)) {
  BigFive_miss[sample(1:nrow(BigFive_miss), 10), col] <- NA
}

# FIMLによる探索的因子分析
mdl_fa <- rtmb_fa(BigFive_miss, nfactors = 5, missing = "fiml")
fit_fa <- mdl_fa$optimize()
```

#### 分析例：相関分析におけるペアワイズ削除（Pairwise）

R標準の [`cor()`](https://rdrr.io/r/stats/cor.html) や
[`cor.test()`](https://rdrr.io/r/stats/cor.test.html)
のペアワイズ削除の挙動と一致させたい場合は、`missing = "pairwise"`
を指定して `classic()` を実行します。

``` r

# ペアワイズ削除による古典的相関分析
mdl_corr <- rtmb_corr(BigFive_miss[, 1:5], method = "pearson", missing = "pairwise")
fit_corr <- mdl_corr$classic()
fit_corr$summary()
```

------------------------------------------------------------------------

## 次に読む記事

BayesRTMB
の全体像や、より細かな使い方を知りたい場合は、次の記事も参照してください。

- [イントロダクション](https://norimune.github.io/BayesRTMB/articles/ja-introduction.md):
  BayesRTMB
  の基本的な考え方と、MCMC、MAP推定、変分推論、頻度主義的分析の位置づけ。
- [クイックスタート](https://norimune.github.io/BayesRTMB/articles/ja-quick_start.md):
  インストール、最初のモデル、MCMC診断、可視化、t検定のベイズファクター。
- [階層モデル・GLMM・分散分析](https://norimune.github.io/BayesRTMB/articles/ja-rtmb_glmer.md):
  [`rtmb_glmer()`](https://norimune.github.io/BayesRTMB/reference/rtmb_glmer.md)
  を使って mixed model / GLMM を詳しく扱う方法。
- [モデルを書く](https://norimune.github.io/BayesRTMB/articles/ja-writing_models.md):
  [`rtmb_code()`](https://norimune.github.io/BayesRTMB/reference/rtmb_code.md)
  を使って独自のモデルを書く方法。
- [内部構造](https://norimune.github.io/BayesRTMB/articles/ja-rtmb_internals.md):
  `setup`、`parameters`、`transform`、`model` など、BayesRTMB
  の内部的なモデル構造。
