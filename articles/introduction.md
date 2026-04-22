# introduction

``` r
library(BayesRTMB)
#> Loading required package: RTMB
```

## BayesRTMB パッケージリファレンスマニュアル

### 1. モデルの定義と構築

#### `Dim(dim = 1, type = NULL, lower = NULL, upper = NULL, random = FALSE)`

パラメータの次元、型、制約、および変量効果（ランダム効果）の指定を行います。

- **引数:**
  - `dim`: パラメータの次元（スカラー=1, ベクトル=長さ, 行列=c(行, 列)
    など）。
  - `type`: パラメータの型（例: `"vector"`, `"matrix"`, `"ordered"`,
    `"simplex"`, `"corr_matrix"`, `"cov_matrix"`, `"CF_corr"`
    など）。省略時は次元から自動判別されます。
  - `lower`: パラメータの下限値。
  - `upper`: パラメータの上限値。
  - `random`: `TRUE`
    の場合、このパラメータを変量効果として扱い、Laplace近似の対象とします（MAP推定・MCMC共通）。
- **返り値:**
  制約された次元長（`length`）や無制約空間での次元長（`unc_length`）を含むリストを返します。

#### `rtmb_model(data, par_list, log_prob, generate = NULL)`

モデルを構築し、`RTMB_Model` インスタンスを生成します。

- **引数:**
  - `data`: モデルで使用するデータのリスト。
  - `par_list`:
    [`Dim()`](https://norimune.github.io/BayesRTMB/reference/Dim.md)
    関数で定義したパラメータのリスト。
  - `log_prob`:
    データとパラメータを引数に取り、対数事後確率（または対数尤度＋対数事前分布）を計算して返す関数。
  - `generate` (オプション): 生成量 (Generated Quantities)
    を計算する関数。
- **返り値:** `RTMB_Model` クラスのインスタンス。

------------------------------------------------------------------------

### 2. 確率分布と数学関数

パッケージには、`log_prob`
関数内で使用するためのヘルパー環境が用意されています。これらは自動微分
(AD) に対応した記述を容易にします。

#### `lpdf` 環境

対数確率密度関数 (log-pdf) および対数確率質量関数 (log-pmf)
を提供します。すべての関数は対数スケールでの合計値
(`sum(..., log=TRUE)`) を返します。 \* **主な分布:** `normal`,
`lognormal`, `exponential`, `beta`, `gamma`, `student_t`, `bernoulli`,
`binomial`, `poisson`, `ordered_logistic`, `multi_normal` など。 \*
**カスタム関数の追加:** `register_lpdf(name, fun, force = FALSE)`
を使用して、ユーザー独自の尤度関数を `lpdf` に登録できます。

#### `math` 環境

モデル構築に役立つ数学関数群を提供します。 \* **主な関数:** `inv_logit`,
`logit`, `log_sum_exp`, `log_mix`, `softmax`, `log_det_chol`,
`quad_form_chol` など。 \* **カスタム関数の追加:**
`register_math(name, fun, force = FALSE)`
を使用して、独自関数を登録できます。

------------------------------------------------------------------------

### 3. モデルクラス (`RTMB_Model`)

[`rtmb_model()`](https://norimune.github.io/BayesRTMB/reference/RTMB_Model.md)
によって生成される R6
クラスです。このオブジェクトから推定やサンプリングを実行します。

#### メソッド

##### `$optimize(laplace = TRUE, init = NULL, control = list())`

事後確率を最大化し、MAP (Maximum A Posteriori) 推定を行います。 \*
**引数:** \* `laplace`: `TRUE` の場合、`random = TRUE`
に設定されたパラメータを周辺化して最適化します。 \* `init`:
初期値のベクトル（省略時はランダム生成）。 \* **返り値:** `MAP_Fit`
クラスのインスタンス。

##### `$sample(sampling=1000, warmup=1000, chains=4, thin=1, seed=..., delta=0.8, max_treedepth=10, parallel=TRUE, laplace=FALSE, init=NULL)`

No-U-Turn Sampler (NUTS) を用いて事後分布から MCMC
サンプリングを行います。 \* **引数:** \* `sampling` / `warmup`:
サンプリングとウォームアップの反復回数。 \* `chains`:
実行するマルコフ連鎖の数。 \* `delta`: HMC/NUTSの目標受容確率。 \*
`parallel`: `TRUE` の場合、`future`
パッケージを利用して複数チェーンを並列実行します。 \* `laplace`: `TRUE`
の場合、変量効果をLaplace近似で周辺化した上でサンプリングを行います。 \*
**返り値:** `MCMC_Fit` クラスのインスタンス。

##### `$print_code()`

モデルの構造（パラメータリストと `log_prob`
関数のコード）をコンソールに出力します。

------------------------------------------------------------------------

### 4. フィット結果オブジェクト

#### `MAP_Fit` クラス

`$optimize()` の結果を格納します。 \* **フィールド:** `par`
(制約スケールのパラメータリスト), `objective` (負の対数事後確率),
`log_ml` (ラプラス近似による対数周辺尤度), `df_fixed` / `random_effects`
(推定値の要約データフレーム)。 \* **メソッド:** \*
`$summary(max_rows = 20, digits = 2)`:
推定値、標準誤差、95%信頼区間をまとめた表を出力します。

#### `MCMC_Fit` クラス

`$sample()` の結果を格納し、事後サンプルの抽出や事後処理を行います。 \*
**メソッド:** \* `$summary(pars = NULL, ...)`: Rhat、有効サンプルサイズ
(ESS)、事後平均、信用区間などの要約統計量をデータフレームで返します。 \*
`$draws(pars = NULL, chains = NULL, inc_random = FALSE, inc_gq = TRUE)`:
事後サンプルを3次元配列 `(iteration, chain, variable)` で抽出します。 \*
`$bridgesampling(method = "warp3", seed = NULL)`: Bridge Sampling
により周辺尤度（Marginal Likelihood）を推定します。 \*
`$generated_quantities(gq_fn = NULL)`:
抽出された事後サンプルに対して生成量（GQ）を計算し、オブジェクトに保存します。
\* `$rotate(target, linked = NULL, overwrite = TRUE)`:
指定した行列パラメータ（例:
因子負荷量）に対して直交プロクルステス回転を適用します。 \*
`$fa_rotate(loadings, scores = NULL, method = "promax", type = "oblique", ...)`:
`GPArotation` を用いて斜交・直交回転を適用します。 \*
`$resolve_switching(target, linked = NULL, overwrite = TRUE)`:
混合モデル等で発生するラベルスイッチング問題を、指定したターゲット変数の順序に基づいて解決・並べ替えます。

------------------------------------------------------------------------

### 5. 可視化・診断関数

MCMCサンプルの診断プロットを提供します。引数 `x` には `$draws()`
で抽出した3次元配列を渡します。

- **`plot_trace(x, mono = FALSE)`**:
  トレースプロットを描画し、チェーンの収束を確認します。
- **`plot_dens(x, mono = FALSE)`**:
  各チェーンの事後カーネル密度推定プロットを描画します。
- **`plot_acf(x, var_idx = 1)`**: 自己相関 (Autocorrelation)
  プロットを描画し、サンプルの自己相関の減衰を確認します。
