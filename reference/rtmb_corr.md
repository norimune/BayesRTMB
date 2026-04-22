# 相関行列 (多変量正規分布) 推定ラッパー

観測データから多変量正規分布を仮定して相関行列（および平均・標準偏差）を推定します。
観測変数が2つの場合は、直接スカラーの相関係数 (\`corr\`)
を推定するよう自動で切り替わります。

## Usage

``` r
rtmb_corr(
  data,
  prior = list(lkj_eta = 1, mu_sd = 10, sigma_rate = 1),
  init = NULL,
  null = NULL
)
```

## Arguments

- data:

  観測データフレームまたは行列 (N x P)。

- prior:

  事前分布のハイパーパラメータのリスト。デフォルトは
  `list(lkj_eta = 1.0, mu_sd = 10, sigma_rate = 1.0)`。

- init:

  初期値のリスト（省略可能）。

- null:

  帰無モデルを作成する際のターゲット（例: `"corr"`）。省略可能。

## Value

`RTMB_Model` クラスのインスタンス。
