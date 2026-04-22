# MCMCフィットオブジェクトの周辺効果を計算する

MCMCフィットオブジェクトの周辺効果を計算する

## Usage

``` r
# S3 method for class 'mcmc_fit'
conditional_effects(fit, effect, resolution = 100, prob = 0.95, ...)
```

## Arguments

- fit:

  \`MCMC_Fit\` クラスのオブジェクト。

- effect:

  効果を可視化したい説明変数の名前（例: "X1" や "X1:X2"）。

- resolution:

  連続変数の場合に計算するグリッドの解像度（デフォルトは100）。

- prob:

  信用区間の確率（デフォルトは0.95）。

- ...:

  その他の引数。
