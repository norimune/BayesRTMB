# ユーザーのコード(AST)を探索し、パッケージの関数に名前空間を自動付与する関数

ユーザーのコード(AST)を探索し、パッケージの関数に名前空間を自動付与する関数

## Usage

``` r
inject_namespace(expr, pkg = "BayesRTMB")
```

## Arguments

- expr:

  評価前のコード(AST)

- pkg:

  対象のパッケージ名
