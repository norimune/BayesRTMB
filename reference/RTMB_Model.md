# RTMB_Model インスタンスを生成するラッパー関数

ユーザーが定義したデータ、およびモデルコード（尤度や事前分布）を結合し、
ベイズ推論（MCMC, VB, MAP推定）を実行するための \`RTMB_Model\`
(R6クラス) インスタンスを生成します。

## Usage

``` r
rtmb_model(
  data,
  code,
  par_names = list(),
  init = NULL,
  view = NULL,
  null_target = NULL
)
```

## Arguments

- data:

  モデルで使用する観測データや定数（サンプルサイズなど）を含む名前付きリスト。

- code:

  \`rtmb_code(...)\` で記述されたモデル定義ブロック（data, parameters,
  model, transform, generate を含む）。

- par_names:

  各パラメータの次元に対応する具体的な変数名のリスト（省略可能）。

- init:

  パラメータの初期値のリストまたは数値ベクトル（省略可能）。指定しない場合は自動でランダムに初期化されます。

- view:

  \`summary()\`
  などで結果を出力する際、優先して上部に表示したいパラメータ名の文字ベクトル（省略可能）。

- null_target:

  Character string.
  帰無モデルを同時に作成する場合に、固定対象のパラメータと無効化する事前分布をフォーミュラ形式の文字列で指定します（例:
  `"delta ~ cauchy(0, r)"`）。指定しない場合は `NULL`（デフォルト）。

## Value

モデルのコンパイルと事前テストが完了した \`RTMB_Model\`
クラスのインスタンス。
