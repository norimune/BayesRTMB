# 安全なRTMBモデルの構築（エラーメッセージ翻訳付き）

\`rtmb_model\` をラップし、内部の \`MakeADFun\` 実行時に発生する
C++/RTMB 由来の
難解なエラーメッセージを、ユーザーが理解しやすい日本語のヒントに翻訳して出力します。

## Usage

``` r
safe_rtmb_model(
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

  \`validate_data\` を通過したデータのリスト。

- code:

  \`model_code()\` で定義された尤度や事前分布のコードブロック。

- par_names:

  各パラメータの次元に対応する変数名のリスト（省略可能）。

- init:

  パラメータの初期値のリストまたは数値ベクトル（省略可能）。

- view:

  結果のsummary出力時に優先して表示したいパラメータ名の文字ベクトル（省略可能）。

- null_target:

  帰無モデルを同時に作成する場合に、固定対象のパラメータと無効化する事前分布を指定する文字列（省略可能）。

## Value

正常にコンパイルが完了した場合は、\`rtmb_model\`
オブジェクト（R6クラスのインスタンス）を返します。
