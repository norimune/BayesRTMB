#' RTMB_Modelインスタンスを生成するラッパー関数
#'
#' @param data モデリングに使用するデータのリスト
#' @param par_list Dim()で定義したパラメータのリスト
#' @param log_prob 対数事後確率（または対数尤度＋事前分布）を計算する関数
#' @param generate 生成量（Generated Quantities）を計算する関数（オプション）
#'
#' @return RTMB_Model クラスのインスタンス
#' @export
rtmb_model <- function(data, par_list, log_prob, generate=NULL) {
  model_instance <- RTMB_Model$new(
    data       = data,
    par_list   = par_list,
    log_prob   = log_prob,
    generate   = generate
  )
  return(model_instance)
}

#' カスタム尤度関数を lpdf 環境に登録する
#'
#' @param name 登録する尤度関数の名前（文字列）
#' @param fun 尤度を計算する関数
#' @param force すでに同名の関数が存在する場合に上書きするかどうか（デフォルトはFALSE）
#'
#' @export
register_lpdf <- function(name, fun, force = FALSE) {
  if (exists(name, envir = lpdf, inherits = FALSE)) {
    if (!force) {
      stop(sprintf("エラー: '%s' はすでに lpdf に登録されています。\n上書きするには force = TRUE を指定してください。", name))
    } else {
      packageStartupMessage(sprintf("注意: 既存の lpdf$%s を上書きしました。", name))
    }
  }

  if (!is.function(fun)) {
    stop("エラー:\n登録しようとしているオブジェクトは関数ではありません。")
  }

  assign(name, fun, envir = lpdf)
  message(sprintf("尤度関数 'lpdf$%s' が登録されました。", name))
}

#' カスタム数学関数を math 環境に登録する
#'
#' @param name 登録する数学関数の名前（文字列）
#' @param fun 登録する関数
#' @param force すでに同名の関数が存在する場合に上書きするかどうか（デフォルトはFALSE）
#'
#' @export
register_math <- function(name, fun, force = FALSE) {
  if (exists(name, envir = math, inherits = FALSE)) {
    if (!force) {
      stop(sprintf("エラー: '%s' はすでに math に登録されています。\n上書きするには force = TRUE を指定してください。", name))
    } else {
      packageStartupMessage(sprintf("注意: 既存の math$%s を上書きしました。", name))
    }
  }

  if (!is.function(fun)) {
    stop("エラー: 登録しようとしているオブジェクトは関数ではありません。")
  }

  assign(name, fun, envir = math)
  message(sprintf("数学関数 'math$%s' が登録されました。", name))
}
