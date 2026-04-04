#' Wrapper function to create an RTMB_Model instance
#'
#' @param data A list of data used for modeling.
#' @param par_list A list of parameters defined by Dim().
#' @param log_prob A function to calculate the log-posterior probability (or log-likelihood + prior).
#' @param generate An optional function to calculate generated quantities.
#'
#' @return An instance of the RTMB_Model class.
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

#' Register a custom likelihood function to the lpdf environment
#'
#' @param name The name of the likelihood function to register (character string).
#' @param fun The function to calculate the likelihood.
#' @param force Logical; whether to overwrite an existing function with the same name (default is FALSE).
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

#' Register a custom mathematical function to the math environment
#'
#' @param name The name of the mathematical function to register (character string).
#' @param fun The function to register.
#' @param force Logical; whether to overwrite an existing function with the same name (default is FALSE).
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
