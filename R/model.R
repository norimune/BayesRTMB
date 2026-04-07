#' Wrapper function to create an RTMB_Model instance
#'
#' @param data A list of data used for modeling.
#' @param parameters A list (or alist) of parameters defined by Dim().
#' @param model An unquoted expression block defining the model. Optional if log_prob is provided.
#' @param transformed An optional function to calculate transformed parameters.
#' @param generate An optional function to calculate generated quantities.
#' @param log_prob An optional function to calculate the log-posterior probability.
#'
#' @return An instance of the RTMB_Model class.
#' @export
rtmb_model <- function(data, parameters,
                       model = NULL,
                       transformed = NULL,
                       generate = NULL,
                       log_prob = NULL) {

  # 1. log_prob の決定 (独立した model_code 関数を再利用)
  raw_model <- substitute(model)

  if (!is.null(log_prob)) {
    # ユーザーが直接 log_prob を指定した場合はそのまま採用
    log_prob_fn <- log_prob
  } else if (!is.null(raw_model)) {
    # 取得したASTを model_code に渡し、log_prob 関数を動的生成
    log_prob_fn <- eval(bquote(model_code(.(raw_model))))
  } else {
    stop("Either 'model' or 'log_prob' must be provided.")
  }

  # 2. parameters の内部評価 (P などの data 内変数を展開)
  build_par_list <- function(dat) {
    eval_env <- list2env(dat, parent = parent.frame())
    lapply(parameters, function(x) {
      if (typeof(x) %in% c("language", "symbol")) {
        eval(x, envir = eval_env)
      } else {
        x
      }
    })
  }

  evaluated_par_list <- build_par_list(data)

  # 3. RTMB_Model クラスのインスタンス化
  model_instance <- RTMB_Model$new(
    data         = data,
    par_list     = evaluated_par_list,
    log_prob     = log_prob_fn,
    transform    = transformed,
    generate     = generate
  )

  return(model_instance)
}
#' Model Code Wrapper for RTMB
#'
#' @param expr A block of code containing model description.
#' @return A log_prob function taking (dat, par).
#' @export
model_code <- function(expr) {
  raw_expr <- substitute(expr)

  # 式を再帰的に走査し、`~` の記述を書き換える
  rewrite_formula <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], as.name("~"))) {
        target <- x[[2]]
        dist_call <- x[[3]]
        dist_name <- as.character(dist_call[[1]])
        dist_args <- as.list(dist_call[-1])

        name_lpdf <- paste0(dist_name, "_lpdf")
        name_lpmf <- paste0(dist_name, "_lpmf")

        # lpdf と lpmf のうち、環境に存在する方を採用
        if (exists(name_lpdf, mode = "function")) {
          actual_name <- name_lpdf
        } else if (exists(name_lpmf, mode = "function")) {
          actual_name <- name_lpmf
        } else {
          stop(sprintf("Distribution function not found: %s or %s", name_lpdf, name_lpmf))
        }

        # lp <- lp + actual_name(target, args...) に書き換え
        new_call <- bquote(lp <- lp + .(as.call(c(as.name(actual_name), target, dist_args))))
        return(new_call)
      }
      # それ以外は再帰的に探索
      x[] <- lapply(x, rewrite_formula)
    }
    return(x)
  }

  processed_expr <- rewrite_formula(raw_expr)

  # 中括弧のネストを防ぐため、ブロック {} の場合は中身を展開する
  if (is.call(processed_expr) && identical(processed_expr[[1]], as.name("{"))) {
    expr_elements <- as.list(processed_expr)[-1]
  } else {
    expr_elements <- list(processed_expr)
  }

  # 最終的な関数の中身をフラットなリストとして結合し、関数ボディを構築
  body_list <- c(
    list(as.name("{")),
    quote(getAll(dat, par)),
    quote(lp <- 0),
    expr_elements,
    quote(return(lp))
  )
  body_expr <- as.call(body_list)

  # 関数オブジェクトとして生成
  log_prob_fn <- function(dat, par) {}
  body(log_prob_fn) <- body_expr

  return(log_prob_fn)
}
