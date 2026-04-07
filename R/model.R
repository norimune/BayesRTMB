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
#' @param env Environment to assign to the generated function.
#' @return A standard R function object taking (dat, par).
#' @export
model_code <- function(expr, env = parent.frame()) {
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

        # 環境に存在する方を優先的に採用
        if (exists(name_lpdf, mode = "function", envir = env)) {
          actual_name <- name_lpdf
        } else if (exists(name_lpmf, mode = "function", envir = env)) {
          actual_name <- name_lpmf
        } else {
          actual_name <- name_lpdf # デフォルト
        }

        # lp <- lp + actual_name(target, args...) の純粋な構文木を作成
        new_call <- as.call(c(
          as.name("<-"),
          as.name("lp"),
          as.call(c(
            as.name("+"),
            as.name("lp"),
            as.call(c(as.name(actual_name), target, dist_args))
          ))
        ))
        return(new_call)
      }
      x[] <- lapply(x, rewrite_formula)
    }
    return(x)
  }

  processed_expr <- rewrite_formula(raw_expr)

  # ブロック {} の場合は中身を展開してフラットにする
  if (is.call(processed_expr) && identical(processed_expr[[1]], as.name("{"))) {
    expr_elements <- as.list(processed_expr)[-1]
  } else {
    expr_elements <- list(processed_expr)
  }

  # 関数の中身（ボディ）を構築: { getAll(dat, par); lp <- 0; ... ; return(lp) }
  body_list <- c(
    list(as.name("{")),
    quote(getAll(dat, par)),
    quote(lp <- 0),
    expr_elements,
    quote(return(lp))
  )
  body_expr <- as.call(body_list)

  # ---------------------------------------------------------
  # ここが最重要ポイント：純粋な R の function オブジェクトを構築する
  # ---------------------------------------------------------
  # 引数リスト `(dat, par)` を作成
  args <- as.pairlist(alist(dat = , par = ))

  # call("function", args, body) で関数定義の構文木を作り、evalで実体化
  fn_expr <- call("function", args, body_expr)
  log_prob_fn <- eval(fn_expr, envir = env)

  return(log_prob_fn)
}
