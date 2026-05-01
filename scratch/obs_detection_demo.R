# BayesRTMB: 独立したデータポイントの自動検知デモ
library(RTMB)

# 1. データ追跡クラス (内部で尤度の要素数と値を記録)
DataTracker <- R6::R6Class("DataTracker",
  public = list(
    log_lik_list = list(),
    
    # 記録をリセット
    reset = function() {
      self$log_lik_list <- list()
    },
    
    # 尤度計算の結果を記録
    add = function(ll_vector) {
      self$log_lik_list[[length(self$log_lik_list) + 1]] <<- as.numeric(ll_vector)
      return(ll_vector)
    },
    
    # 合計観測数 (n_obs) を取得
    get_n_obs = function() {
      sum(sapply(self$log_lik_list, length))
    },
    
    # ポイントごとの対数尤度を取得 (WAIC/LOO用)
    get_pointwise_ll = function() {
      unlist(self$log_lik_list)
    }
  )
)

tracker <- DataTracker$new()

# 2. チルダ (~) の動作を模した関数
# BayesRTMBの内部変換では、LHS ~ RHS(params) をこれに変換するイメージ
model_link <- function(ll) {
  tracker$add(ll)
  return(-sum(ll)) # 負の対数尤度を返す
}

# --- デモ開始 ---

cat("--- Case 1: Vectorized (Normal) ---\n")
tracker$reset()
Y_vec <- rnorm(10) # 10個のデータ
mu <- 0; sigma <- 1
# model = { Y_vec ~ normal(mu, sigma) } の変換後イメージ
nll <- model_link(dnorm(Y_vec, mu, sigma, log = TRUE))
cat("Detected n_obs:", tracker$get_n_obs(), "\n") # 10であるべき


cat("\n--- Case 2: Loop (Long format) ---\n")
tracker$reset()
Y_long <- rnorm(5)
# model = { for(i in 1:5) Y_long[i] ~ normal(mu, sigma) } の変換後イメージ
nll_total <- 0
for (i in 1:5) {
  nll_total <- nll_total + model_link(dnorm(Y_long[i], mu, sigma, log = TRUE))
}
cat("Detected n_obs:", tracker$get_n_obs(), "\n") # 5であるべき


cat("\n--- Case 3: Multivariate (Single Unit) ---\n")
tracker$reset()
Y_mv <- matrix(rnorm(3), 1, 3) # 3次元のベクトル1つ
Mu_mv <- rep(0, 3); Sigma_mv <- diag(3)
# model = { Y_mv ~ mvnormal(Mu_mv, Sigma_mv) } の変換後イメージ
# dmvnorm はベクトルに対して「1つの」対数尤度を返す
nll_mv <- model_link(dmvnorm(Y_mv, Mu_mv, Sigma_mv, log = TRUE))
cat("Detected n_obs:", tracker$get_n_obs(), "\n") # 1であるべき (独立した単位は1つ)


cat("\n--- ポイントごとの対数尤度 (WAIC/CV用) ---\n")
# これらが取れることで、ユーザーがコードを書かなくてもWAICなどが計算可能になる
print(tracker$get_pointwise_ll())
