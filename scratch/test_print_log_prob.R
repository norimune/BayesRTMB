library(BayesRTMB)
devtools::load_all(".")

# テストデータ
set.seed(123)
N <- 5
X <- rnorm(N)
Y <- 2 + 0.5 * X + rnorm(N, 0, 1)

mdl <- rtmb_model(
  data = list(Y = Y, X = X),
  code = rtmb_code(
    model = {
      mu <- alpha + beta * X
      Y ~ normal(mu, sigma)
    },
    parameters = {
      alpha = Dim(0); beta = Dim(0); sigma = Dim(lower=0)
    }
  )
)

# 新しく追加したメソッドを呼び出し
mdl$print_log_prob()
