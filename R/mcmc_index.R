map_est <- function(z) {
  if (abs(max(z) - min(z)) < 1e-10) return(z[1])
  d <- density(z)
  d$x[which.max(d$y)]
}

quantile95 <- function(x){
  quantile(x, probs = c(0.025, 0.975), names = TRUE)
}

ess_bulk <- function(x){
  posterior::ess_bulk(x)
}

ess_tail95 <- function(x){
  min(posterior::ess_quantile(x, probs = c(0.025, 0.975)))
}

r_hat <- posterior::rhat
