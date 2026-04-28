pkgload::load_all(".")

cat("Fitting model with fixed effect interaction (wt * hp)...\n")
fit <- rtmb_glmer(mpg ~ wt * hp, data = mtcars, family = "gaussian")

cat("Fixed effect names:\n")
print(colnames(fit$raw_data)) # Wait, this is not right. 
# Let's check the AST parameters
print(fit$code$parameters)

cat("\nOptimizing...\n")
res <- fit$optimize()
res$summary()
