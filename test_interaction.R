pkgload::load_all(".")

mtcars$gear_factor <- as.factor(mtcars$gear)
mtcars$cyl_factor <- as.factor(mtcars$cyl)

cat("Fitting model with interaction using * operator...\n")
fit <- rtmb_glmer(mpg ~ wt * hp + (1|cyl_factor*gear_factor), data = mtcars, family = "gaussian")

cat("Generated AST:\n")
print(fit$code)

cat("\nOptimizing...\n")
res <- fit$optimize()
res$summary()
