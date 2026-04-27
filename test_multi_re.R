pkgload::load_all(".")

# Create a dummy variable for testing multiple random effects
mtcars$gear_factor <- as.factor(mtcars$gear)

# Fit model with two random effects: cyl and gear
cat("Fitting model...\n")
fit <- rtmb_glmer(mpg ~ wt + (1|cyl) + (1|gear_factor), data = mtcars, family = "gaussian")

cat("Generated AST:\n")
print(fit$code)

cat("\nOptimizing...\n")
res <- fit$optimize()
res$summary()
