  \donttest{
  # --- 1. Binary Data (e.g., correct/incorrect answers) ---
  # Simulate binary response data for 100 persons and 5 items
  set.seed(123)
  bin_data <- matrix(rbinom(500, size = 1, prob = 0.6), nrow = 100, ncol = 5)
  colnames(bin_data) <- paste0("Item", 1:5)

  # Introduce some missing values (NA) to demonstrate automatic handling
  #bin_data[sample(1:500, 10)] <- NA

  # Fit a 2-Parameter Logistic (2PL) model
  fit_2pl <- rtmb_irt(data = bin_data, model = "2PL", type = "binary")

  # Maximum A Posteriori (MAP) estimation
  map_2pl <- fit_2pl$optimize()
  map_2pl$summary()

  # --- 2. Ordered Data (e.g., Likert scale) ---
  # Simulate ordered response data (categories 1 to 5) with a true latent structure
  set.seed(123)
  N <- 100
  J <- 5
  K <- 4

  # True parameters
  theta <- rnorm(N, 0, 1)          # Person abilities (latent traits)
  a <- runif(J, 0.8, 2.0)          # Item discriminations

  # Item-specific thresholds (must be strictly increasing)
  b <- matrix(NA, nrow = J, ncol = K - 1)
  for (j in 1:J) {
    b[j, ] <- sort(c(-1.5, 0, 1.5) + rnorm(3, 0, 0.2))
  }

  ord_data <- matrix(NA, nrow = N, ncol = J)
  colnames(ord_data) <- paste0("Item", 1:J)

  # Generate responses based on the Graded Response Model
  for (i in 1:N) {
    for (j in 1:J) {
      # Latent continuous response (eta + standard logistic noise)
      eta <- a[j] * theta[i]
      y_star <- eta + rlogis(1)

      # Apply thresholds to determine the observed category (1 to 4)
      category <- 1
      for (k in 1:(K - 1)) {
        if (y_star > b[j, k]) category <- k + 1
      }
      ord_data[i, j] <- category
    }
  }

  # Fit a Graded Response Model (2PL for ordered data)
  fit_ord <- rtmb_irt(data = ord_data, model = "2PL", type = "ordered")
  fit_ord$print_code()

  map_ord <- fit_ord$optimize()
  map_ord$summary()

  # Note: For complex models like the Graded Response Model, the Wald CI from optimize()
  # may become extremely wide due to parameter transformations.
  # MCMC sampling is recommended for reliable interval estimation.

  # MCMC sampling for the ordered model (chains and iterations reduced)
  mcmc_ord <- fit_ord$sample(sampling = 500, warmup = 500, chains = 2)
  mcmc_ord$summary()
 }
