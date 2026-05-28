
  # Create a contingency table
  tab <- matrix(c(10, 20, 30, 40), nrow = 2)
  dimnames(tab) <- list(A = c("A1", "A2"), B = c("B1", "B2"))
  
  # Fit a log-linear model (independence model: ~ A + B)
  fit_log <- rtmb_loglinear(~ A + B, data = tab)
  
  # MAP estimation
  map_log <- fit_log$optimize()
  map_log$summary()
