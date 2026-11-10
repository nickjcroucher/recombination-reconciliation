# Function to calculate predictions from fixed effects
predict_fixed <- function(beta, newdata) {
  
  # Extract design matrices for each structural parameter
  X_i <- model.matrix(~ 1, newdata)
  X_j <- model.matrix(~ 1, newdata)
  X_k <- model.matrix(~ 1, newdata)
  
  # Parameter estimates
  i <- as.vector(X_i %*% beta[grep("^i\\.", names(beta)) |
                                names(beta) == "i"])
  j <- as.vector(X_j %*% beta[grep("^j\\.", names(beta)) |
                                names(beta) == "j"])
  k <- as.vector(X_k %*% beta[grep("^k\\.", names(beta)) |
                                names(beta) == "k"])
  
  # Nonlinear response
  i - (i - j) * exp(k * newdata$Core)
}

# Prediction function using the fitted model's fixed effects
pred_fun <- function(beta, newdata, fit_mixed) {
  
  names(beta) <- names(fixef(fit_mixed))
  
  i <- beta["i"]
  j <- beta["j"]
  k <- beta["k"]
  
  i - (i - j) * exp(k * newdata$Core)
}

fit_nlme_model <- function(in_df) {
  
  start_vals <- c(
    i = 0.1,
    j = 0.01,
    k = -100
  )
  
  fit_mixed <- nlme(
    model = nl_expr,
    data = in_df,
    
    # FIXED EFFECTS: Define the linear models for each structural parameter
    fixed = list(
      i ~ 1, 
      j ~ 1, 
      k ~ 1
    ),
    
    # RANDOM EFFECTS: Species-dependent random deviations allowed on i and j
    random = i + j ~ 1 | species,
    
    start = start_vals,
    control = nlmeControl(
      maxiter = 200, 
      msMaxIter = 200, 
      pnlsMaxIter = 50,    # Give the PNLS step more attempts
      minScale = 1e-5,     # Lower the step-halving limit floor
      tolerance = 1e-3     # Slightly loosen precision required for convergence
    )
  )
  
}

make_nlme_prediction <- function(fit_mixed, in_df) {
  
  pred_grid <- expand_grid(
    Core = seq(
      min(in_df$Core, na.rm = TRUE),
      max(in_df$Core, na.rm = TRUE),
      length.out = 200
    )
  )
  
  # Population-level predictions
  pred_grid$response <- predict(
    fit_mixed,
    newdata = pred_grid,
    level = 0
  )
  
  # Fixed effects and covariance matrix
  beta <- fixef(fit_mixed)
  V <- vcov(fit_mixed)
  
  # Numerical Jacobian
  J <- jacobian(
    func = function(b) pred_fun(b, pred_grid, fit_mixed),
    x = beta
  )
  
  # Prediction variance
  pred_var <- diag(J %*% V %*% t(J))
  
  # Standard errors
  pred_se <- sqrt(pred_var)
  
  # 95% confidence intervals
  pred_grid <- pred_grid %>%
    mutate(
      asymp.LCL = response - 1.96 * pred_se,
      asymp.UCL = response + 1.96 * pred_se
    )
  
  return(pred_grid)
}