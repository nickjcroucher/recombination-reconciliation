# Test for significance of differences between models
parameter_pval <- function(model1 = NULL,model2 = NULL,param = NULL) {
  # Extract fixed effects coefficients and standard errors from both models
  beta1 <- fixef(model1)[param]
  beta2 <- fixef(model2)[param]
  
  se1 <- summary(model1)$tTable[param, "Std.Error"]
  se2 <- summary(model2)$tTable[param, "Std.Error"]
  
  # Calculate the Z-statistic
  z_stat <- (beta1 - beta2) / sqrt(se1^2 + se2^2)
  
  # Get the two-tailed p-value
  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  
  # Print results
  return(as.numeric(p_value))
}

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
  names(beta) <- names(nlme::fixef(fit_mixed))
  
  i <- beta["i"]
  j <- beta["j"]
  k <- beta["k"]
  
  as.vector(i - (i - j) * exp(k * newdata$Core))
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

make_nlme_prediction <- function(fit_mixed, in_df,
                                 prediction_level = 0,
                                 n_points = 200) {
  
  # Create prediction grid
  if (prediction_level == 0) {
    pred_grid <- data.frame(
      Core = seq(
        min(in_df$Core, na.rm = TRUE),
        max(in_df$Core, na.rm = TRUE),
        length.out = n_points
      )
    )
  } else {
    pred_grid <- expand_grid(
      Core = seq(
        min(in_df$Core, na.rm = TRUE),
        max(in_df$Core, na.rm = TRUE),
        length.out = 200
      ),
      species = unique(in_df$species)
    )
  }
  
  # Fixed-effect estimates and covariance matrix
  beta <- nlme::fixef(fit_mixed)
  V <- vcov(fit_mixed)
  
  # Population-level predictions
  predictions <- predict(
    fit_mixed,
    newdata = pred_grid,
    level = prediction_level
  )
  
  # Jacobian of predictions with respect to fixed effects
  J <- numDeriv::jacobian(
    func = function(b) {
      pred_fun(b, pred_grid, fit_mixed)
    },
    x = beta
  )
  
  # Delta-method standard errors
  pred_se <- sqrt(
    pmax(0, diag(J %*% V %*% t(J)))
  )
  
  # 95% confidence intervals
  pred_grid$fit <- predictions
  pred_grid$se <- pred_se
  pred_grid$lower <- predictions - 1.96 * pred_se
  pred_grid$upper <- predictions + 1.96 * pred_se
  
  pred_grid
}

fit_beta_model <- function(df,x_var,y_var,zi_var="0",disp_var="1") {
  
  glm_mod<-glmmTMB::glmmTMB(
    as.formula(paste(y_var,"~",x_var)),
    data = df,
    ziformula = as.formula(paste("~",zi_var)),
    dispformula = as.formula(paste("~",disp_var)),
    family = glmmTMB::beta_family(link = "logit")
  )
  
}

get_model_pred <- function(in_df, model, x_var) {
  
  # Remove missing/non-finite predictor values
  x_values <- in_df[[x_var]]
  x_values <- x_values[is.finite(x_values)]
  
  # Create prediction data across observed range
  pred_df <- data.frame(
    seq(
      min(x_values),
      max(x_values),
      length.out = 200
    )
  )
  
  names(pred_df) <- x_var
  
  # Predict fitted values and standard errors
  pred <- predict(
    model,
    newdata = pred_df,
    type = "response",
    se.fit = TRUE,
    re.form = NA
  )
  
  # Add predictions and confidence intervals
  pred_df <- pred_df %>%
    dplyr::mutate(
      fit = as.numeric(pred$fit),
      se = as.numeric(pred$se.fit),
      lower = pmax(0, fit - 1.96 * se),
      upper = pmin(1, fit + 1.96 * se)
    )
  
  return(pred_df)
}