
bootstrap_inference <- function(X, Y, eval_points, hal_fit, y_type, B = 200){
  
  basis_list <- hal_fit$basis_list
  
  y_hat_B <- matrix(ncol = length(eval_points))
  
  for (b in 1:B) {
    #--------------data--------------
    idx <- sample(1:length(Y), length(Y), replace = T)
    Xb <- X[idx,]
    Yb <- Y[idx]
    
    if (length(basis_list) > 0) {
      # generate basis matrix
      x_basis <- hal9001::make_design_matrix(as.matrix(Xb), basis_list)
      
      # maybe
      
      #--------------fit--------------
      lasso_fit <- tryCatch({
        lasso_fit <- glmnet::glmnet(x = x_basis, y = Yb, 
                                    family = y_type, 
                                    lambda = hal_fit$lambda_star,
                                    intercept = FALSE, standardize = FALSE)
      },
      error = function(){
        lasso_fit <- NA
      })
      
      y_hat_b <- c()
      for (i in 1:length(eval_points)) {
        X_new <- Xb
        X_new$A = eval_points[i]
        
        x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit$basis_list)
        
        #--------------prediction--------------
        if (y_type == "binomial") {
          preds <- predict(lasso_fit, x_basis_a, type = "response")
        } else {
          preds <- predict(lasso_fit, x_basis_a)
        }
        
        y_hat_b[i] = mean(preds)
      }
    } else {
      y_hat_b = rep(NA, length(eval_points))
    }
    
    y_hat_B <- rbind(y_hat_B, y_hat_b)
  }
  
  #--------------confidence bounds--------------
  y_hat_B <- y_hat_B[-1,]
  lower_bd <- apply(y_hat_B, 2, quantile, probs = 0.05/2, na.rm = T)
  upper_bd <- apply(y_hat_B, 2, quantile, probs = 1-0.05/2, na.rm = T)
  SE <- sqrt(apply(y_hat_B, 2, var, na.rm = T))
  
  return(list(lower_bd=lower_bd, upper_bd=upper_bd, SE=SE))
  
###############################################################################

bootstrap_inference_u_l <- function(X, Y, eval_points, hal_fit, y_type, lambda_u_l_idx, B = 200){
  
  y_hat_B <- matrix(ncol = length(eval_points))
  
  for (b in 1:B) {
    #--------------data--------------
    idx <- sample(1:length(Y), length(Y), replace = T)
    Xb <- X[idx,]
    Yb <- Y[idx]
    
    
    x_basis <- list()
    lasso_fit <- list()
    
    for (i in 1:length(hal_fit)) {
      
      # generate basis matrix
      basis_list <- hal_fit[[i]]$basis_list
      
      # --------------fit glmnet--------------
      if (length(basis_list) > 0) {
        
        x_basis_l <- hal9001::make_design_matrix(as.matrix(Xb), basis_list)
        
        lasso_fit_l <- tryCatch({
          lasso_fit_l <- glmnet::glmnet(x = x_basis_l, 
                                        y = Yb, 
                                        family = y_type, 
                                        lambda = hal_fit[[i]]$lambda_star,
                                        intercept = FALSE, standardize = FALSE)
        },
        error = function(){
          lasso_fit_l <- NA
        })
        
      } else {
        
        x_basis_l <- NA
        lasso_fit_l <- NA
        
      } 
      
      x_basis[[i]] = x_basis_l
      lasso_fit[[i]] = lasso_fit_l
      
    }
    
    #--------------predictions, 95% lower and upper bounds--------------
    y_hat_b <- c()
    for (i in 1:length(eval_points)) {
      X_new <- Xb
      X_new$A = eval_points[i]
      u_l_idx = lambda_u_l_idx[i]
      
      if (length(basis_list) > 0){
        
        # generate basis matrix
        x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit[[u_l_idx]]$basis_list)
        
        # prediction
        if(any(!is.na(lasso_fit[[u_l_idx]]))){
          
          if (y_type == "binomial") {
            preds <- predict(lasso_fit[[u_l_idx]], x_basis_a, type = "response")
          } else {
            preds <- predict(lasso_fit[[u_l_idx]], x_basis_a)
          }
          
          y_hat_b[i] = mean(preds)
          
        } else {
          y_hat_b[i] = NA
        }
        
      } else {
        y_hat_b[i] = NA
      }
      
    }
    y_hat_B <- rbind(y_hat_B, y_hat_b)
  }
  
  y_hat_B <- y_hat_B[-1,]
  lower_bd <- apply(y_hat_B, 2, quantile, probs = 0.05/2, na.rm = T)
  upper_bd <- apply(y_hat_B, 2, quantile, probs = 1-0.05/2, na.rm = T)
  SE <- sqrt(apply(y_hat_B, 2, var, na.rm = T))
  
  return(list(lower_bd=lower_bd, upper_bd=upper_bd, SE=SE))
}


