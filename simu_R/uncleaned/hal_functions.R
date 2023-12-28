###############################################################################
#'  fit undersoomthed HAL function with global criterion
#'
#' @details The procedure select the lambda that satisfies the global 
#' undersmoothing criterion (\eqn{ P_n(\phi_{s,i}(Y-\bar{Q}_{n,\lambda}))\leq \freq{\sigma_n}{\sqrt{n}log(n)} }).
#' It performs the following steps:
#'     1). get all the basis functions from cv-hal 
#'     2). calculate a grid of new lambdas by scaling the cv-lambda from
#'          cv-hal with seq(from=0, to=-3, length=Nlam)
#'     3). refit lasso using \code{\link[glmnet]{glmnet}}) with all the basis 
#'         functions and with the grid of new lambdas
#'     4).calculate the scores (\eqn{ P_n(\phi_{s,i}(Y-\bar{Q}_{n,\lambda})) })of each lasso fit     
#'     5).find the biggest lambda that max score (\eqn{ \leq \freq{\sigma_n}{\sqrt{n}log(n)} })
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param Y A \code{numeric} vector of observations of the outcome variable.
#' @param fit_init The initial HAL fit object from the output list of \code{undersmooth_init}.
#' @param Nlam Number of lambda candidates. The sequence ranges from \code{fit_init$lambda_star} to
#' \code{fit_init$lambda_star*10^(-3)}.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.

undersmooth_hal <- function(X,
                            Y,
                            fit_init,
                            Nlam = 20,
                            family = c("gaussian", "binomial", "poisson", "cox")){
  
  n = length(Y)
  
  nonzero_col_init = which(fit_init$coefs[-1] != 0)
  if (length(nonzero_col_init) == 0){
    res <- list("lambda_init" = fit_init$lambda_star,
                "lambda_under" = fit_init$lambda_star)
    return(res)
  }
  
  #  refit on a grid of lambda (scaled on the lambda from cv-hal)
  us_lambda <- fit_init$lambda_star*10^seq(from=0, to=-3, length=Nlam)
  us_fit <- glmnet(fit_init$x_basis, Y, lambda=us_lambda, family = family, standardize = FALSE)
  
  if(identical(us_fit$df, 0)){
    res <- list("lambda_init" = fit_init$lambda_star,
                "lambda_under" = fit_init$lambda_star)
    return(res)
  }
  
  
  preds_init <- predict(fit_init, new_data = X)
  resid_init <- preds_init - Y
  
  if (family != "binomial"){
    pred_mat <- predict(us_fit, fit_init$x_basis)
  }else {
    pred_mat <- predict(us_fit, fit_init$x_basis, type = "response")
  }
  resid_mat <- pred_mat - Y
  

  ## estimates of sd in each direction using initial fit
  basis_mat_init <- as.matrix(fit_init$x_basis)
  basis_mat_init <- as.matrix(basis_mat_init[, nonzero_col_init])
  
  sd_est  <- apply(basis_mat_init, 2, function(phi) sd(resid_init*phi))
  
  
  ## calculate scores
  max_score <- get_maxscore(basis_mat = basis_mat_init,
                            resid_mat = resid_mat,
                            sd_est = sd_est,
                            Nlam = Nlam, us_fit = us_fit)
  
  ## get the first lambda that satisfies the criteria
  lambda_under <- us_lambda[max_score <= 1/(sqrt(n)*log(n))][1] # over under-smoothing 
  
  if (is.na(lambda_under)){
    res <- list("lambda_init" = fit_init$lambda_star,
                "lambda_under" = fit_init$lambda_star)
    return(res)
  }
  
  # collect results
  coef_mat <- as.matrix(us_fit$beta)
  
  spec_under <- list("lambda" = us_lambda,
                     "l1_norm" = NA,
                     "n_coef" = NA)
  
  spec_under$l1_norm <- apply(coef_mat, 2, function(x){sum(abs(x))})
  spec_under$n_coef <- apply(coef_mat, 2, function(x){sum(x != 0)})
  
  res <- list("lambda_init" = fit_init$lambda_star,
              "lambda_under" = lambda_under,
              "spec_under" = spec_under)
  return(res)
}


###############################################################################
#'  undersoomthed HAL helper function for global criterion
#'
#' @details For each candidate lambda, do:
#'     1). standardize the score formed by each basis.
#'     2). calculate the mean of the standardized scores for each basis.
#' Select the max of the mean.
#' @param basis_mat The selected basis matrix from initial fit for undersmoothing,
#'  obtained from the output list of \code{undersmooth_init}.
#' @param resid_mat The residual matrix with each column the residuals correspongding to a lambda.
#' @param sd_est A numeric vector containing the sd of each column of \code{basis_mat}.
#' @param Nlam Number of lambda candidates.
#' @param us_fit The \code{glmnet} fit of the sequence of candidate lambdas.

get_maxscore <- function(basis_mat, resid_mat, sd_est, Nlam, us_fit){
  
  basis_mat_sd <- sweep(basis_mat, 2, sd_est, FUN = '/')
  score_all <- apply(basis_mat_sd, 2, function(u) {
    score_mat <- resid_mat * u
    score_mean <- apply(score_mat, 2, mean)
  })
  # absolute value
  max_score <- apply(abs(score_all), 1, max, na.rm=T)
  return(max_score)
}



###############################################################################
#'  With given fitted HAL object and evaluation points, return the empirical SE

IC_based_se <- function(X, Y, hal_fit, eval_points, family = "binomial", X_unpenalized = NULL ){
  
  n = length(Y)
  coef <- hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
  
  nonzero_idx <- which(coef != 0)
  
  if(length(nonzero_idx) > 0) {
    coef_nonzero <- coef[nonzero_idx]
    basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
    
    if(family == "binomial"){
      if(is.null(X_unpenalized)){
        Y_hat = predict(hal_fit, new_data = X, type = "response")
      } else {
        Y_hat = predict(hal_fit, new_data = X, new_X_unpenalized = as.matrix(X_unpenalized), type = "response")
      }
      

    } else {
      if(is.null(X_unpenalized)){
        Y_hat = predict(hal_fit, new_data = X)
      } else {
        Y_hat = predict(hal_fit, new_X_unpenalized = as.matrix(X_unpenalized), new_data = X)
      }
      
    }
    
    IC_beta <- cal_IC_for_beta(X = basis_mat_nonzero, 
                               Y = Y, 
                               Y_hat =Y_hat,
                               beta_n = coef_nonzero
    )
    
    if(any(! is.na(IC_beta))){
      se <- c()
      
      for (i in 1:length(eval_points)) {
        X_new <- X
        X_new$A = eval_points[i]
        
        # efficient influence curve
        x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
        if(is.null(X_unpenalized)){
          x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
        } else {
          x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a, as.matrix(X_unpenalized))[, nonzero_idx])
        }
        
        
        IC_EY <- cal_IC_for_EY(X_new = x_basis_a_nonzero, 
                               beta_n = coef_nonzero, IC_beta = IC_beta)
        
        # empirical SE
        se[i] <- sqrt(var(IC_EY)/n)
        
      }
    } else {
      se <- NA
    }
    
  } else {
    se <- NA
  }
  
  return(se)
}

###############################################################################

IC_based_se_u_l <- function(X, Y, hal_fit, single_eval_point){
  
  coef <- hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
  
  nonzero_idx <- which(coef != 0)
  
  if(length(nonzero_idx) > 0) {
    coef_nonzero <- coef[nonzero_idx]
    basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
    
    IC_beta <- cal_IC_for_beta(X = basis_mat_nonzero, 
                               Y = Y, 
                               Y_hat = predict(hal_fit, new_data = X, type = "response"),
                               beta_n = coef_nonzero
    )
    
    X_new <- X
    X_new$A = single_eval_point
    
    # efficient influence curve
    x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
    x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
    
    IC_EY <- cal_IC_for_EY(X_new = x_basis_a_nonzero, 
                           beta_n = coef_nonzero, IC_beta = IC_beta)
    
    # empirical SE
    se <- sqrt(var(IC_EY)/n)
    
  } else {
    se <- NA
  }
  
  return(se)
}



###############################################################################
# calculating efficient influence curves

cal_IC_for_beta <- function(X, Y, Y_hat, beta_n, family = 'binomial', X_unpenalized = NULL){
  n <- dim(X)[1] 
  p <- length(beta_n)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # 1. calculate score: X'(Y - phi(X))
  res <- Y-Y_hat
  score <- sweep(t(X), 2, res, `*`)
  
  # 2. calculate the derivative of phi:
  if(family == 'binomial'){
    d_phi_scaler <- as.vector(exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2)) # exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2))
    d_phi <- sweep(X, 1, d_phi_scaler, `*`)
  } else {
    d_phi = - X
  }
  
  # 3. -E_{P_n}(X d_phi)^(-1)
  tmat <- t(X) %*% d_phi / n
  if(! is.matrix(try(solve(tmat), silent = TRUE))){
    return(NA)
  }
  tmat <- -solve(tmat)
  
  # 4. calculate influence curves
  IC <- tmat %*% score
  
  return(IC)
}


cal_IC_for_EY <- function(X_new, beta_n, IC_beta, family = 'binomial'){
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  if(family == 'binomial'){
    d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X_new)) / ((1 + exp(- beta_n %*% t(X_new)))^2))
    d_phi_new <- sweep(X_new, 1, d_phi_scaler_new, `*`)
  } else {
    d_phi_new = X_new
  }
  
  IC = diag(d_phi_new %*% IC_beta)
  
  return(IC)
}


cal_IC_for_ATE <- function(X_new_a, X_new_0, beta_n, IC_beta, family = 'binomial'){
  
  if (!is.matrix(X_new_a)) X_new_a <- as.matrix(X_new_a)
  if (!is.matrix(X_new_0)) X_new_0 <- as.matrix(X_new_0)
  
  if (family == 'binomial') {
    d_phi_scaler_new_a <- as.vector(exp(- beta_n %*% t(X_new_a)) / ((1 + exp(- beta_n %*% t(X_new_a)))^2))
    d_phi_new_a <- sweep(X_new_a, 1, d_phi_scaler_new_a, `*`)
    
    d_phi_scaler_new_0 <- as.vector(exp(- beta_n %*% t(X_new_0)) / ((1 + exp(- beta_n %*% t(X_new_0)))^2))
    d_phi_new_0 <- sweep(X_new_0, 1, d_phi_scaler_new_0, `*`)
    
    d_phi_new <- d_phi_new_a - d_phi_new_0
  } else {
    d_phi_new <- X_new_a - X_new_0
  }
  
  IC = diag(d_phi_new %*% IC_beta)
  
  return(IC)
}


###############################################################################

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
  #   
  #   basis_list <- hal_fit$basis_list
  #   # copy_map <- hal_fit$copy_map
  #   
  #   y_hat_B <- matrix(ncol = length(eval_points))
  #   for (b in 1:B) {
  #     #--------------data--------------
  #     idx <- sample(1:length(Y), length(Y), replace = T)
  #     Xb <- X[idx,]
  #     Yb <- Y[idx]
  #     
  #     # generate basis matrix
  #     if (length(basis_list) > 0) {
  #       x_basis <- hal9001::make_design_matrix(as.matrix(Xb), basis_list)
  #       # unique_columns <- as.numeric(names(copy_map))
  #       # x_basis <- x_basis[, unique_columns]
  #     } else {
  #       x_basis <- matrix(1, ncol = 2, nrow = nrow(Xb))
  #     }
  #     x_basis <- as.matrix(x_basis)
  #     
  #     #--------------fit--------------
  #     is_glmnet = T
  #     if (dim(x_basis)[2] <= 1) {
  #       # dim of X_basis < 2. make it larger
  #       x_basis <- cbind(matrix(1, ncol = 1, nrow = nrow(X)), x_basis)
  #       x_basis <- cbind(matrix(0, ncol = 1, nrow = nrow(X)), x_basis)
  #       lasso_fit <- glm(Yb ~ x_basis, x = FALSE, y = FALSE, family = y_type)
  #       is_glmnet = F
  #     } else {
  #       lasso_fit <- tryCatch({
  #         lasso_fit <- glmnet::glmnet(x = x_basis, y = Yb, 
  #                                     family = y_type, 
  #                                     lambda = hal_fit$lambda_star,
  #                                     intercept = FALSE, standardize = FALSE)
  #       },
  #       error = function(){
  #         lasso_fit <- glm.fit(x = x_basis, y = Yb, family = y_type)
  #         lasso_fit <- glm(Yb ~ x_basis, x = FALSE, y = FALSE, family = y_type)
  #         is_glmnet = F
  #       })
  #     }
  #     
  #     #--------------predictions, 95% lower and upper bounds--------------
  #     y_hat_b <- c()
  #     for (i in 1:length(eval_points)) {
  #       X_new <- Xb
  #       X_new$A = eval_points[i]
  #       
  #       # generate basis matrix
  #       if (length(basis_list) > 0){
  #         x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit$basis_list)
  #         # x_basis_a <- hal9001::apply_copy_map(x_basis_a, hal_fit$copy_map)
  #         
  #         if(dim(x_basis)[2] <= 1){
  #           x_basis_a <- cbind(matrix(1, ncol = 1, nrow = nrow(X)), x_basis_a)
  #           x_basis_a <- cbind(matrix(0, ncol = 1, nrow = nrow(X)), x_basis_a)
  #         }
  #       } else {
  #         x_basis_a <- matrix(1, ncol = 2, nrow = nrow(X_new))
  #       }
  # 
  #       # prediction
  #       #-----
  #       if (y_type == "binomial") {
  #         preds <- predict(lasso_fit, x_basis_a, type = "response")
  #       } else {
  #         preds <- predict(lasso_fit, x_basis_a)
  #       }
  #       #-----
  #       # beta_hat <- stats::coef(lasso_fit)
  #       # beta_hat[is.na(beta_hat)] <- 0
  #       # beta_hat <- as.matrix(beta_hat)
  #       # preds <- as.vector(
  #       #   Matrix::tcrossprod(x = x_basis_a, y = beta_hat[-1]) + beta_hat[1]
  #       # )
  #       # if (y_type == "binomial") preds <- stats::plogis(preds)
  #       #-----
  #       y_hat_b[i] = mean(preds)
  #     }
  #     y_hat_B <- rbind(y_hat_B, y_hat_b)
  #   }
  #   
  # y_hat_B <- y_hat_B[-1,]
  # lower_bd <- apply(y_hat_B, 2, quantile, probs = 0.05/2)
  # upper_bd <- apply(y_hat_B, 2, quantile, probs = 1-0.05/2)
  # 
  # return(list(lower_bd=lower_bd, upper_bd=upper_bd))
  
}


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
  #   
  #   y_hat_B <- matrix(ncol = length(eval_points))
  #   for (b in 1:B) {
  #     #--------------data--------------
  #     idx <- sample(1:length(Y), length(Y), replace = T)
  #     Xb <- X[idx,]
  #     Yb <- Y[idx]
  #     
  #     
  #     x_basis <- list()
  #     is_glmnet <- list()
  #     lasso_fit <- list()
  #     
  #     for (i in 1:length(hal_fit)) {
  #       
  #       # generate basis matrix
  #       basis_list <- hal_fit[[i]]$basis_list
  #       if (length(basis_list) > 0) {
  #         x_basis_l <- hal9001::make_design_matrix(as.matrix(Xb), basis_list)
  #       } else {
  #         x_basis_l <- matrix(1, ncol = 2, nrow = nrow(Xb))
  #       }
  #       x_basis[[i]] <- as.matrix(x_basis_l)
  #       
  #       #--------------fit--------------
  #       is_glmnet_l = T
  #       
  #       if (dim(x_basis_l)[2] <= 1) {
  #         # dim of X_basis < 2. make it larger
  #         x_basis_l <- cbind(matrix(1, ncol = 1, nrow = nrow(X)), x_basis_l)
  #         x_basis_l <- cbind(matrix(0, ncol = 1, nrow = nrow(X)), x_basis_l)
  #         lasso_fit_l <- glm(Yb ~ x_basis_l, x = FALSE, y = FALSE, family = y_type)
  #         is_glmnet_l = F
  #       } else {
  #         lasso_fit_l <- tryCatch({
  #           lasso_fit_l <- glmnet::glmnet(x = x_basis_l, y = Yb, 
  #                                       family = y_type, 
  #                                       lambda = hal_fit[[i]]$lambda_star,
  #                                       intercept = FALSE, standardize = FALSE)
  #         },
  #         error = function(){
  #           lasso_fit_l <- glm.fit(x = x_basis_l, y = Yb, family = y_type)
  #           lasso_fit_l <- glm(Yb ~ x_basis_l, x = FALSE, y = FALSE, family = y_type)
  #           is_glmnet_l = F
  #         })
  #       }
  #       is_glmnet[[i]] = is_glmnet_l
  #       lasso_fit[[i]] = lasso_fit_l
  #       
  #     }
  #     
  #     
  #     #--------------predictions, 95% lower and upper bounds--------------
  #     y_hat_b <- c()
  #     for (i in 1:length(eval_points)) {
  #       X_new <- Xb
  #       X_new$A = eval_points[i]
  #       u_l_idx = lambda_u_l_idx[i]
  #       
  #       # generate basis matrix
  #       if (length(basis_list) > 0){
  #         x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit[[u_l_idx]]$basis_list)
  # 
  #         if(dim(x_basis[[u_l_idx]])[2] <= 1){
  #           x_basis_a <- cbind(matrix(1, ncol = 1, nrow = nrow(X)), x_basis_a)
  #           x_basis_a <- cbind(matrix(0, ncol = 1, nrow = nrow(X)), x_basis_a)
  #         }
  #       } else {
  #         x_basis_a <- matrix(1, ncol = 2, nrow = nrow(X_new))
  #       }
  #       
  #       # prediction
  #       if (y_type == "binomial") {
  #         preds <- predict(lasso_fit[[u_l_idx]], x_basis_a, type = "response")
  #       } else {
  #         preds <- predict(lasso_fit[[u_l_idx]], x_basis_a)
  #       }
  # 
  #       
  #       y_hat_b[i] = mean(preds)
  #     }
  #     y_hat_B <- rbind(y_hat_B, y_hat_b)
  #   }
  #   
  #   y_hat_B <- y_hat_B[-1,]
  #   lower_bd <- apply(y_hat_B, 2, quantile, probs = 0.05/2)
  #   upper_bd <- apply(y_hat_B, 2, quantile, probs = 1-0.05/2)
  #   
  #   return(list(lower_bd=lower_bd, upper_bd=upper_bd))
  #   
}


fit_undersmoothed_smoothness_adaptive_HAL <- function(df, y_name, x_names, family){
  
  Y <- as.numeric(df[, y_name])
  
  X <- df %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  # 1. --- SL
  task <- make_sl3_Task(
    data = df,
    outcome = y_name,
    covariates = x_names
  )
  
  # num_knots = c(200, 100,  50)
  lrn_hal0 <- Lrnr_hal9001$new(smoothness_orders = 0, family = family, return_x_basis = TRUE)
  lrn_hal1 <- Lrnr_hal9001$new(smoothness_orders = 1, family = family, return_x_basis = TRUE)
  lrn_hal2 <- Lrnr_hal9001$new(smoothness_orders = 2, family = family, return_x_basis = TRUE)
  lrn_hal3 <- Lrnr_hal9001$new(smoothness_orders = 3, family = family, return_x_basis = TRUE)
  
  # num_knots = c(20,10,5)
  lrn_hal0_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 0,
                                            family = family, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 0,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal1_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 1,
                                            family = family, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 1,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal2_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 2,
                                            family = family, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 2,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal3_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 3,
                                            family = family, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 3,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  
  learners <- c(lrn_hal0, lrn_hal0_smallerKnots, lrn_hal1, lrn_hal1_smallerKnots, 
                lrn_hal2, lrn_hal2_smallerKnots, lrn_hal3, lrn_hal3_smallerKnots)
  names(learners) <- c("halcv0", "halcv0_sKnots", "halcv1", "halcv1_sKnots", 
                       "halcv2", "halcv2_sKnots", "halcv3", "halcv3_sKnots")
  stack <- make_learner(Stack, learners)
  
  if(family == "binomial"){
    cv_selector <- Lrnr_cv_selector$new(eval_function = loss_loglik_binomial) # https://tlverse.org/sl3/reference/loss_functions.html
  } else {
    cv_selector <- Lrnr_cv_selector$new(eval_function = loss_squared_error) # https://tlverse.org/sl3/reference/loss_functions.html
  }
  
  dSL <- Lrnr_sl$new(learners = stack, metalearner = cv_selector)
  
  dSL_fit <- dSL$train(task)
  
  #--- 2. SL pick undersmooth
  smooth_orders = c(0,0,1,1,2,2,3,3)
  num_knots = rep(c("default", "smaller"),4)
  
  # undersmo0th hal based on the one sl picked
  sl_pick_idx = which(dSL_fit$coefficients==1)
  
  hal_fit_sl_pick = dSL_fit$learner_fits[[sl_pick_idx]]$fit_object
  SO_pick = smooth_orders[sl_pick_idx]
  n_knots_default_pick = as.numeric(num_knots[sl_pick_idx] == "default")
  
  hal_undersmooth <- undersmooth_hal(X, Y,
                                     fit_init = hal_fit_sl_pick,
                                     family = family)
  
  lambda_u_g = hal_undersmooth$lambda_under
  
  if(is.na(lambda_u_g) | hal_undersmooth$lambda_under == hal_undersmooth$lambda_init){
    lambda_u_g <- lambda_CV
    print(sprintf('  globally u lambdas: %f', lambda_u_g))
    
    hal_u_g <- hal_fit_sl_pick
  } else {
    if(n_knots_default_pick == 1){
      hal_u_g <- fit_hal(X = X, Y = Y, 
                         family = family,
                         smoothness_orders = SO_pick,
                         return_x_basis = TRUE,
                         fit_control = list(cv_select = FALSE),
                         lambda = lambda_u_g)
    } else {
      hal_u_g <- fit_hal(X = X, Y = Y, 
                         family = family,
                         smoothness_orders = SO_pick,
                         return_x_basis = TRUE,
                         fit_control = list(cv_select = FALSE),
                         lambda = lambda_u_g,
                         num_knots = hal9001:::num_knots_generator(
                           max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                           smoothness_orders = SO_pick,
                           base_num_knots_0 = 20,
                           base_num_knots_1 = 20  
                         ))
    }
    
  }
  
  num_basis_u_g <- sum(hal_u_g$coefs[-1] != 0)
  
  returns = list(hal_fit_object = hal_u_g,
                 hal_fit_info = list(smooth_order = SO_pick,
                                     n_knots_default = n_knots_default_pick,
                                     X = X,
                                     Y = Y,
                                     lambda = lambda_u_g,
                                     lambda_scaler = lambda_u_g / hal_undersmooth$lambda_init, 
                                     num_basis = num_basis_u_g))
  
  return(returns) 
}


predict_curve = function(X, Y, points_curve, hal_fit, family = "binomial"){
 
  # if (!is.matrix(X)) X <- as.matrix(X)
  psi_hat_u_g <- sapply(points_curve, function(a){ X_new <- X
  X_new[,"A"] = a
  if(hal_fit$unpenalized_covariates == 0){
    mean(predict(hal_fit, new_data = X_new, type = "response"))
  } else {
    mean(predict(hal_fit, new_data = X_new, new_X_unpenalized = as.matrix(X_new), type = "response"))
  }
  } )
  
  # IC-based inference
  if(hal_fit$unpenalized_covariates == 0){
    psi_hat_pnt_u_g_se <- IC_based_se(X, Y, hal_fit, points_curve)
  } else {
    psi_hat_pnt_u_g_se <- IC_based_se(X, Y, hal_fit, points_curve, X_unpenalized = X)
  }
  
  
  psi_hat_df <- cbind(points_curve, matrix(psi_hat_u_g, ncol=1), psi_hat_pnt_u_g_se)
  
  colnames(psi_hat_df) <- c("a", "y_hat", "SE")
  
  psi_hat_df <- as.data.frame(psi_hat_df) %>% 
    mutate(ci_lwr = y_hat - 1.96 * SE,
           ci_upr = y_hat + 1.96 * SE)
  
  if(family == "binomial"){
    bounds <- c(0, 1)
  } else {
    bounds <- c(min(Y), max(Y))
  }
  psi_hat_df[,"y_hat"] <- pmax(bounds[1], psi_hat_df[,"y_hat"])
  psi_hat_df[,"y_hat"] <- pmin(psi_hat_df[,"y_hat"], bounds[2])
  
  psi_hat_df[,"ci_lwr"] <- pmax(bounds[1], psi_hat_df[,"ci_lwr"])
  psi_hat_df[,"ci_lwr"] <- pmin(psi_hat_df[,"ci_lwr"], bounds[2])
  
  psi_hat_df[,"ci_upr"] = pmax(bounds[1], psi_hat_df[,"ci_upr"])
  psi_hat_df[,"ci_upr"] <- pmin(psi_hat_df[,"ci_upr"], bounds[2])
  
  return(psi_hat_df)
    
}



summary_hal9001_unpenalized <- function(object,
                            lambda = NULL,
                            only_nonzero_coefs = TRUE,
                            include_redundant_terms = FALSE,
                            round_cutoffs = 3,
                            ...) {
  abs_coef <- basis_list_idx <- coef_idx <- dup <- NULL
  
  # retain coefficients corresponding to lambda
  if (!is.null(lambda)) {
    if (length(lambda) > 1) {
      stop("Cannot summarize over multiple values of lambda.")
    }
    if (lambda != object$lambda_star) {
      if (is.null(object$lasso_fit)) {
        stop(
          "Coefficients for specified lambda do not exist, or are not ",
          "accessible since the fit of the lasso model was not returned ",
          "(i.e., return_lasso was set to FALSE in `hal_fit()`)."
        )
      } else {
        if (!(lambda %in% object$lasso_fit$lambda)) {
          stop("Coefficients for the specified lambda do not exist.")
        } else {
          lambda_idx <- which(object$lasso_fit$lambda == lambda)
          coefs <- object$lasso_fit$glmnet.fit$beta[, lambda_idx]
        }
      }
    } else {
      lambda_idx <- which(object$lambda_star == lambda)
      coefs <- object$coefs[, lambda_idx]
    }
  }
  
  if (is.null(lambda)) {
    lambda <- object$lambda_star
    coefs <- object$coefs
    if (length(lambda) > 1) {
      warning(
        "Coefficients for many lambda exist --\n",
        "Summarizing coefficients corresponding to minimum lambda."
      )
      lambda_idx <- which.min(lambda)
      coefs <- object$coefs[, lambda_idx]
    }
  }
  
  # cox model has no intercept
  if (object$family != "cox") {
    coefs_no_intercept <- coefs[-1]
  } else {
    coefs_no_intercept <- coefs
  }
  
  # subset to non-zero coefficients
  if (only_nonzero_coefs) {
    coef_idxs <- which(coefs_no_intercept != 0)
  } else {
    coef_idxs <- seq_along(coefs_no_intercept)
  }
  copy_map <- object$copy_map[coef_idxs]
  
  if (object$unpenalized_covariates > 0){
    copy_map <- copy_map[1:(length(copy_map) - object$unpenalized_covariates)]
  }
  
  # summarize coefficients with respect to basis list
  coefs_summ <- data.table::rbindlist(
    lapply(seq_along(copy_map), function(map_idx) {
      coef_idx <- coef_idxs[map_idx]
      coef <- coefs_no_intercept[coef_idx]
      
      basis_list_idxs <- copy_map[[map_idx]] # indices of duplicates
      basis_dups <- object$basis_list[basis_list_idxs]
      
      data.table::rbindlist(
        lapply(seq_along(basis_dups), function(i) {
          coef_idx <- ifelse(object$family != "cox", coef_idx + 1, coef_idx)
          dt <- data.table::data.table(
            coef_idx = coef_idx, # coefficient index
            coef, # coefficient
            basis_list_idx = basis_list_idxs[i], # basis list index
            col_idx = basis_dups[[i]]$cols, # column idx in X
            col_cutoff = basis_dups[[i]]$cutoffs, # cutoff
            col_order = basis_dups[[i]]$orders # smoothness order
          )
          return(dt)
        })
      )
    })
  )
  
  if (!include_redundant_terms) {
    coef_idxs <- unique(coefs_summ$coef_idx)
    coefs_summ <- data.table::rbindlist(lapply(coef_idxs, function(idx) {
      # subset to matching coefficient index
      coef_summ <- coefs_summ[coef_idx == idx]
      
      # label duplicates (i.e. basis functions with identical col & cutoff)
      dups_tbl <- coef_summ[, c("col_idx", "col_cutoff", "col_order")]
      if (!anyDuplicated(dups_tbl)) {
        return(coef_summ)
      } else {
        # add col indicating whether or not there is a duplicate
        coef_summ[, dup := (duplicated(dups_tbl) |
                              duplicated(dups_tbl, fromLast = TRUE))]
        
        # if basis_list_idx contains redundant duplicates, remove them
        redundant_dups <- coef_summ[dup == TRUE, "basis_list_idx"]
        if (nrow(redundant_dups) > 1) {
          # keep the redundant duplicate term that has the shortest length
          retain_idx <- which.min(apply(redundant_dups, 1, function(idx) {
            nrow(coef_summ[basis_list_idx == idx])
          }))
          idx_keep <- unname(unlist(redundant_dups[retain_idx]))
          coef_summ <- coef_summ[basis_list_idx == idx_keep]
        }
        return(coef_summ[, -"dup"])
      }
    }))
  }
  
  # summarize with respect to x column names:
  x_names <- data.table::data.table(
    col_idx = 1:length(object$X_colnames),
    col_names = object$X_colnames
  )
  summ <- merge(coefs_summ, x_names, by = "col_idx", all.x = TRUE)
  
  # combine name, cutoff into 0-order basis function (may include interaction)
  summ$zero_term <- paste0(
    "I(", summ$col_names, " >= ", round(summ$col_cutoff, round_cutoffs), ")"
  )
  summ$higher_term <- ifelse(
    summ$col_order == 0, "",
    paste0(
      "(", summ$col_names, " - ",
      round(summ$col_cutoff, round_cutoffs), ")"
    )
  )
  summ$higher_term <- ifelse(
    summ$col_order < 1, summ$higher_term,
    paste0(summ$higher_term, "^", summ$col_order)
  )
  summ$term <- ifelse(
    summ$col_order == 0,
    paste0("[ ", summ$zero_term, " ]"),
    paste0("[ ", summ$zero_term, "*", summ$higher_term, " ]")
  )
  
  term_tbl <- data.table::as.data.table(stats::aggregate(
    term ~ basis_list_idx,
    data = summ, paste, collapse = " * "
  ))
  
  # no longer need the columns or rows that were incorporated in the term
  redundant <- c(
    "term", "col_cutoff", "col_names", "col_idx", "col_order", "zero_term",
    "higher_term"
  )
  summ <- summ[, -..redundant]
  summ_unique <- unique(summ)
  summ <- merge(
    term_tbl, summ_unique,
    by = "basis_list_idx", all.x = TRUE, all.y = FALSE
  )
  
  # summarize in a list
  coefs_list <- lapply(unique(summ$coef_idx), function(this_coef_idx) {
    coef_terms <- summ[coef_idx == this_coef_idx]
    list(coef = unique(coef_terms$coef), term = t(coef_terms$term))
  })
  
  # summarize in a table
  coefs_tbl <- data.table::as.data.table(stats::aggregate(
    term ~ coef_idx,
    data = summ, FUN = paste, collapse = "  OR  "
  ))
  redundant <- c("term", "basis_list_idx")
  summ_unique_coefs <- unique(summ[, -..redundant])
  coefs_tbl <- data.table::data.table(merge(
    summ_unique_coefs, coefs_tbl,
    by = "coef_idx", all = TRUE
  ))
  coefs_tbl[, "abs_coef" := abs(coef)]
  coefs_tbl <- data.table::setorder(coefs_tbl[, -"coef_idx"], -abs_coef)
  coefs_tbl <- coefs_tbl[, -"abs_coef", with = FALSE]
  
  # incorporate intercept
  if (object$family != "cox") {
    intercept <- list(data.table::data.table(
      coef = coefs[1], term = "(Intercept)"
    ))
    
    
    if (object$unpenalized_covariates > 0){
      coefs_unpenalized_tbl <- list(data.table::data.table(
        coef = object$coefs[(nrow(object$coefs) - object$unpenalized_covariates + 1) : nrow(object$coefs),], 
        term = names(object$coefs[(nrow(object$coefs) - object$unpenalized_covariates + 1) : nrow(object$coefs),])
      ))
      
      coefs_tbl <- data.table::rbindlist(
        c(intercept, list(coefs_tbl), c(coefs_unpenalized_tbl)),
        fill = TRUE
      )
      
      intercept <- list(coef = coefs[1], term = "(Intercept)")
      coefs_list <- c(list(intercept), list(coefs_unpenalized_tbl), coefs_list)
    } else {
      coefs_tbl <- data.table::rbindlist(
        c(intercept, list(coefs_tbl)),
        fill = TRUE
      )
      intercept <- list(coef = coefs[1], term = "(Intercept)")
      coefs_list <- c(list(intercept), coefs_list)
    }
    
  }
  
  out <- list(
    table = coefs_tbl,
    list = coefs_list,
    lambda = lambda,
    only_nonzero_coefs = only_nonzero_coefs
  )
  class(out) <- "summary.hal9001"
  return(out)
}

