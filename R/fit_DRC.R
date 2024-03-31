
###############################################################################
#'  fitting HAL and learn the dose response curve
#'
#' @details The procedure fits undersmoothed HAL with provided data and HAL fitting 
#' hyperparameters, then it learns predicts dose response curve using the fitted 
#' HAL as a working model and provide delta-method based CI.
#' @param data An input \code{data.frame} with dimensions number of observations -by-
#'  number of covariates plus the outcome variable. 
#' @param y_var_name A \code{character} of the outcome variable name
#' @param trt_var_name A \code{character} of the treatment variable name
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.
#' @param curvePoints A \code{numeric vector} with points on the dose response 
#' curve that wants to be estimated.
#' @param smoothOrderAdapt A \code{boolean} indicating whether data-adaptivly
#' selects the smoothness order and base number knots for HAL fitting. 
#' @param smoothOrder A \code{numeric} smoothness order number of HAL fit. 
#' @param baseNumKnots A \code{numeric} base_num_knots of HAL fit. 
#' @param boundResults A \code{boolean} indicating whether bound the predicted 
#' values and confidence intervals by the minimum and maximum of the provided Y.

fit_UHAL_DRC <- function(dat, y_var_name, trt_var_name, family, curvePoints = NA, smoothOrderAdapt = FALSE, smoothOrder = 1,  baseNumKnots = 20, boundResults = TRUE){
  
  n <- length(Y)

  Y <- as.numeric(as.matrix(dat %>% select(all_of(y_var_name))))
  
  
  x_names = names(dat)[names(dat) != y_var_name]
  
  X <- dat %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  names(X)[names(X) == trt_var_name] = "A"

  if(is.na(curvePoints)){
    curvePoints =  seq(min(X$A), max(X$A), by = 0.1)
  }
  
  # 1. fit undersmoothed HAL
  ## fit CV HAL
  if(smoothOrderAdapt){
    dSL_fit <- fit_SL_smoothness_adaptive_HAL(dat, X, Y, x_names, y_name, family)

    sl_pick_idx = which(dSL_fit$coefficients==1)
    
    hal_CV = dSL_fit$learner_fits[[sl_pick_idx]]$fit_object
    
    smoothOrder = smooth_orders[sl_pick_idx]
    n_knots_default_pick = as.numeric(num_knots[sl_pick_idx] == "default")
    
  } else {
    
    hal_CV <- fit_hal(X = X, Y = Y, family = family,
                      return_x_basis = TRUE,
                      num_knots = hal9001:::num_knots_generator(
                        max_degree = ifelse(ncol(X) >= 20, 2, 3),
                        smoothness_orders = smoothOrder,
                        base_num_knots_0 = baseNumKnots,
                        base_num_knots_1 = baseNumKnots # max(100, ceiling(sqrt(n)))
                      )
    )
  }
  
  ## Unsersmoothing
  CV_nonzero_col <- which(hal_CV$coefs[-1] != 0)
  if (length(CV_nonzero_col) == 0){
    hal_fit = hal_CV
  }else{
    CV_basis_mat <- as.matrix(hal_CV$x_basis)
    CV_basis_mat <- as.matrix(CV_basis_mat[, CV_nonzero_col])
    
    hal_undersmooth <- undersmooth_hal(X, Y,
                                       fit_init = hal_CV,
                                       family = family)
    
    lambda_u_g = hal_undersmooth$lambda_under
    
    if(is.na(lambda_u_g)){
      hal_fit <- hal_CV
    } else {
      if(smoothOrderAdapt == TRUE & n_knots_default_pick == 1) {
        hal_fit <- fit_hal(X = X, Y = Y, family = family,
                           smoothness_orders = smoothOrder,
                           return_x_basis = TRUE,
                           fit_control = list(cv_select = FALSE),
                           lambda = lambda_u_g)
      } else {
        hal_fit <- fit_hal(X = X, Y = Y, family = family,
                           smoothness_orders = smoothOrder,
                           return_x_basis = TRUE,
                           fit_control = list(cv_select = FALSE),
                           lambda = lambda_u_g,
                           num_knots = hal9001:::num_knots_generator(
                             max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                             smoothness_orders = smoothOrder,
                             base_num_knots_0 = baseNumKnots,
                             base_num_knots_1 = baseNumKnots  
                           ))
      }
      
    }
  }
  
  # 3. prediction
  psi_hat <- sapply(curvePoints, function(a){ X_new <- X
  X_new$A = a
  mean(predict(hal_fit, new_data = X_new)) } )
  
  # 4. IC-based inference
  psi_hat_pnt_se <- IC_based_se(X, Y, hal_fit, curvePoints)
  
  # returns
  psi_hat_pnt <- cbind(curvePoints, matrix(psi_hat, ncol=1), psi_hat_pnt_se)
  
  colnames(psi_hat_pnt) <- c("a", "y_hat", "SE")
  
  psi_hat_pnt <- as.data.frame(psi_hat_pnt) %>% 
    mutate(ci_lwr = y_hat - 1.96 * SE,
           ci_upr = y_hat + 1.96 * SE)
  
  
  if(boundResults){
    # Setting bounds based on family
    if(family == "binomial") {
      bounds <- c(0, 1)
    } else {
      bounds <- c(min(Y), max(Y))
    }
    
    ## Applying bounds to each column
    psi_hat_pnt <- apply_bounds(psi_hat_pnt, "y_hat", bounds)
    psi_hat_pnt <- apply_bounds(psi_hat_pnt, "ci_lwr", bounds)
    psi_hat_pnt <- apply_bounds(psi_hat_pnt, "ci_upr", bounds)
  }


  return(psi_hat_pnt)
}



# Function to apply bounds
apply_bounds <- function(data, column, bounds) {
  data[, column] <- pmax(bounds[1], data[, column])
  data[, column] <- pmin(data[, column], bounds[2])
  return(data)
}
