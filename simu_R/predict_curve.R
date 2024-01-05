###############################################################################
#'  Estimate the marginal dose response curve for continuous treatment
#'
#' @details \code{HAL_DRC} is used to estimate the mean outcomes in a population had all subjects received given levels of a continuous (unconfounded) treatment.
#' It performs the following steps:
#'     1). fits HAL with selected tunning parameters
#'     2). obtain curve estimations 
#'     3). calculate the influence curve based confidence interval
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param Y A \code{numeric} vector of observations of the outcome variable.
#' @param pts A \code{numeric} vector of treatment values 
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.
#' @param HAL_smooth A \code{character} $\in$ {"adaptive", "zero", "first"}, defining the smoothness order used in HAL
#' @param HAL_undersmooth A \code{boolean}

HAL_DRC <- function(X, Y, points_curve, family = "binomial", HAL_smooth = "adaptive", HAL_undersmooth = T){
  y_name <- "Y"
  x_names <- names(X)
  
  Y <- as.numeric(as.matrix(obs %>% select(all_of(y_name))))
  X <- obs %>% 
    select(all_of(x_names)) %>% 
    mutate_if(sapply(., is.factor), as.numeric)
  
  if(HAL_smooth == "adaptive"){
    dSL_fit <- fit_SL_smoothness_adaptive_HAL(X, Y, y_type)
  } else if (HAL_smooth == "zero"){
    
  } else if (HAL_smooth == "first"){
    
  }
  
  
  
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


