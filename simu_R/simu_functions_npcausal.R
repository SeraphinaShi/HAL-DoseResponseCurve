

##############################################################
# This function runs the simulation with given:
#       - data generating function
#       - sample size n
# returns the estimated ATE and empirical 95% CI
##############################################################
run_simu_1round_npcausal <- function(simu.num, eval_points, y_type, n){
  
  obs <-  DGS(simu.num, n)
  
  l_names <-  names(obs)[! names(obs) %in% c('Y', 'A')]
  x = obs %>% select(all_of(l_names)) %>% mutate_if(sapply(., is.factor), as.numeric) 
  a = obs$A
  y = obs$Y
  
  if (dim(x)[2]) {
    x$ones <- rep(1, nrow(x))
  }

  ce.res <- ctseff(y, a, x, bw.seq = seq(.2, 2, length.out = 20), eval_points)
  
  
  psi_hat_pnt <- ce.res$res
  names(psi_hat_pnt) <- c("a", "y_hat", "SE", "ci_lwr", "ci_upr")
  
  if(y_type == "binomial"){
    bounds <- c(0, 1)
  } else {
    bounds <- c(min(Y), max(Y))
  }
  psi_hat_pnt[,"y_hat"] <- pmax(bounds[1], psi_hat_pnt[,"y_hat"])
  psi_hat_pnt[,"y_hat"] <- pmin(psi_hat_pnt[,"y_hat"], bounds[2])
  
  psi_hat_pnt[,"ci_lwr"] <- pmax(bounds[1], psi_hat_pnt[,"ci_lwr"])
  psi_hat_pnt[,"ci_lwr"] <- pmin(psi_hat_pnt[,"ci_lwr"], bounds[2])
  
  psi_hat_pnt[,"ci_upr"] = pmax(bounds[1], psi_hat_pnt[,"ci_upr"])
  psi_hat_pnt[,"ci_upr"] <- pmin(psi_hat_pnt[,"ci_upr"], bounds[2])

  return(psi_hat_pnt)
}



##############################################################
# This function runs the simulation for B rounds with given:
#       - data generating function
#       - sample size n
#       - number of simulations: B

# returns the estimated ATE, empirical 95% CI, and coverage rates
##############################################################
run_simu_npcausal_rep <- function(simu.num, eval_points, y_type, n, rounds){
  
  result_list <- list()
  
  for(r in 1:rounds){
    print(paste0("round ", r))
    result <- tryCatch({
      run_simu_1round_npcausal(simu.num, eval_points, y_type, n=n)
    }, error = function(e) {
      print(paste0("Error: ", e$message))
      NULL
    })
    
    while(is.null(result)) {
      print('retry with a new generated data')
      result <- tryCatch({
        run_simu_1round_npcausal(simu.num, eval_points, y_type, n=n)
      }, error = function(e) {
        print(paste0("Error: ", e$message))
        NULL
      })
    }
    
    result_list[[r]] <- result
  }
  
  result_all <-  do.call("rbind", result_list) %>% as.data.frame()
  result_all <- merge(as.data.frame(psi0_pnt), result_all, by=c("a"))
  
  result_summary <- result_all %>% 
    filter(SE != 0, ! is.na(y_hat)) %>% 
    mutate(SE = SE/sqrt(n)) %>%
    mutate(bias = abs(y_hat - psi0),
           bias_se_ratio = bias / SE,
           cover_rate = as.numeric(ci_lwr <= psi0 & psi0 <= ci_upr),
           MSE =(y_hat - psi0)^2 ) %>% 
    group_by(a) %>% 
    mutate(oracal_SE = sqrt(var(y_hat)),
           oracal_bias_se_ratio = bias / oracal_SE,
           oracal_ci_lwr = y_hat - 1.96 * oracal_SE,
           oracal_ci_upr = y_hat + 1.96 * oracal_SE,
           oracal_cover_rate = as.numeric(oracal_ci_lwr <= psi0 & psi0 <= oracal_ci_upr)) %>%
    summarise(across(where(is.numeric), mean)) %>% 
    ungroup() %>%
    mutate(method = "npcausal")
  
  results <- list(result_summary = result_summary,
                  all_results = result_list)
  
  return(results)
}


