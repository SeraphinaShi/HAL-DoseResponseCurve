
simu.num <- Sys.getenv("SIMU_NUM")
n <- Sys.getenv("SAMPLE_N")

half1 <- Sys.getenv("HALF1")
half2 <- Sys.getenv("HALF2")
grid_extra <- Sys.getenv("GRID_EXTRA")
  
print(simu.num)
print(n)



library(here)
library(dplyr)
library(glmnet)
library(hal9001)

R.files.sources <- c(list.files("R", pattern="*.R$", full.names=TRUE),
                     list.files("simu_R", pattern="*.R$", full.names=TRUE))
print(R.files.sources)

sapply(R.files.sources, source)


set.seed(123)


rep.num = 500

eval_points = seq(0,5,0.25)

set.seed(123)

##========True curve
psi0_all <- true_curve(simu.num = simu.num, N = 1e+07)
psi0_pnt <- psi0_all[psi0_all$a %in% eval_points,] 

save.image(file=here("data", "rdata", paste("simu", simu.num, "psi0.Rdata", sep="_")))

if(half1){
  ## ------1st order HAL
  set.seed(123)
  
  results_1 <- run_simu_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num)
  save.image(file=here("data", "rdata", paste("simu", simu.num, n, "first.Rdata", sep="_")))
  
  rm(results_1)
  
  ## ------0 order HAL
  set.seed(123)
  
  results_0 <- run_simu_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num, defualt_setting = T)
  save.image(file=here("data", "rdata", paste("simu", simu.num, n, "zero.Rdata", sep="_")))
  
  rm(results_0)
  
  ## ------gird
  set.seed(123)
  results_grid <- run_simu_scaled_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num)
  
  save.image(file=here("data", "rdata", paste("simu", simu.num, n, "grid.Rdata", sep="_")))
  
  rm(results_grid)
}

if(half2){
  ## -------apadt smoothness HAL
  set.seed(123)
  results_adapt <- run_simu_smoothness_adaptive_HAL_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num)
  save.image(file=here("data", "rdata", paste("simu", simu.num, n, "adapt.Rdata", sep="_")))
  
  rm(results_adapt)
  
  ## ------gam
  set.seed(123)
  results_gam <- run_simu_gam_poly_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num, method = "GAM")
  save.image(file=here("data", "rdata", paste("simu", simu.num, n, "GAM.Rdata", sep="_")))
  
  rm(results_gam)
  
  ## -------poly
  set.seed(123)
  results_poly <- run_simu_gam_poly_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num, method = "POLY")
  save.image(file=here("data", "rdata", paste("simu", simu.num, n, "poly.Rdata", sep="_")))
  
  rm(results_poly)
  
  ## -------poly
  set.seed(123)
  results_npcausal <- run_simu_npcausal_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num)
  save.image(file=here("data", "rdata", paste("simu", simu.num, n, "npcausal.Rdata", sep="_")))
  
  rm(results_npcausal)
}


if(grid_extra){
  set.seed(123)
  results_grid_extra <- run_simu_scaled_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num, grid_extra = T)
  
  save.image(file=here("data", "rdata", paste("simu", simu.num, n, "grid_extra.Rdata", sep="_")))
}
