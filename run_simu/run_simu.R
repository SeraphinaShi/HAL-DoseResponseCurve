
simu.num <- Sys.getenv("SIMU_NUM")
n <- Sys.getenv("SAMPLE_N")

print(simu.num)
print(n)

library(here)
library(dplyr)
library(glmnet)
library(hal9001)

R.files.sources <- c(list.files("R", pattern="*.R$", full.names=TRUE),
                     list.files("simu_R", pattern="*.R$", full.names=TRUE))

sapply(R.files.sources, source)


set.seed(123)


rep.num = 2

eval_points = seq(0,5,0.25)

set.seed(123)

generate_data_func <- DGS[[simu.num]]


##========True curve
psi0_all <- true_curve(simu.num = simu.num, N = 1000)
psi0_pnt <- psi0_all[psi0_all$a %in% eval_points,] 

save.image(file=here("data", "rdata", paste("simu", simu.num, "psi0.Rdata", sep="_")))

## ------1st order HAL
set.seed(123)

results_1 <- run_simu_rep(generate_data_func, eval_points, y_type = "binomial", n=n, rounds=rep.num)
save.image(file=here("data", "rdata", paste("simu", simu.num, n, "first.Rdata", sep="_")))

## ------0 order HAL
rm(results_1)

set.seed(123)

results_0 <- run_simu_rep(generate_data_func, eval_points, y_type = "binomial", n=n, rounds=rep.num, defualt_setting = T)
save.image(file=here("data", "rdata", paste("simu", simu.num, n, "zero.Rdata", sep="_")))

## ------gird
rm(results_0)

set.seed(123)
results_grid <- run_simu_scaled_rep(generate_data_func, eval_points, y_type = "binomial", n=n, rounds=rep.num)

save.image(file=here("data", "rdata", paste("simu", simu.num, n, "grid.Rdata", sep="_")))

## -------apadt smoothness HAL
rm(results_grid)

set.seed(123)
results_adapt <- run_simu_smoothness_adaptive_HAL_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=rep.num)
save.image(file=here("data", "rdata", paste("simu", simu.num, n, "adapt.Rdata", sep="_")))


## ------gam
set.seed(123)
results_gam <- run_simu_gam_poly_rep(generate_data_func, eval_points, y_type = "binomial", n=n, rounds=rep.num, method = "GAM")
save.image(file=here("data", "rdata", paste("simu", simu.num, n, "GAM.Rdata", sep="_")))

## -------poly
set.seed(123)
results_poly <- run_simu_gam_poly_rep(generate_data_func, eval_points, y_type = "binomial", n=n, rounds=rep.num, method = "POLY")
save.image(file=here("data", "rdata", paste("simu", simu.num, n, "poly.Rdata", sep="_")))

## -------poly
set.seed(123)
results_npcausal <- run_simu_npcausal_rep(generate_data_func, eval_points, y_type = "binomial", n=n, rounds=rep.num)
save.image(file=here("data", "rdata", paste("simu", simu.num, n, "npcausal.Rdata", sep="_")))
