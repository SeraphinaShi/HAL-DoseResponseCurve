

library(here)
library(dplyr)
library(glmnet)
library(hal9001)

R.files.sources <- c(list.files("R", pattern="*.R$", full.names=TRUE),
                     list.files("simu_R", pattern="*.R$", full.names=TRUE))
print(R.files.sources)

sapply(R.files.sources, source)

####################################################################
simu.num = 4
nn = 200
load(here("data", "rdata", paste0("simu_", simu.num, "_", nn, "_npcausal.RData")))
sum(is.na(results_npcausal$result_summary$y_hat))

na_rounds <- which(sapply(1:500, function(x) any(is.na(results_npcausal$all_results[[x]]$y_hat))))
na_rounds

for(i in 1:na_rounds) {
  print(results_npcausal$all_results[[na_rounds[i]]])
}

# ====================
nn = 500
load(here("data", "rdata", paste0("simu_", simu.num, "_", nn, "_npcausal.RData")))
sum(is.na(results_npcausal$result_summary$y_hat))

na_rounds <- which(sapply(1:500, function(x) any(is.na(results_npcausal$all_results[[x]]$y_hat))))
na_rounds

for(i in 1:na_rounds) {
  print(results_npcausal$all_results[[na_rounds[i]]])
}

# ====================
nn = 1000
load(here("data", "rdata", paste0("simu_", simu.num, "_", nn, "_npcausal.RData")))
sum(is.na(results_npcausal$result_summary$y_hat))

na_rounds <- which(sapply(1:500, function(x) any(is.na(results_npcausal$all_results[[x]]$y_hat))))
na_rounds
# ====================
nn = 500
load(here("data", "rdata", paste0("simu_", simu.num, "_", nn, "_npcausal.RData")))
sum(is.na(results_npcausal$result_summary$y_hat))

na_rounds <- which(sapply(1:500, function(x) any(is.na(results_npcausal$all_results[[x]]$y_hat))))
na_rounds



####################################################################
n = 200

set.seed(123)
results_npcausal <- run_simu_npcausal_rep(simu.num, eval_points, y_type = "binomial", n=n, rounds=5)

# round 6

obs <-  DGS(simu.num, n)

l_names <-  names(obs)[! names(obs) %in% c('Y', 'A')]
x = obs %>% select(all_of(l_names)) %>% mutate_if(sapply(., is.factor), as.numeric) 
a = obs$A
y = obs$Y

if (dim(x)[2]) {
  x$ones <- rep(1, nrow(x))
}

ce.res <- ctseff(y, a, x, bw.seq = seq(.2, 2, length.out = 20), eval_points)
warnings()

# line 114 in the ctseff_from_npcausal gives the NA prediction
# no warnings
# 113 # estimate effect curve with optimal bandwidth
# 114 est <- approx(locpoly(a, pseudo.out, bandwidth = h.opt), xout = a.vals)$y

psi_hat_pnt <- ce.res$res
names(psi_hat_pnt) <- c("a", "y_hat", "SE", "ci_lwr", "ci_upr")
print(psi_hat_pnt)


hist(a)
summary(a)
sort(a)

plot(a, y)

