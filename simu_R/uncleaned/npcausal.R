library(here)
library(data.table)
library(dplyr)

library(SuperLearner)

generate_data_1 <- function(n, a=NA){
  # exogenous variables
  U_W <- rnorm(n, 0, 1)
  U_A <- rnorm(n, 0, 2)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <-  2 - 0.5*W + U_A
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }
  
  
  Y <- as.numeric(U_Y < plogis(-5 + W + 2.25*A - 0.5 * W * A ))
  
  # data frame
  O <- data.frame(W, A, Y)
  return(O)
}

load(here("data", "rdata", "02_simu_V5_sys1_psi0.RData"))

n = 200

set.seed(123)

obs <- generate_data_1(n)

## -----------------------------------------------------------------------------------------------------------------------


###INPUT:l is an n*p matrix,a and y are vectors of length n 
### l=matrixofcovariates 
### a=vectoroftreatmentvalues 
### y=vectorofobservedoutcomes
l_names <-  names(obs)[! names(obs) %in% c('Y', 'A')]
l = obs %>% select(all_of(l_names)) %>% mutate_if(sapply(., is.factor), as.numeric) 
a = obs$A
y = obs$Y

if (dim(l)[2]) {
  l$ones <- rep(1, nrow(l))
}

source(here("npcausal_R", "ctseff.R"))
source(here("npcausal_R", "plot.ctseff.R"))

x = l
ce.res <- ctseff(y, a, x, bw.seq = seq(.2, 2, length.out = 20))

plot.ctseff(ce.res)
plot(ce.res$bw.risk$bw, ce.res$bw.risk$risk)



library(ggplot2)
results_df <- data.frame("a" = ce.res$res$a.vals, "est" = ce.res$res$est, "ci_lwr" = ce.res$res$ci.ll, "ci_upr" = ce.res$res$ci.ul)
results_df_1 <- merge(results_df, psi0_line, by = "a", all.x = T)

p_est_avg <- ggplot(data=results_df_1, aes(x=a)) +
    geom_line(data = psi0_line ,aes(x = a, y=psi0), alpha = 0.5, color="darkgrey") +
    geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr), width=0.7, alpha=0.1) +
    geom_point(aes(y=psi0), color = "black") +
    geom_point(aes(y=est), shape=17, size=2, alpha= 0.7) +
    ylim(c(0,1)) +
    labs(x="a", y="E[Y(a)]", title = "Estimation") +
    theme_bw() +
    theme(legend.box = "horizontal")
p_est_avg
