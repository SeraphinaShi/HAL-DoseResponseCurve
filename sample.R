# This file give an example of how to use the HAL-based plugin estimator 
# for estimating marginal (causal) dose response curve.

# 1. Load the libraries 

## Install needed packages
# install.packages("here")
# install.packages("dplyr")
# install.packages("glmnet")
# install.packages("hal9001")

library(here)
library(dplyr)
library(glmnet)
library(hal9001)


# 2. Load the functions downloaded from 
# https://github.com/SeraphinaShi/HAL-DoseResponseCurve/tree/main/R
R.files.sources <- c(list.files("R", pattern="*.R$", full.names=TRUE))
print(R.files.sources)

sapply(R.files.sources, source)


# 3.  load data
# The input data should be numeric, containing a treatment column, a outcome
# column, and all other columns that will be used in fitting the data generation
# distribution.
# Here we simulate a data. W is a baseline (confounder) variable; A is a
# treatment variable, depending on W; and Y is a binary outcome, depending on
# both A and W.
set.seed(145)

generate_data <- function(n, a = NA){
  
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
  
  Y <- as.numeric(U_Y < plogis(-3 + 0.5*W + 1.25*A - 0.5 * W * A ))
  
  # data frame
  O <- data.frame(W, A, Y)
  return(O)
}

obs <- generate_data(n = 500)

head(obs)

# 3. Estimate the Dose-Response Curve with HAL-based plugin estimator.
## 3.1. fit undersmoothed smoothness-order adaptive HAL
  # Note that fitting undersmoothed smoothness-order adaptive HAL takes the 
  # longest time to run.
DRC_fit_UAdaptive <- fit_UHAL_DRC(
  dat=obs, y_var_name="Y", trt_var_name="A", family="binomial")

## 3.2. fit zero-smoothness order HAL
DRC_fit_zero <- fit_UHAL_DRC(
  dat=obs, y_var_name="Y", trt_var_name="A", family="binomial", 
  Unsersmoothing = FALSE,
  smoothOrderAdapt = FALSE, smoothOrder = 0)

## 3.3. fit undersmoothed 1st order smoothness HAL with baseNumKnots = 20
##  and provide the points on the curve that we want to estimate at
DRC_fit_custom <- fit_UHAL_DRC(
  dat=obs, y_var_name="Y", trt_var_name="A", family="binomial", 
  curvePoints = c(0.05, 0.1, 0.3, 0.5, 0.6, 0.75, 0.8, 0.9,0.99,1))
  

# 4. Visual the estimated curve
library(ggplot2)

df_ests <- DRC_fit_UAdaptive$curve_est
head(df_ests)

ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_point() + 
  geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr),  alpha=0.3) +
  theme_bw()


# 5. compare with the true curve
# Note: this is just for the simulation

## get true curve
a.vals <- seq(0, 5, length.out = 20)
psi0_a <- c()
for (i in 1:length(a.vals)) {
  a_i <- a.vals[i]
  data_a <- generate_data(n=1e+7, a=a_i) 
  psi0_a[i] <- mean(data_a$Y)
}
psi0 <- data.frame(a=a.vals, psi0 = psi0_a)

## Compare
ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_point() + 
  geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr),  alpha=0.3) + 
  geom_point(data = psi0, aes(x=a, y=psi0), color = "red") + 
  geom_line(data = psi0, aes(x=a, y=psi0), color = "red") + 
  theme_bw()


#=============================================================
# 6. Another sample with continuous outcome
generate_data_2 <- function(n, a = NA) {
  
  # Exogenous variables
  U_W <- rnorm(n, mean = 0, sd = 1)
  U_A <- rnorm(n, mean = 0, sd = 0.8)
  U_Y <- rnorm(n, mean = 0, sd = 0.5)
  
  # Endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <- 0.5 * W + U_A
  } else {
    A <- rep(a, n)
  }
  
  # Modify Y to be continuous between 0 and 1 using a logistic function
  logit_Y <- 0.5 * A^2 + A + 0.3 * W + U_Y
  Y <- 1 / (1 + exp(-logit_Y))
  
  # Data frame
  O <- data.frame(W, A, Y)
  return(O)
}

## 6.1. Generate data with 500 samples
set.seed(123)
obs <- generate_data_2(n = 500)

head(obs)
### Plot W vs A
plot(obs$W, obs$A, 
     main = "W vs A", xlab = "W (confounder)", ylab = "A (treatment)", 
     pch = 19, col = rgb(0, 0, 1, 0.5))

### Plot A vs Y
plot(obs$A, obs$Y, 
     main = "A vs Y", 
     xlab = "A (treatment)", ylab = "Y (outcome)", 
     pch = 19, col = rgb(1, 0, 0, 0.5))

## 6.2. Estimate the Dose-Response Curve with HAL-based plugin estimator.
### 6.2.1 fit undersmoothed smoothness-order adaptive HAL
  # Note that fitting undersmoothed smoothness-order adaptive HAL takes the 
  # longest time to run.
DRC_fit_UAdaptive <- fit_UHAL_DRC(
  dat=obs, y_var_name="Y", trt_var_name="A", family="gaussian")

## 6.2.2. fit zero-smoothness order HAL
DRC_fit_zero <- fit_UHAL_DRC(
  dat=obs, y_var_name="Y", trt_var_name="A", family="gaussian", 
  Unsersmoothing = FALSE,
  smoothOrderAdapt = FALSE, smoothOrder = 0)

## 6.2.3. fit undersmoothed 1st order smoothness HAL with baseNumKnots = 20
##  and provide the points on the curve that we want to estimate at
DRC_fit_custom <- fit_UHAL_DRC(
  dat=obs, y_var_name="Y", trt_var_name="A", family="gaussian", 
  curvePoints = c(-3,-2.5,-2,-1.5,-1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 1.5, 2, 3))


# 6.3. Visual the estimated curve
library(ggplot2)

df_ests <- DRC_fit_UAdaptive$curve_est
head(df_ests)

ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_point() + 
  geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr),  alpha=0.3) +
  theme_bw()


## 6.4. compare with the true curve
  # Note: this is just for the simulation

### 6.4.1. Get true dose-response curve
a.vals <- seq(min(obs$A), max(obs$A), length.out = 100)
Y.values <- sapply(a.vals, function(a_i) {
  data_a <- generate_data_2(n = 1e+7, a = a_i)
  mean(data_a$Y)
})
psi0 <- data.frame(a = a.vals, psi0 = Y.values)

## Compare
ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_point() + 
  geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr),  alpha=0.3) + 
  geom_point(data = psi0, aes(x=a, y=psi0), color = "red") + 
  geom_line(data = psi0, aes(x=a, y=psi0), color = "red") + 
  theme_bw()

