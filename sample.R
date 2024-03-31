library(here)
library(dplyr)
library(glmnet)
library(hal9001)

set.seed(123)

# data

n = 200

W <- rnorm(n, 0, 1)

A <-  0.5 + 0.3*W + rnorm(n, 0, 0.2)
A[A<=0] = 0
A[A>=1] = 1


U_Y <- runif(n, 0, 1)
Y <- as.numeric(U_Y < plogis(-0.5 + 0.5*W + 1.25*A - 0.25 * W * A ))

obs <- data.frame(W,A,Y) 
# NOTE: the treatment variable should be "A" and as a column in the X matrix for codes to work



# fit HAL

R.files.sources <- c(list.files("R", pattern="*.R$", full.names=TRUE))
print(R.files.sources)

sapply(R.files.sources, source)

DRC_fit <- fit_UHAL_DRC(dat=obs, y_var_name="Y", trt_var_name="A", family="binomial", 
                        curvePoints = NA, 
                        smoothOrderAdapt = FALSE, 
                        smoothOrder = 1,  baseNumKnots = 20, boundResults = TRUE)

family="binomial"
curvePoints = NA
smoothOrderAdapt = FALSE
smoothOrder = 1
baseNumKnots = 20
boundResults = TRUE

head(DRC_fit)
  


