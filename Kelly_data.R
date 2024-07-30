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
library(readr)
obs <- read_csv("~/Downloads/Kelly_sample_data.csv")
obs <- rename(obs, W = X)

head(obs)

### Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

### Plotting the correlation matrix
pairs(obs,
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth) # Smoothed regression lines


## 4. Estimate the Dose-Response Curve with HAL-based plugin estimator.
### 4.1 fit undersmoothed smoothness-order adaptive HAL
  # Note that fitting undersmoothed smoothness-order adaptive HAL takes the 
  # longest time to run.
DRC_fit_UAdaptive <- fit_UHAL_DRC(
  dat=obs, y_var_name="Y", trt_var_name="A", family="gaussian")

# Preview the fitted curve data
DRC_fit_UAdaptive$curve_est |> head()
## No CI outputted

# Try one without undersmoothing
# 4.2 Fit 1st smoothness order HAL-based DUC
DRC_fit <- fit_UHAL_DRC(
  dat=obs, y_var_name="Y", trt_var_name="A", family="gaussian", 
  Unsersmoothing = FALSE,
  smoothOrderAdapt = FALSE, smoothOrder = 1)

## 5. Visual the estimated curve
library(ggplot2)

curve_fit <- DRC_fit

df_ests <- curve_fit$curve_est
head(df_ests)

ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_point() + 
  geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr),  alpha=0.3) +
  theme_bw()

### Compare with GAM
library(gam)
model_fit <- gam(as.formula("Y ~ s(W) + s(A)"),
                 data = obs)

eval_points <- df_ests$a
psi_hat_pnt <- matrix(NA, nrow = length(eval_points), ncol = 5)
colnames(psi_hat_pnt) = c("a", "y_hat", "SE", "ci_lwr", "ci_upr")
for (i in 1:length(eval_points)) {
  psi_hat_pnt[i,1] = eval_points[i]
  
  X_new <- obs[x_names]
  X_new[, colnames(X_new)=='A'] = eval_points[i]
  
  predictions <- predict(model_fit, newdata = X_new, se.fit = TRUE)
  
  y_hat <-  mean(predictions)
  # SE <-  mean(predictions$se.fit)
  SE <- NA
  
  psi_hat_pnt[i,2] = y_hat
  psi_hat_pnt[i,3] = SE
  psi_hat_pnt[i,4] = y_hat - 1.96 * SE
  psi_hat_pnt[i,5] = y_hat + 1.96 * SE
} 

df_ests <- df_ests |>
  mutate(method = "HAL") |>
  rbind(psi_hat_pnt |>
          as.data.frame() |>
          mutate(method = "GAM"))

ggplot(df_ests, aes(x = a, y = y_hat, color = method)) +
  geom_line() +
  geom_point() + 
  geom_ribbon(aes(ymin=ci_lwr, ymax=ci_upr),  alpha=0.3) +
  theme_bw() 