## -----------------------------------------------------------------------------------------------------------------------
# define data generating systems
library(mvtnorm)

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


generate_data_2 <- function(n, a=NA){
  # exogenous variables
  U_W <- rnorm(n, 0, 1)
  U_A <- rnorm(n, 0, 1.3)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <- 2.5 - 0.5*W + U_A
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }
  
  
  Y <- as.numeric(U_Y < plogis(-7 + 3*W + 5*sin(1.25*A^1.5) + 5*A + 3 * W * A ))
  
  # data frame
  O <- data.frame(W, A, Y)
  return(O)
}

generate_data_3 <- function(n, a=NA){
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
  
  Y <- as.numeric(U_Y < plogis(-10 - 3*W + 4*A + 5*sin((0.8*A)^2-2.6)*as.numeric(A > 2)) )
  
  # data frame
  O <- data.frame(W, A, Y)
  return(O)
}


generate_data_4 <- function(n, a=NA){
  # exogenous variables
  U_W <- rnorm(n, 0, 1)
  U_A <- rnorm(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- U_W
  
  if(is.na(a)){
    A <-  2.5 - 0.5*W + U_A
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }
  
  Y <- as.numeric(U_Y < plogis(-6 + W + 3.5*A*as.numeric(A >= 2) - 4*A*as.numeric(A >= 4) - 0.5 * W * A ))
  
  # data frame
  O <- data.frame(W, A, Y)
  return(O)
}


generate_data_5 <- function(n, a=NA){
  # exogenous variables
  U_W_sigma = c(1, 1, 1.5, 0.5, 0.3, 1, 1, 5, 0.5, 2)
  
  U_W_r = matrix(c(1,    0.8,  0.9, 0.2, 0.1, 0,   0,   0.3, 0.1, 0.96,
                   0.8,  1,    0.8, 0.3, 0.1, 0,   0,   0.8, 0.1, 0.8,   
                   0.9,  0.8,  1,   0.1, 0.5, 0.7, 0.1, 0.7, 0.8, 0.9, 
                   0.2,  0.3,  0.1, 1,   0.5, 0.3, 0,   0,   1,   1,   
                   0.1,  0.1,  0.5, 0.5, 1,   0,   0,   0,   0.1, 0.1,
                   0,    0,    0.7, 0.3, 0,   1,   0.8, 0,   0.4, 0.1,
                   0,    0,    0.1, 0,   0,   0.8, 1,   0.9, 0.8, 0.5,
                   0.3,  0.8,  0.7, 0,   0,   0,   0.9, 1,   0.7, 0.9,
                   0.1,  0.1,  0.8, 1,   0.1, 0.4, 0.8, 0.7, 1,   0.8,
                   0.96, 0.8,  0.9, 1,   0.1, 0.1, 0.5, 0.9, 0.8, 1), nrow = 10)
  
  U_W_var = matrix(nrow = 10, ncol = 10)
  
  for (i in 1:10) {
    for (j in 1:10) {
      U_W_var[i,j] = U_W_r[i,j] * U_W_sigma[i] * U_W_sigma[j]
    }
  }
  
  U_W <- rmvnorm(n, mean = c(0, 0, 5, 1, 1, 0.3, 1, 0, 3, 2), sigma = U_W_var)
  U_A <- rnorm(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- data.frame(U_W)
  colnames(W) <- paste0('W', 1:ncol(W))
  
  if(is.na(a)){
    A <-  (0.1*W$W1 + 0.2*W$W2  + 0.5*W$W3 + 0.15*W$W4 - 0.05*W$W5 - 0.01*W$W3*W$W5 + U_A)
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }
  
  Y <- as.numeric(U_Y < plogis(-7 - W$W1 + 2*W$W2 - 1*W$W1 * W$W3 - 0.5*W$W4 - 0.2 * W$W6^2 + 0.01 * W$W8 *W$W10 + 4*A + 0.5*A*W$W2))
  mean(Y)
  
  # data frame
  O <- W
  O$A = A
  O$Y = Y
  return(O)
}


generate_data_6 <- function(n, a=NA){
  # exogenous variables
  U_W_sigma = c(1, 1, 1.5, 0.5, 0.3, 1, 1, 5, 0.5, 2)
  
  U_W_r = matrix(c(1,    0.8,  0.9, 0.2, 0.1, 0,   0,   0.3, 0.1, 0.96,
                   0.8,  1,    0.8, 0.3, 0.1, 0,   0,   0.8, 0.1, 0.8,   
                   0.9,  0.8,  1,   0.1, 0.5, 0.7, 0.1, 0.7, 0.8, 0.9, 
                   0.2,  0.3,  0.1, 1,   0.5, 0.3, 0,   0,   1,   1,   
                   0.1,  0.1,  0.5, 0.5, 1,   0,   0,   0,   0.1, 0.1,
                   0,    0,    0.7, 0.3, 0,   1,   0.8, 0,   0.4, 0.1,
                   0,    0,    0.1, 0,   0,   0.8, 1,   0.9, 0.8, 0.5,
                   0.3,  0.8,  0.7, 0,   0,   0,   0.9, 1,   0.7, 0.9,
                   0.1,  0.1,  0.8, 1,   0.1, 0.4, 0.8, 0.7, 1,   0.8,
                   0.96, 0.8,  0.9, 1,   0.1, 0.1, 0.5, 0.9, 0.8, 1), nrow = 10)
  
  U_W_var = matrix(nrow = 10, ncol = 10)
  
  for (i in 1:10) {
    for (j in 1:10) {
      U_W_var[i,j] = U_W_r[i,j] * U_W_sigma[i] * U_W_sigma[j]
    }
  }
  
  U_W <- rmvnorm(n, mean = c(0, 0, 5, 1, 1, 0.3, 1, 0, 3, 2), sigma = U_W_var)
  U_A <- rnorm(n, 0, 1)
  U_Y <- runif(n, 0, 1)
  
  # endogenous variables
  W <- data.frame(U_W)
  colnames(W) <- paste0('W', 1:ncol(W))
  
  if(is.na(a)){
    A <-  (0.1*W$W1 + 0.2*W$W2  + 0.5*W$W3 + 0.15*W$W4 - 0.05*W$W5 - 0.01*W$W3*W$W5 + U_A)
    A[A<=0] = 0
    A[A>=5] = 5
  } else {
    A <- rep(a, n)
  }
  
  Y <- as.numeric(U_Y < plogis(-9 - W$W1 + 0.7*W$W2 - 0.3*W$W1 * W$W3 - 0.4*W$W4 - 0.2 * W$W6^2 + 0.01 * W$W8 *W$W10 + 4*A + 3*sin((0.8*A)^2)) )
  mean(Y)
  
  # data frame
  O <- W
  O$A = A
  O$Y = Y
  return(O)
}

expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

generate_data_7 <- function(n, a=NA){
  W1 <- runif(n)
  W2 <- rbinom(n, size=1, prob =0.7)
  W3 <- rnorm(W1, 0.25*exp(2*W1))
  
  v_W <- exp(1 + 2*W1*expit(W3))
  mu_W <- expit(0.03 - 0.8*log(1+W2) + 0.9*exp(W1)*W2 - 0.4*atan(W3+2)*W2*W1)
  if(is.na(a)){
    A <- rbeta(n, v_W*mu_W, v_W*(1-mu_W))
  } else {
    A <- a
  }
  
  Q <- expit(-2 + 1.5*A + 5*A^3 - 2.5*W1 + 0.5*A*W2 - log(A)*W1*W2 + 0.5*(A^(3/4))*W1*W3)
  Y <- rbinom(n, size = 1, prob = Q)
  
  O <- data.frame(W1, W2,W3, A, Y)
  return(O)
}


DGS <- list("simu1" = generate_data_1,
            "simu2" = generate_data_2,
            "simu3" = generate_data_3,
            "simu4" = generate_data_4,
            "simu5" = generate_data_5,
            "simu6" = generate_data_6,
            "simu7" = generate_data_7)


gene_data <- function(simu.num =1, n=100, a=NA){
  data = DGS[[simu.num]](n=n, a=a)
  return(data)
}

## ----true_psi_sys1-------------------------------------------------------------------------------------------------------------------
# Getting trul value of psi
# 

true_curve <- function(simu.num = 1, a.pts.dist = 0.05, N = 1e+07){
  
  psi0_a <- c()
  a.min = 0
  
  if(simu.num %in% 1:6){
    a.max = 5
  } else {
    a.max = 1
  }
  
  a.vals <- seq(a.min, a.max, a.pts.dist)
  
  for (i in 1:length(a.vals)) {
    a <- a.vals[i]

    data_a <- DGS[[simu.num]](n=N, a=a) # gene_data(simu.num=simu.num, n=N, a=a)
    psi0_a[i] <- mean(data_a$Y)
  }
  
  psi0 <- data.frame(a=a.vals, psi0 = psi0_a)
  return(psi0)
}

# psi0 <- true_curve_0(simu.num = 5)
# plot(psi0$a, psi0$psi0)
# 
# 
# psi1 <- true_curve(simu.num = 5)
# plot(psi1$a, psi1$psi0)

# 
# true_curve <- function(simu.num = 1, a.pts.dist = 0.05){
#   
#   psi0_a <- c()
#   a.min = 0
#   
#   if(simu.num %in% 1:6){
#     a.max = 5
#   } else {
#     a.max = 1
#   }
# 
#   a.vals <- seq(a.min, a.max, a.pts.dist)
#     
#   for (i in 1:length(a.vals)) {
#     a <- a.vals[i]
#     EW = 0
#     
#     if(simu.num == 1){
#       psi0_a[i] = plogis(-5 + EW + 2.25*a - 0.5 * EW * a)
#     }
# 
#     if(simu.num %in% c(2, 3)){
#       data_a <- gene_data(simu.num=simu.num, n=1e+07, a=a)
#       psi0_a[i] <- mean(data_a$Y)
#     }
#     
#     if(simu.num == 4){
#       if(a < 2 | a >= 4){
#         psi0_a[i] = 0
#       } else {
#         psi0_a[i] = plogis(-6 + EW + 3.5*a*as.numeric(a >= 2) - 4*a*as.numeric(a >= 4) - 0.5 * EW * a )
#       }
#     }
#     
#     if(simu.num == 5){
#       EWs <- c(0, 0, 5, 1, 1, 0.3, 1, 0, 3, 2)
#       psi0_a[i] = plogis(-7 - EWs[1] + 2*EWs[2] - 1*EWs[1] * EWs[3] - 0.5*EWs[4] - 0.2 * EWs[6]^2 + 0.01 * EWs[8] *EWs[10] + 4*a + 0.5*a*EWs[2])
#     }
#     
#     if(simu.num == 6){
#       EWs <- c(0, 0, 5, 1, 1, 0.3, 1, 0, 3, 2)
#       psi0_a[i] =  plogis(-9 - EWs[1] + 0.7*EWs[2] - 0.3*EWs[1] * EWs[3] - 0.4*EWs[4] - 0.2 * EWs[6]^2 + 0.01 * EWs[8] *EWs[10] + 4*a + 3*sin((0.8*a)^2))
#     }
#     
#     if(simu.num == 7){
#       data_a <- gene_data(simu.num=7, n=1e+07, a=a)
#       psi0_a[i] <- mean(data_a$Y)
#     }
#   }
#   
#   psi0 <- data.frame(a=a.vals, psi0 = psi0_a)
#   return(psi0)
# }
# 
