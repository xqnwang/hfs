#----------------------------------------------------------------------
# Toy example: Subset selection with shrinkage under unbiasedness assumptions
#----------------------------------------------------------------------
library(fable)
library(dplyr)
library(tidyr)
library(tsibble)
library(tsibbledata)
library(lubridate)
library(purrr)
library(ggplot2)

source('R/sourceDir.R')
sourceDir("R", recursive = TRUE)
MinT_G <- function(S, W){
  solve(t(S) %*% solve(W) %*% S) %*% t(S) %*% solve(W)
}
obj_MP <- function(fc, S, W, G){
  t(fc - S %*% G %*% fc) %*% solve(W) %*% (fc - S %*% G %*% fc) |> as.numeric()
}
obj_MinT <- function(S, W, G){
  S %*% G %*% W %*% t(G) %*% t(S) |> diag() |> sum()
}

S <- rbind(rep(1, 4),
           c(1, 1, 0, 0), c(0, 0, 1, 1),
           diag(1, 4))

W_I <- diag(7)
G_ols <- MinT_G(S, W_I)

W_s <- S %*% matrix(1, nrow = NCOL(S), ncol = 1) |> as.vector() |> diag()
G_s <- MinT_G(S, W_s)

G_bu <- cbind(rep(0, 4), rep(0, 4), rep(0, 4), diag(1, 4))
G_td <- cbind(c(0.2, 0.4, 0.2, 0.2), matrix(0, nrow = 4, ncol = 6))

#----------------
# example 1: the first base forecast in the middle level performs poorly
fc1 <- c(10, 4, 4, 2, 4, 2, 2)
# OLS?
rec1 <- mip_l0(fc1, S, W_I, G_bench = NULL, 
              lambda_0 = 0, lambda_1 = 0, lambda_2 = 0, M = 10, 
              solver = "gurobi")
rec1$G
## reconciled forecast
S %*% rec1$G %*% fc1
S %*% G_ols %*% fc1
## objective value
### (\hat{y} - SG\hat{y})' W^{-1} (\hat{y} - SG\hat{y})
t(fc1 - S %*% rec1$G %*% fc1) %*% solve(W_I) %*% (fc1 - S %*% rec1$G %*% fc1)
t(fc1 - S %*% G_ols %*% fc1) %*% solve(W_I) %*% (fc1 - S %*% G_ols %*% fc1)
### tr(GWG')
rec1$G %*% W_I %*% t(rec1$G) |> diag() |> sum()
G_ols %*% W_I %*% t(G_ols) |> diag() |> sum()
### SGWG'S'
S %*% rec1$G %*% W_I %*% t(rec1$G) %*% t(S) |> diag() |> sum()
S %*% G_ols %*% W_I %*% t(G_ols) %*% t(S) |> diag() |> sum()

# No shrinkage
rec11 <- mip_l0(fc1, S, W_I, G_bench = NULL, 
               lambda_0 = 0.1, lambda_1 = 0, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec11$G

# Shrink toward G_bu
rec12 <- mip_l0(fc1, S, W_I, G_bench = G_bu, 
               lambda_0 = 0.1, lambda_1 = 0.2, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec12$G

rec13 <- mip_l0(fc1, S, W_I, G_bench = G_bu, 
               lambda_0 = 0.1, lambda_1 = 0, lambda_2 = 0.2, M = NULL, 
               solver = "gurobi")
rec13$G

# Shrink toward G_td
rec14 <- mip_l0(fc1, S, W_I, G_bench = G_td, 
               lambda_0 = 0.1, lambda_1 = 0.1, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec14$G

rec15 <- mip_l0(fc1, S, W_I, G_bench = G_td, 
               lambda_0 = 0.1, lambda_1 = 1, lambda_2 = 1, M = NULL, 
               solver = "gurobi")
rec15$G

rec16 <- mip_l0(fc1, S, W_I, G_bench = G_td, 
               lambda_0 = 10, lambda_1 = 20, lambda_2 = 20, M = NULL, 
               solver = "gurobi")
rec16$G

rec16 <- mip_l0(fc1, S, W_I, G_bench = G_td, 
                lambda_0 = 0, lambda_1 = 20, lambda_2 = 20, M = NULL, 
                solver = "gurobi")
rec16$G

#----------------
# example 2: the first two base forecasts in the bottom level perform poorly
fc2 <- c(10, 6, 4, 1, 2, 2, 2)
# OLS?
rec2 <- mip_l0(fc2, S, W_I, G_bench = NULL, 
              lambda_0 = 0, lambda_1 = 0, lambda_2 = 0, M = 1, 
              solver = "gurobi")
rec2$G

# No shrinkage
rec21 <- mip_l0(fc2, S, W_I, G_bench = NULL, 
               lambda_0 = 0.1, lambda_1 = 0, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec21$G

# Shrink toward G_bu
rec22 <- mip_l0(fc2, S, W_I, G_bench = G_bu, 
               lambda_0 = 0.1, lambda_1 = 0.2, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec22$G

rec23 <- mip_l0(fc2, S, W_I, G_bench = G_bu, 
               lambda_0 = 0.1, lambda_1 = 0, lambda_2 = 0.2, M = NULL, 
               solver = "gurobi")
rec23$G

# Shrink toward G_td
rec24 <- mip_l0(fc2, S, W_I, G_bench = G_td, 
               lambda_0 = 0.1, lambda_1 = 0.1, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec24$G

rec25 <- mip_l0(fc2, S, W_I, G_bench = G_td, 
               lambda_0 = 1, lambda_1 = 1, lambda_2 = 1, M = NULL, 
               solver = "gurobi")
rec25$G

rec26 <- mip_l0(fc2, S, W_I, G_bench = G_td, 
               lambda_0 = 1, lambda_1 = 10, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec26$G

