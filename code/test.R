#----------------------------------------------------------------------
# Comparison with L0Group
#----------------------------------------------------------------------
# import numpy as np
# import scipy as sc
# import math
# from bnb import gen_synthetic
# 
# 
# # Set the data generation parameters.
# rho = 0.3 # Correlation parameter.
# n = 1000 # Number of observations.
# p = 100 # Number of features. 
# k0 = 5 # Number of nonzero groups.
# SNR = 10 # Signal-to-noise ratio.
# Gsize = 10 # Group size.
# 
# x, y, group_indices, true_coefs = gen_synthetic(rho,n,p,k0,SNR,Gsize)

library(dplyr)
source('R/sourceDir.R')
sourceDir("R", recursive = TRUE)

x <- read.csv("data/x.csv", header = FALSE) %>% 
  as.matrix()
y <- read.csv("data/y.csv", header = FALSE) %>% 
  {as.vector(.$V1)}
group_indices <- read.csv("data/group_indices.csv", header = FALSE) %>% 
  {.+1} %>%
  {lapply(seq.int(NROW(.)), function(g) as.numeric(.[g,]))}

model <- mip_l0l2(x = x, y = y, group_indices = group_indices, 
                  lambda_0 = 20, lambda_2 = 2, M = 50, 
                  solver = "gurobi")

b <- model$coefficients
z <- model$binary


#----------------------------------------------------------------------
# Subset selection with shrinkage under unbiasedness assumptions
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

fc <- c(10, 4, 4, 2, 4, 2, 2)
S <- rbind(rep(1, 4),
           c(1, 1, 0, 0), c(0, 0, 1, 1),
           diag(1, 4))
W <- diag(7)
G_bu <- cbind(rep(0, 4), rep(0, 4), rep(0, 4), diag(1, 4))
G_td <- cbind(c(0.2, 0.4, 0.2, 0.2), matrix(0, nrow = 4, ncol = 6))

# No shrinkage
rec1 <- mip_l0(fc, S, W, G0 = NULL, 
               lambda_0 = 0.1, lambda_1 = 0, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec1$z
rec1$G

# Shrink toward G_bu
rec2 <- mip_l0(fc, S, W, G0 = G_bu, 
               lambda_0 = 0.1, lambda_1 = 0.1, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec2$z
rec2$G

rec3 <- mip_l0(fc, S, W, G0 = G_bu, 
               lambda_0 = 0.1, lambda_1 = 0, lambda_2 = 0.1, M = NULL, 
               solver = "gurobi")
rec3$z
rec3$G

# Shrink toward G_td
rec4 <- mip_l0(fc, S, W, G0 = G_td, 
               lambda_0 = 0.1, lambda_1 = 0.1, lambda_2 = 0, M = NULL, 
               solver = "gurobi")
rec4$z
rec4$G

rec5 <- mip_l0(fc, S, W, G0 = G_td, 
               lambda_0 = 10, lambda_1 = 20, lambda_2 = 20, M = NULL, 
               solver = "gurobi")
rec5$z
rec5$G



