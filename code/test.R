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
#----------------------------------------------------------------------

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

{b <- model$coefficients}
{z <- model$binary}


