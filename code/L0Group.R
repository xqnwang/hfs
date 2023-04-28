################################################
#   minimize   1/n ||y-xb||_{2}^{2} +
#   (b,z,s)    \lambda_{0} \sum_{g=1}^{q}z_g +
#              \lambda_{2} \sum_{g=1}^{q}s_g
# subject to  ||b_g||_{2}^{2} <= s_g z_g
#             ||b_g||_{1} <= M z_g
#             z_g \in {0, 1}
#             s_g >= 0
#             g \in [q]
#
################################################
#   minimize   1/n v'v +
# (b,z,s,v,u)  \lambda_{0} \sum_{g=1}^{q}z_g +
#              \lambda_{2} \sum_{g=1}^{q}s_g
# subject to  y - xb = v
#             ||b_g||_{2}^{2} <= s_g z_g
#             sum(u_g) <= M z_g
#             u_j >= b_j
#             u_j >= -b_j
#             z_g \in {0, 1}
#             s_g >= 0
#             g \in [q]
#             j \in [p]

library(ROI)

dbind <- function(...) {
  ## sparse matrices construction
  .dbind <- function(x, y) {
    A <- slam::simple_triplet_zero_matrix(NROW(x), NCOL(y))
    B <- slam::simple_triplet_zero_matrix(NROW(y), NCOL(x))
    rbind(cbind(x, A), cbind(B, y))
  }
  Reduce(.dbind, list(...))
}

mip_l0l2 <- function(x, y, group_indices, lambda_0, lambda_2, M,
                     solver = "gurobi"){
  # x: n*p
  # y: n
  # group_indices: list of length q (q groups)
  # lambda_0: lambda_{0}
  # lambda_2: lambda_{2}
  # M: Big-M vaule
  
  stzm <- slam::simple_triplet_zero_matrix
  stdm <- slam::simple_triplet_diag_matrix
  stm <- slam::simple_triplet_matrix
  
  n <- NROW(x); p <- NCOL(x)
  q <- length(group_indices)
  
  # varnames <- c(paste0("b", seq.int(p)),
  #               paste0("z", seq.int(q)),
  #               paste0("s", seq.int(q)),
  #               paste0("v", seq.int(n)),
  #               paste0("u", seq.int(p)))
  
  # Obj
  Q0 <- dbind(stzm(p), stzm(q), stzm(q), stdm(2/n, n), stzm(p))
  a0 <- stm(i = rep(1, 2*q),
            j = (p+1):(p+2*q),
            v = c(rep(lambda_0, q), rep(lambda_2, q)),
            nrow = 1,
            ncol = p+2*q+n+p)
  model <- OP(objective = Q_objective(Q = Q0, L = a0))
  
  # C1: xb + v = y
  A1 <- cbind(x, stzm(n, 2*q), stdm(1, n), stzm(n, p))
  LC1 <- L_constraint(A1, eq(n), y)
  
  # C2: ||b_g||_{2}^{2} - s_g z_g <= 0
  Q2 <- lapply(1:q, function(g){
    x_indices <- group_indices[[g]]
    Q_q <- stm(i = x_indices,
               j = x_indices,
               v = rep(2, length(x_indices)),
               nrow = p+2*q+n+p,
               ncol = p+2*q+n+p)
    Q_q[p+g, p+q+g] <- -2
    Q_q
  })
  L2 <- stzm(nrow = q, ncol = p+2*q+n+p)
  QC2 <- Q_constraint(Q = Q2,
                      L = L2,
                      dir = rep("<=", q),
                      rhs = rep(0, q))
  
  # C3: sum(u_g) - M z_g <= 0
  A3 <- sapply(1:q, function(g){
    x_indices <- group_indices[[g]]
    L_q <- stm(i = rep(1, length(x_indices)),
               j = p+2*q+n+x_indices,
               v = rep(1, length(x_indices)),
               nrow = 1,
               ncol = p+2*q+n+p)
    L_q[1, p+g] <- -M
    as.matrix(L_q)
  }) %>% t()
  LC3 <- L_constraint(A3, rep("<=", q), rep(0, q))
  
  # C4: b_j - u_j <= 0
  A4 <- stm(i = c(1:p, 1:p),
            j = c(1:p, p+2*q+n+(1:p)),
            v = c(rep(1, p), rep(-1, p)),
            nrow = p,
            ncol = p+2*q+n+p)
  LC4 <- L_constraint(A4, rep("<=", p), rep(0, p))
  
  # C5: b_j + u_j >= 0
  A5 <- stm(i = c(1:p, 1:p),
            j = c(1:p, p+2*q+n+(1:p)),
            v = c(rep(1, p), rep(1, p)),
            nrow = p,
            ncol = p+2*q+n+p)
  LC5 <- L_constraint(A5, rep(">=", p), rep(0, p))
  
  constraints(model) <- rbind(LC1, QC2, LC3, LC4, LC5)
  # z_g \in {0, 1}; s_g >= 0
  types(model) <- c(rep("C", p), rep("B", q), rep("C", q+n+p))
  bounds(model) <- V_bound(li = c(1:p, (p+2*q+1):(p+2*q+n)), lb = rep.int(-Inf, p+n), nobj = p+2*q+n+p) # default of lower bound is 0
  
  # Optimal solution
  model.solver <- ROI_solve(model, solver)
  
  # Output
  model.slt <- solution(model.solver)
  b <- model.slt[seq.int(p)] %>% 
    {sapply(1:q, function(g) .[group_indices[[g]]])} %>%
    t() %>%
    `rownames<-`(paste0("group", seq.int(q)))
  z <- model.slt[(p+1):(p+q)]
  return(list(model = model, opt = model.solver, solution = model.slt,
              coefficients = b, binary = z))
}


#----------------------------------------------------------------------
# Comparison with L0Group https://github.com/hazimehh/L0Group/blob/main/Demo.ipynb
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