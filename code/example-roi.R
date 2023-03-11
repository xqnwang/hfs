Sys.setenv(ROI_LOAD_PLUGINS = "FALSE")
library(ROI)
library(slam)
library(magrittr)
library(ROI.plugin.neos)
library(ROI.plugin.gurobi)

dbind <- function(...) {
  ## sparse matrices construction
  .dbind <- function(x, y) {
    A <- simple_triplet_zero_matrix(NROW(x), NCOL(y))
    B <- simple_triplet_zero_matrix(NROW(y), NCOL(x))
    rbind(cbind(x, A), cbind(B, y))
  }
  Reduce(.dbind, list(...))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%
# Example data
#%%%%%%%%%%%%%%%%%%%%%%%%%%
set.seed(123)
N <- 100
S <- rbind(matrix(c(rep(1, 5), c(rep(1, 3), rep(0, 2)), c(rep(0, 3), rep(1, 2))), 
                  nrow = 3, byrow = TRUE),
           diag(1, nrow = 5))
B <- matrix(c(rnorm(N, 1, 1), 
              rnorm(N, 4, 1), 
              rnorm(N, 0, 1),
              rnorm(N, 2, 1),
              rnorm(N, 5, 1)), 
            byrow = FALSE, nrow = N)
Y <- B %*% t(S)
Y_hat <- matrix(c(rnorm(N, 10, 1), 
                  rnorm(N, 6, 1), 
                  rnorm(N, 5, 1), 
                  rnorm(N, 1, 1), 
                  rnorm(N, 4, 1), 
                  rnorm(N, 0, 1),
                  rnorm(N, 2, 1),
                  rnorm(N, 5, 1)), 
                byrow = FALSE, nrow = N)
VecY <- as.vector(Y)
D <- kronecker(S, Y_hat)


#%%%%%%%%%%%%%%%%%%%%%%%%%%
# MIQCP - Mixed Integer Quadratically Constrained Program
#%%%%%%%%%%%%%%%%%%%%%%%%%%
# If a MIP has any constraints containing a quadratic term, regardless of the objective function, the problem is termed a MIQCP.
gp_op <- function(x, y, n, m, lambda, M) {
  ## x: kronecker(S, Y_hat)
  ## y: vec(Y)
  ## lambda: Lagrange multiplier
  ## n: number of all series
  ## m: number of bottom-level series
  ## M: big-M
  stzm <- simple_triplet_zero_matrix
  stdm <- simple_triplet_diag_matrix
  
  Nn <- NROW(x); mn <- NCOL(x)
  Q0 <- dbind(stzm(mn), stdm(2/Nn, Nn), stzm(n))
  a0 <- c(g = double(mn), ga = double(Nn), z = lambda * rep(1, n))
  op <- OP(objective = Q_objective(Q = Q0, L = a0))
  
  ## y - X %*% g = gamma  <=>  X %*% g + gamma = y
  A1 <- cbind(x, stdm(1, Nn), stzm(Nn, n))
  LC1 <- L_constraint(A1, eq(Nn), y)
  ## \sum_{i=0}^{m-1} (g_{j+in})^2 - M z_j <= 0 for j = 1,...,n
  QNULL <- diag(0, mn + Nn + n)
  LNULL <- c(double(mn), double(Nn), double(n))
  Q <- lapply(1:n, function(j){
    i <- 0:(m - 1)
    Q_j <- QNULL
    diag(Q_j)[j + i*n] <- 2
    Q_j
  })
  L <- sapply(1:n, function(j){
    L_j <- LNULL
    L_j[mn + Nn + j] <- -M
    L_j
  }) %>% t()
  QC1 <- Q_constraint(Q = Q,
                      L = L,
                      dir = rep("<=", n),
                      rhs = rep(0, n))
  
  constraints(op) <- rbind(LC1, QC1)
  bounds(op) <- V_bound(li = 1:(mn + Nn), lb = rep.int(-Inf, mn + Nn), nobj = mn + Nn + n)
  types(op) <- c(rep("C", mn + Nn), rep("B", n))
  op
}
op <- gp_op(x = D, y = VecY, n = NROW(S), m = NCOL(S), lambda = 0.1, M = 100)

# Optimal solution - solver = "neos"
job_neos <- ROI_solve(op, "neos", email = "xiaoqian.wang@monash.edu")
## str(job_neos)
slt_neos <- solution(job_neos)
(z_neos <- tail(slt_neos, NROW(S)))
(G_neos <- matrix(slt_neos[1:(NCOL(S)*NROW(S))], 
                  nrow = NROW(S), ncol = NCOL(S), byrow = FALSE) %>% 
    t() %>% 
    round(digits = 3))

# Optimal solution - solver = "gurobi"
job_gurobi <- ROI_solve(op, "gurobi")
## Register a Gurobi account as an academic user, request for a license, and download the current version of Gurobi optimizer
## str(job_gurobi)
slt_gurobi <- solution(job_gurobi)
(z_gurobi <- tail(slt_gurobi, NROW(S)))
(G_gurobi <- matrix(slt_gurobi[1:(NCOL(S)*NROW(S))], 
                    nrow = NROW(S), ncol = NCOL(S), byrow = FALSE) %>% 
    t() %>% 
    round(digits = 3))


#%%%%%%%%%%%%%%%%%%%%%%%%%%
# MIQP - Mixed Integer Quadratic Program
#%%%%%%%%%%%%%%%%%%%%%%%%%%
# If a MIP is a quadratic term in the objective function and all the constraints in the model are linear, the problem is termed a MIQP.
gp_abs_op <- function(x, y, n, m, lambda, M) {
  ## x: kronecker(S, Y_hat)
  ## y: vec(Y)
  ## lambda: Lagrange multiplier
  ## n: number of all series
  ## m: number of bottom-level series
  ## M: big-M
  stzm <- simple_triplet_zero_matrix
  stdm <- simple_triplet_diag_matrix
  
  Nn <- NROW(x); mn <- NCOL(x)
  Q0 <- dbind(stzm(mn), stzm(mn), stdm(2/Nn, Nn), stzm(n))
  A0 <- c(g = double(mn), p = double(mn), ga = double(Nn), z = lambda * rep(1, n))
  op <- OP(objective = Q_objective(Q = Q0, L = A0))
  
  ## y - X %*% g = gamma  <=>  X %*% g + gamma = y
  A1 <- cbind(x, stzm(Nn, mn), stdm(1, Nn), stzm(Nn, n))
  LC1 <- L_constraint(A1, eq(Nn), y)
  ## \sum_{i=0}^{m-1} p_{j+in} - M z_j <= 0 for j = 1,...,n
  LNULL <- c(double(mn), double(mn), double(Nn), double(n))
  A2 <- sapply(1:n, function(j){
    i <- 0:(m - 1)
    L_j <- LNULL
    L_j[mn + j + i*n] <- 1
    L_j[mn + mn + Nn + j] <- -M
    L_j
  }) %>% t()
  LC2 <- L_constraint(A2, rep("<=", n), rep(0, n))
  ## p_k >= g_k for k = 1,...,mn  <=>  g_k - p_k <= 0
  A3 <- sapply(1:mn, function(k){
    L_k <- LNULL
    L_k[k] <- 1
    L_k[mn + k] <- -1
    L_k
  }) %>% t()
  LC3 <- L_constraint(A3, rep("<=", mn), rep(0, mn))
  ## p_k >= -g_k for k = 1,...,mn  <=>  g_k + p_k >= 0
  A4 <- sapply(1:mn, function(k){
    L_k <- LNULL
    L_k[k] <- 1
    L_k[mn + k] <- 1
    L_k
  }) %>% t()
  LC4 <- L_constraint(A4, rep(">=", mn), rep(0, mn))
  
  constraints(op) <- rbind(LC1, LC2, LC3, LC4)
  bounds(op) <- V_bound(li = 1:(2*mn + Nn), lb = rep.int(-Inf, 2*mn + Nn), nobj = mn + mn + Nn + n)
  types(op) <- c(rep("C", mn + mn + Nn), rep("B", n))
  op
}
op_abs <- gp_abs_op(x = D, y = VecY, n = NROW(S), m = NCOL(S), lambda = 0.1, M = 100)

# Optimal solution - solver = "neos"
job_abs_neos <- ROI_solve(op_abs, "neos", email = "xiaoqian.wang@monash.edu")
## str(job_neos)
slt_abs_neos <- solution(job_abs_neos)
(z_abs_neos <- tail(slt_abs_neos, NROW(S)))
(G_abs_neos <- matrix(slt_abs_neos[1:(NCOL(S)*NROW(S))], 
                      nrow = NROW(S), ncol = NCOL(S), byrow = FALSE) %>% 
    t() %>% 
    round(digits = 3))

# Optimal solution - solver = "gurobi"
job_abs_gurobi <- ROI_solve(op_abs, "gurobi")
## Register a Gurobi account as an academic user, request for a license, and download the current version of Gurobi optimizer
## str(job_gurobi)
slt_abs_gurobi <- solution(job_abs_gurobi)
(z_abs_gurobi <- tail(slt_abs_gurobi, NROW(S)))
(G_abs_gurobi <- matrix(slt_abs_gurobi[1:(NCOL(S)*NROW(S))], 
                        nrow = NROW(S), ncol = NCOL(S), byrow = FALSE) %>% 
    t() %>% 
    round(digits = 3))


