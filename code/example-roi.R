Sys.setenv(ROI_LOAD_PLUGINS = "FALSE")
library(ROI)
library(slam)
library(magrittr)
library(ROI.plugin.neos)

dbind <- function(...) {
  ## sparse matrices construction
  .dbind <- function(x, y) {
    A <- simple_triplet_zero_matrix(NROW(x), NCOL(y))
    B <- simple_triplet_zero_matrix(NROW(y), NCOL(x))
    rbind(cbind(x, A), cbind(B, y))
  }
  Reduce(.dbind, list(...))
}

# Example data
set.seed(123)
S <- rbind(matrix(c(rep(1, 5), c(rep(1, 3), rep(0, 2)), c(rep(0, 3), rep(1, 2))), 
                  nrow = 3, byrow = TRUE),
           diag(1, nrow = 5))
B <- matrix(c(rnorm(10, 1, 1), 
              rnorm(10, 4, 1), 
              rnorm(10, 0, 1),
              rnorm(10, 2, 1),
              rnorm(10, 5, 1)), 
            byrow = FALSE, nrow = 10)
Y <- B %*% t(S)
Y_hat <- matrix(c(rnorm(10, 10, 1), 
                  rnorm(10, 6, 1), 
                  rnorm(10, 5, 1), 
                  rnorm(10, 1, 1), 
                  rnorm(10, 4, 1), 
                  rnorm(10, 0, 1),
                  rnorm(10, 2, 1),
                  rnorm(10, 5, 1)), 
                byrow = FALSE, nrow = 10)
VecY <- as.vector(Y)
D <- kronecker(S, Y_hat)

# Quadratic optimization problem
gp_op <- function(x, y, lambda, n, m, M) {
  ## x: kronecker(S, Y_hat)
  ## y: vec(Y)
  ## lambda: Lagrange multiplier
  ## n: number of all series
  ## m: number of bottom-level series
  ## M: big-M
  stzm <- simple_triplet_zero_matrix
  stdm <- simple_triplet_diag_matrix
  
  Nn <- NROW(x); mn <- NCOL(x)
  Q0 <- dbind(stzm(mn), stdm(1, Nn), stzm(n))
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

# Optimal solution
op <- gp_op(x = D, y = VecY, lambda = 3, n = NROW(S), m = NCOL(S), M = 100)
neos_job <- ROI_solve(op, "neos", email = "xiaoqian.wang@monash.edu")
 ## Error in raise_licence_error(job$password) : 
 ##   Error: A valid email address is required for all NEOS jobs In some cases the solver licence does not permit connections via XML-RPC. Therefore these solvers can not be accessed by 'ROI.plugin.neos' directly. An alternative option is to write the problem out via write.op(model, 'my_op.gms', 'gams') and commit it via the web-interface. Or just use a alternative solver.
str(neos_job)
(slt <- solution(neos_job))
(G <- matrix(slt[1:(NCOL(S)*NROW(S))], nrow = NROW(S), ncol = NCOL(S), byrow = FALSE) %>% 
  t() %>% 
  round(digits = 5))


