# least squares optimization problem with a linear equality constraint


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

lsop <- function(fc, S, W, solver = "gurobi"){
  # fc: an n-dimensional vector of base forecasts
  # S: an n*n_b summing matrix
  # W: an n*n covariance matrix of the base forecast errors
  
  # sparse functions
  stzm <- slam::simple_triplet_zero_matrix
  stdm <- slam::simple_triplet_diag_matrix
  stm <- slam::simple_triplet_matrix
  as.matrix <- Matrix::as.matrix
  t <- Matrix::t
  solve <- Matrix::solve
  cov2cor <- Matrix::cov2cor
  
  # check input dimension
  if(length(unique(c(length(fc), NROW(S), NROW(W), NCOL(W)))) != 1L){
    stop("The dimensions of the inputs do not match each other")
  }
  
  # dimension info
  n <- NROW(S); n_b <- NCOL(S)
  
  # parameters in MIP
  # parameters <- c(paste0("y_tilde", seq.int(n)), # vec(G)
  #                 paste0("e_check", seq.int(n)))
  n_parameters <- n + n
  
  # optimization problem construction
  ## Obj: 1/2 t(e_check) W^{-1} e_check
  Q0 <- dbind(stzm(n), solve(W))
  model <- OP(objective = Q_objective(Q = Q0))
  
  ## constraints
  ### C1: fc - y_tilde = e_check
  ### <==> y_tilde + e_check = fc
  A1 <- cbind(stdm(1, n), stdm(1, n))
  LC1 <- L_constraint(A1, eq(n), fc)
  
  ### C2: U' %*% y_tilde = 0
  ### <==> [I_{n-n_b} | -C] %*% y_tilde = 0
  A2 <- cbind(cbind(stdm(1, n - n_b), -S[seq.int(n - n_b), ]), stzm(n - n_b, n))
  LC2 <- L_constraint(A2, eq(n - n_b), rep(0, n - n_b))

  constraints(model) <- rbind(LC1, LC2)
  
  ### z_g \in {0, 1}
  types(model) <- c(rep("C", n), rep("C", n))
  bounds(model) <- V_bound(li = 1:(2*n), lb = rep.int(-Inf, 2*n), nobj = n_parameters) # default of lower bound is 0
  
  # optimal solution
  model.solver <- ROI_solve(model, solver)
  
  # output
  model.slt <- solution(model.solver)
  y_tilde <- model.slt[seq.int(n)]
  
  return(list(model = model, solver = model.solver, solution = model.slt,
              y_tilde = y_tilde))
}