# Subset selection achieved via adjustment of the MinT estimates

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

mint_adjust <- function(fc, S, W, G, lambda_0, M = NULL, solver = "gurobi"){
  # fc: an n-dimensional vector of base forecasts
  # S: an n*n_b summing matrix
  # W: an n*n covariance matrix of the base forecast errors
  # lambda_0: lambda_{0}
  # M: a Big-M value
  
  # sparse functions
  stzm <- slam::simple_triplet_zero_matrix
  stdm <- slam::simple_triplet_diag_matrix
  stm <- slam::simple_triplet_matrix
  as.matrix <- Matrix::as.matrix
  t <- Matrix::t
  solve <- Matrix::solve

  # check input dimension
  if(length(unique(c(length(fc), NROW(S), NROW(W), NCOL(W)))) != 1L){
    stop("The dimensions of the inputs do not match each other")
  }
  
  # Big-M value
  if(is.null(M)){
    M <- 100
  }
  
  # dimension info
  n <- NROW(S); n_b <- NCOL(S)
  
  # parameters in MIP
  # parameters <- c(paste0("a", seq.int(n*n)), # vec(A)
  #                 paste0("z", seq.int(n)),
  #                 paste0("e_check", seq.int(n)))
  n_parameters <- n*n + n + n
  
  X <- kronecker(t(fc), S %*% G)
  # optimization problem construction
  ## Obj: 1/2 t(e_check) W^{-1} e_check + lambda_0 \sum_{j=1}^{n} z_j
  Q0 <- dbind(stzm(n*n), stzm(n), solve(W))
  A0 <- stm(i = rep(1, n),
            j = c((n*n + 1):(n*n + n)),
            v = rep(lambda_0, n),
            nrow = 1, ncol = n_parameters)
  model <- OP(objective = Q_objective(Q = Q0, L = A0))
  
  ## constraints
  ### C1: fc - X %*% a = e_check
  ### <==> X %*% a + e_check = fc
  A1 <- cbind(X, stzm(n), stdm(1, n))
  LC1 <- L_constraint(A1, eq(n), fc)
  
  ### C2: kronecker(t(S), G) %*% a = vec(diag(1, n_b))
  #A2 <- cbind(kronecker(t(S), G), stzm(n_b*n_b, n), stzm(n_b*n_b, n))
  #LC2 <- L_constraint(A2, eq(n_b*n_b), as.vector(diag(1, n_b)))
  
  ### C3: a_{ij} = 0, i ≠ j
  A3 <- stm(i = seq.int(n*n - n),
            j = seq.int(n*n)[-(1 + (n + 1) * seq.int(0, n-1, 1))],
            v = rep(1, n*n - n),
            nrow = n*n - n,
            ncol = n_parameters)
  LC3 <- L_constraint(A3, eq(n*n - n), rep(0, n*n - n))
  
  ### C4: A1 - M z_j <= 0 for j = 1,...,n
  ### C5: A1 + M z_j >= 0 for j = 1,...,n
  
  ### C4: kronecker(t(rep(1, n)), I_n) a - M z <= 0
  A4 <- cbind(kronecker(t(rep(1, n)), diag(1, n)), stdm(-M, n), stzm(n))
  LC4 <- L_constraint(A4, rep("<=", n), rep(0, n))
  
  ### C5: kronecker(t(rep(1, n)), I_n) a + M z >= 0
  A5 <- cbind(kronecker(t(rep(1, n)), diag(1, n)), stdm(M, n), stzm(n))
  LC5 <- L_constraint(A5, rep(">=", n), rep(0, n))
  
  if(lambda_0 == 0L){
    constraints(model) <- rbind(LC1, LC3)
  } else {
    constraints(model) <- rbind(LC1, LC3, LC4, LC5)
  }
  
  ### z_g \in {0, 1}
  types(model) <- c(rep("C", n*n), rep("B", n), rep("C", n))
  bounds(model) <- V_bound(li = c(1:(n*n), (n*n + n + 1):(n*n + n + n)),
                           lb = rep.int(-Inf, n*n + n),
                           #ui = c(1:(n*n))[-seq.int(1, by = n, length.out = n)],
                           #ub = rep.int(0, n*n - n),
                           nobj = n_parameters) # default of lower bound is 0
  
  # optimal solution
  model.solver <- ROI_solve(model, solver)
  
  # output
  model.slt <- solution(model.solver)
  A <- model.slt[seq.int(n*n)] %>%
    matrix(nrow = n, ncol = n, byrow = FALSE)
  G_adjust <- G %*% A
  z <- model.slt[(n*n+1):(n*n + n)]
  
  return(list(model = model, solver = model.solver, solution = model.slt,
              G = G_adjust, z = z))
}