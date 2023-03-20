# Subset selection with shrinkage under unbiasedness assumptions
## Weight matrix G is directly calculated using forecasts in test set
## Cross-validation may be still needed to select the best value of lambda


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

mip_l0 <- function(fc, S, W, G0 = NULL, 
                   lambda_0, lambda_1 = 0, lambda_2 = 0, M = NULL, solver = "gurobi"){
  # fc: an n-dimensional vector of base forecasts
  # S: an n*n_b summing matrix
  # W: an n*n covariance matrix of the base forecast errors
  # G0: a n_b*n benchmark weight matrix used for shrinkage
  # lambda_0: lambda_{0}
  # lambda_1: lambda_{1}
  # lambda_2: lambda_{2}
  # M: a Big-M vaule
  
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
  if(is.null(G0)){
    lambda_1 <- lambda_2 <- 0;
  } else if (!identical(NCOL(G0), NROW(S))){
    stop("The dimensions of the inputs do not match each other")
  } else if (!identical(NROW(G0), NCOL(S))){
    stop("The dimensions of the inputs do not match each other")
  }
  
  # Big-M value
  if(is.null(M)){
    M <- NCOL(S)*3
  }
  
  # dimension info
  n <- NROW(S); n_b <- NCOL(S)
  
  # parameters in MIP
  # parameters <- c(paste0("g", seq.int(n_b*n)), # vec(G)
  #                 paste0("z", seq.int(n)),
  #                 paste0("e_check", seq.int(n)),
  #                 paste0("d_positive", seq.int(n_b*n)),
  #                 paste0("g_positive", seq.int(n_b*n)))
  n_parameters <- n_b*n + n + n + n_b*n + n_b*n
  
  # optimization problem construction
  ## Obj: 1/2 t(e_check) W^{-1} e_check + lambda_0 \sum_{j=1}^{n} z_j + 
  ##                                    lambda_1 \sum d_positive +
  ##                                    lambda_2 t(d_positive) d_positive
  Q0 <- dbind(stzm(n_b*n), stzm(n), solve(W), stdm(2*lambda_2, n_b*n), stzm(n_b*n))
  A0 <- stm(i = rep(1, n + n_b*n),
            j = c((n_b*n + 1):(n_b*n + n), (n_b*n + n + n + 1):(n_b*n + n + n + n_b*n)),
            v = c(rep(lambda_0, n), rep(lambda_1, n_b*n)),
            nrow = 1, ncol = n_parameters)
  model <- OP(objective = Q_objective(Q = Q0, L = A0))
  
  ## constraints
  ### C1: fc - kronecker(t(fc), S) %*% g <= e_check
  ### <==> kronecker(t(fc), S) %*% g + e_check >= fc
  A1 <- cbind(kronecker(t(fc), S), stzm(n), stdm(1, n), stzm(n, n_b*n), stzm(n, n_b*n))
  LC1 <- L_constraint(A1, rep(">=", n), fc)
  
  ### C2: kronecker(t(S), diag(1, n_b)) %*% g = vec(diag(1, n_b))
  A2 <- cbind(kronecker(t(S), diag(1, n_b)), stzm(n_b*n_b, n), stzm(n_b*n_b, n), stzm(n_b*n_b, n_b*n), stzm(n_b*n_b, n_b*n))
  LC2 <- L_constraint(A2, eq(n_b*n_b), as.vector(diag(1, n_b)))
  
  ### C3: \sum g_positive_{(1+(j-1)*n_b):(j*n_b)} - M z_j <= 0 for j = 1,...,n
  A3 <- sapply(1:n, function(j){
    L_j <- stm(i = rep(1, 1 + n_b),
               j = c(n_b*n + j, n_b*n + n + n + n_b*n + (1+(j-1)*n_b):(j*n_b)),
               v = c(-M, rep(1, n_b)),
               nrow = 1, ncol = n_parameters)
    as.matrix(L_j)
  })  %>% t()
  LC3 <- L_constraint(A3, rep("<=", n), rep(0, n))
  
  ### C4: g - g_positive <= 0
  A4 <- cbind(stdm(1, n_b*n), stzm(n_b*n, n), stzm(n_b*n, n), stzm(n_b*n), stdm(-1, n_b*n))
  LC4 <- L_constraint(A4, rep("<=", n_b*n), rep(0, n_b*n))
  
  ### C5: g + g_positive >= 0
  A5 <- cbind(stdm(1, n_b*n), stzm(n_b*n, n), stzm(n_b*n, n), stzm(n_b*n), stdm(1, n_b*n))
  LC5 <- L_constraint(A5, rep(">=", n_b*n), rep(0, n_b*n))
  
  ### C6: g - d_positive <= vec(G0)
  A6 <- cbind(stdm(1, n_b*n), stzm(n_b*n, n), stzm(n_b*n, n), stdm(-1, n_b*n), stzm(n_b*n))
  LC6 <- L_constraint(A6, rep("<=", n_b*n), if(is.null(G0)){rep(0, n_b*n)}else{as.vector(G0)})
  
  ### C7: g + d_positive >= vec(G0)
  A7 <- cbind(stdm(1, n_b*n), stzm(n_b*n, n), stzm(n_b*n, n), stdm(1, n_b*n), stzm(n_b*n))
  LC7 <- L_constraint(A7, rep(">=", n_b*n), if(is.null(G0)){rep(0, n_b*n)}else{as.vector(G0)})
  
  constraints(model) <- rbind(LC1, LC2, LC3, LC4, LC5, LC6, LC7)
  
  ### z_g \in {0, 1}
  types(model) <- c(rep("C", n_b*n), rep("B", n), rep("C", n), 
                    rep("C", n_b*n), rep("C", n_b*n))
  bounds(model) <- V_bound(li = c(1:(n_b*n), (n_b*n + n + 1):(n_b*n + n + n)), lb = rep.int(-Inf, n_b*n + n), nobj = n_parameters) # default of lower bound is 0
  
  # optimal solution
  model.solver <- ROI_solve(model, solver)
  
  # output
  model.slt <- solution(model.solver)
  G <- model.slt[seq.int(n_b*n)] %>%
    matrix(nrow = n_b, ncol = n, byrow = FALSE)
  z <- model.slt[(n_b*n+1):(n_b*n + n)]
  
  return(list(model = model, opt = model.solver, solution = model.slt,
              G = G, z = z))
}
