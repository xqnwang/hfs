#' Subset Selection for Forecast Reconciliation via A diagonal matrix (without unbiasedness assumption)
#' 
#' Estimate diagonal matrix A to zero out some columns of G, thus producing G*A
#' 
#' @param fc A vector contains the base forecasts.
#' @param S Summing or structural matrix.
#' @param W Estimated variance-covariance matrix of the h-step-ahead base forecast errors.
#' @param G Estimated G using MinT methods.
#' @param G_bench Target G matrix shrinks towards.
#' @param lambda_0 A user supplied `lambda_0` value.
#' @param lambda_1 A user supplied `lambda_1` value.
#' @param lambda_2 A user supplied `lambda_2` value.
#' @param M The value of the Big M.
#' @param solver A character vector specifying the solver to use. If missing, then `gurobi` is used.
#' 
#' @import ROI
#' 
#' @export
library(ROI)
mip_l0_diag <- function(fc, S, W, G, G_bench = NULL, unbiased = TRUE,
                        lambda_0 = 0, lambda_1 = 0, lambda_2 = 0, 
                        M = NULL, solver = "gurobi"){
  # Check input dimension
  if(length(unique(c(length(fc), NROW(S), dim(W)))) != 1L){
    stop("The dimensions of the inputs do not match each other")
  }
  
  # Dimension info
  n <- NROW(S); n_b <- NCOL(S)
  
  # Set Big-M value
  if(is.null(M)){
    M <- n_b*10
  }
  
  # Parameters in the MIP problem
  # parameters <- c(paste0("a", seq.int(n)), # diag(A)
  #                 paste0("z", seq.int(n)),
  #                 paste0("e_check", seq.int(n)))
  n_parameters <- n + n + n
  diag_index <- (1 + (n + 1) * seq.int(0, n-1, 1))
  SG <- S %*% G
  solveW <- solve(W)
  X <- kronecker(t(fc), SG)[, diag_index]
  
  # Optimization problem construction
  ## Obj: 1/2 t(e_check) W^{-1} e_check + lambda_0 \sum_{j=1}^{n} z_j + 
  ##                                    lambda_2 t(a) %*% a
  Q0 <- Matrix::bdiag(diag(x = 2*lambda_2, nrow = n),
                      matrix(0, n, n),
                      solveW) |>
    as.matrix()
  A0 <- Matrix::sparseVector(
    i = (n + 1):(n + n),
    x = rep(lambda_0, n),
    length = n_parameters
  ) |> as.vector()
  model <- ROI::OP(objective = ROI::Q_objective(Q = Q0, L = A0))
  
  ## Constraints
  ### C1: fc - X %*% a = e_check
  ### <==> X %*% a + e_check = fc
  A1 <- cbind(X, matrix(0, n, n), diag(x = 1, nrow = n))
  LC1 <- ROI::L_constraint(A1, ROI::eq(n), fc)
  
  ### C2: a_{ij} = 0, i ≠ j
  # A2 <- Matrix::sparseMatrix(
  #   i = seq.int(n*n - n),
  #   j = seq.int(n*n)[-(1 + (n + 1) * seq.int(0, n-1, 1))],
  #   x = rep(1, n*n - n),
  #   dims = c(n*n - n, n_parameters)
  # ) |> as.matrix()
  # LC2 <- ROI::L_constraint(A2, ROI::eq(n*n - n), rep(0, n*n - n))
  
  ### C3: A1 - Mz <= 0
  A3 <- cbind(diag(x = 1, nrow = n),
              diag(x = -M, nrow = n),
              matrix(0, n, n))
  LC3 <- L_constraint(A3, rep("<=", n), rep(0, n))
  
  ### C4: A1 + M z >= 0
  A4 <- cbind(diag(x = 1, nrow = n),
              diag(x = M, nrow = n),
              matrix(0, n, n))
  LC4 <- L_constraint(A4, rep(">=", n), rep(0, n))
  
  ### C5: \sum_{i=1}^{n} zj >= n_b
  A5 <- Matrix::sparseVector(
    i = (n + 1):(n + n),
    x = rep(1, n),
    length = n_parameters
  ) |> as.vector()
  LC5 <- L_constraint(A5, rep(">=", 1), rep(n_b, 1))
  
  if (unbiased){
    ### C6: (unbiasedness constraints) kronecker(t(S), G) %*% vec(A) = vec(I_nb)
    A6 <- cbind(kronecker(t(S), G)[, diag_index],
                matrix(0, n_b*n_b, n),
                matrix(0, n_b*n_b, n))
    LC6 <- L_constraint(A6, ROI::eq(n_b*n_b), as.vector(diag(x = 1, nrow = n_b)))
    
    constraints(model) <- rbind(LC1, LC3, LC4, LC5, LC6)
  } else{
    constraints(model) <- rbind(LC1, LC3, LC4, LC5)
  }
  
  ### z_g \in {0, 1}
  types(model) <- c(rep("C", n), rep("B", n), rep("C", n))
  bounds(model) <- ROI::V_bound(li = c(1:n, (n + n + 1):(n + n + n)),
                                lb = c(rep(-M, n), rep(-Inf, n)),
                                ui = 1:n,
                                ub = rep(M, n),
                                nobj = n_parameters) # default of lower bound is 0
  
  # optimal solution
  model.solver <- ROI::ROI_solve(model, solver)
  
  # output
  model.slt <- ROI::solution(model.solver)
  z <- model.slt[(n + 1):(n + n)]
  a <- model.slt[seq.int(n)]
  A <- diag(a, nrow = n)
  G <- G %*% A
  
  return(list(model = model, solver = model.solver, solution = model.slt,
              G = G, z = z, a = a))
}