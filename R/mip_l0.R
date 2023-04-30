#' Subset Selection for Forecast Reconciliation
#' 
#' Estimate G under a trace minimization framework, while zeroing out some columns 
#' of G to enable subset selection in the context of forecast reconciliation.
#' 
#' @param fc A vector contains the base forecasts.
#' @param S Summing or structural matrix.
#' @param W Estimated variance-covariance matrix of the h-step-ahead base forecast errors.
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
mip_l0 <- function(fc, S, W, G_bench = NULL, 
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
    M <- ifelse(is.null(G_bench), n_b*10, abs(G_bench) |> colSums() |> max()*10)
  }
  
  # G_bench
  if(is.null(G_bench)){
    G_bench <- matrix(0, nrow = n_b, ncol = n) # G shrink towards zero matrix
  }else if(NROW(G_bench) != n_b | NCOL(G_bench) != n){
    stop("The dimensions of G_bench do not match that of S")
  }
  
  # Parameters in the MIP problem
  # parameters <- c(paste0("g", seq.int(n_b*n)), # vec(G)
  #                 paste0("z", seq.int(n)),
  #                 paste0("e_check", seq.int(n)),
  #                 paste0("d_positive", seq.int(n_b*n)),
  #                 paste0("g_positive", seq.int(n_b*n)))
  n_parameters <- n_b*n + n + n + n_b*n + n_b*n
  p <- n_b*n
  solveW <- solve(W)
  
  # Optimization problem construction
  ## Obj: 1/2 t(e_check) W^{-1} e_check + lambda_0 \sum_{j=1}^{n} z_j + 
  ##                                    lambda_1 \sum d_positive +
  ##                                    lambda_2 t(d_positive) d_positive
  Q0 <- Matrix::bdiag(matrix(0, p, p), matrix(0, n, n), solveW,
                      diag(x = 2*lambda_2, p), matrix(0, p, p)) |> 
    as.matrix()
  A0 <- Matrix::sparseVector(
    i = c((p + 1):(p + n), (p + n + n + 1):(p + n + n + p)),
    x = c(rep(lambda_0, n), rep(lambda_1, p)),
    length = n_parameters
  ) |> as.vector()
  model <- ROI::OP(objective = ROI::Q_objective(Q = Q0, L = A0))
  
  ## Constraints
  ### C1: fc - kronecker(t(fc), S) %*% g = e_check
  ### <==> kronecker(t(fc), S) %*% g + e_check = fc
  A1 <- cbind(kronecker(t(fc), S), matrix(0, n, n),
              diag(n), matrix(0, n, 2*p))
  LC1 <- ROI::L_constraint(A1, ROI::eq(n), fc)
  
  ### C2: kronecker(t(S), diag(1, n_b)) %*% g = vec(diag(1, n_b))
  A2 <- cbind(kronecker(t(S), diag(n_b)), matrix(0, n_b*n_b, 2*n + 2*p))
  LC2 <- ROI::L_constraint(A2, ROI::eq(n_b*n_b), as.vector(diag(n_b)))
  
  ### C3: \sum g_positive_{(1+(j-1)*n_b):(j*n_b)} - M z_j <= 0 for j = 1,...,n
  A3 <- sapply(1:n, function(j){
    L_j <- Matrix::sparseVector(
      i = c(p + j, p + n + n + p + (1+(j-1)*n_b):(j*n_b)),
      x = c(-M, rep(1, n_b)),
      length = n_parameters
    )
    return(as.matrix(L_j))
  }) |> t()
  LC3 <- ROI::L_constraint(A3, rep("<=", n), rep(0, n))
  
  ### C4: g - g_positive <= 0
  A4 <- cbind(diag(p), matrix(0, p, 2*n + p), diag(x = -1, p))
  LC4 <- ROI::L_constraint(A4, rep("<=", p), rep(0, p))
  
  ### C5: g + g_positive >= 0
  A5 <- cbind(diag(p), matrix(0, p, 2*n + p), diag(p))
  LC5 <- ROI::L_constraint(A5, rep(">=", p), rep(0, p))
  
  ### C6: g - d_positive <= vec(G_bench)
  A6 <- cbind(diag(p), matrix(0, p, 2*n), diag(x = -1, p), matrix(0, p, p))
  LC6 <- ROI::L_constraint(A6, rep("<=", p), as.vector(G_bench))
  
  ### C7: g + d_positive >= vec(G_bench)
  A7 <- cbind(diag(p), matrix(0, p, 2*n), diag(p), matrix(0, p, p))
  LC7 <- ROI::L_constraint(A7, rep(">=", p), as.vector(G_bench))
  
  if(lambda_0 == 0L & lambda_1 == 0L & lambda_2 == 0L){
    constraints(model) <- rbind(LC1, LC2)
  }else if(lambda_0 == 0L){
    constraints(model) <- rbind(LC1, LC2, LC6, LC7)
  }else if(lambda_1 == 0L & lambda_2 == 0L){
    constraints(model) <- rbind(LC1, LC2, LC3, LC4, LC5)
  } else {
    constraints(model) <- rbind(LC1, LC2, LC3, LC4, LC5, LC6, LC7)
  }

  ### z_g \in {0, 1}
  types(model) <- c(rep("C", p), rep("B", n), rep("C", n + 2*p))
  bounds(model) <- ROI::V_bound(li = c(1:p, (p + n + 1):(p + n + n)), lb = rep.int(-Inf, p + n), nobj = n_parameters) # default of lower bound is 0
  
  # optimal solution
  model.solver <- ROI::ROI_solve(model, solver)
  
  # output
  model.slt <- ROI::solution(model.solver)
  G <- model.slt[seq.int(p)] |>
    matrix(nrow = n_b, ncol = n, byrow = FALSE)
  if(lambda_0 == 0L){
    z <- NA
  }else{
    z <- model.slt[(p + 1):(p + n)]
  }
  
  return(list(model = model, solver = model.solver, solution = model.slt,
              G = G, z = z))
}