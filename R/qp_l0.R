#' G Estimation for Forecast Reconciliation with Zero-column Constraints
#' 
#' Estimate G under a trace minimization framework with zero columns specified.
#' 
#' @param fc A vector contains the base forecasts.
#' @param S Summing or structural matrix.
#' @param W Estimated variance-covariance matrix of the h-step-ahead base forecast errors.
#' @param z A vector contains values of integer variables returned by MIP.
#' @param solver A character vector specifying the solver to use. If missing, then `gurobi` is used.
#' 
#' @import ROI
#' 
#' @export
library(ROI)
qp_l0 <- function(fc, S, W, z, solver = "gurobi"){
  # Check input dimension
  if(length(unique(c(length(fc), NROW(S), dim(W), length(z)))) != 1L){
    stop("The dimensions of the inputs do not match each other")
  }
  
  # Dimension info
  n <- NROW(S); n_b <- NCOL(S)
  
  # Parameters in the QP problem
  # parameters <- c(paste0("g", seq.int(n_b*n)), # vec(G)
  #                 paste0("e_check", seq.int(n)))
  n_parameters <- n_b*n + n
  p <- n_b*n
  zero_index <- which(z == 0) |> purrr::map(\(x) (x-1)*n_b + 1:n_b) %>% do.call(c, .)
  n_zero <- length(zero_index)
  solveW <- solve(W)
  
  # Optimization problem construction
  ## Obj: 1/2 t(e_check) W^{-1} e_check
  Q0 <- Matrix::bdiag(matrix(0, p, p), solveW) |> as.matrix()
  model <- ROI::OP(objective = ROI::Q_objective(Q = Q0))
  
  ## Constraints
  ### C1: fc - kronecker(t(fc), S) %*% g = e_check
  ### <==> kronecker(t(fc), S) %*% g + e_check = fc
  A1 <- cbind(kronecker(t(fc), S), diag(n))
  LC1 <- ROI::L_constraint(A1, ROI::eq(n), fc)
  
  ### C2: (unbiasedness constraints) kronecker(t(S), diag(1, n_b)) %*% g = vec(diag(1, n_b))
  A2 <- cbind(kronecker(t(S), diag(n_b)), matrix(0, n_b*n_b, n))
  LC2 <- ROI::L_constraint(A2, ROI::eq(n_b*n_b), as.vector(diag(n_b)))
  
  # constraints(model) <- LC1
  constraints(model) <- rbind(LC1, LC2)
  
  ### z_g \in {0, 1}
  types(model) <- rep("C", n_parameters)
  bounds(model) <- ROI::V_bound(li = c((1:p)[-zero_index], (p + 1):(p + n)),
                                lb = rep(-Inf, (p - n_zero) + n), 
                                ui = zero_index,
                                ub = rep.int(0, n_zero),
                                nobj = n_parameters) # default of lower bound is 0
  
  # optimal solution
  model.solver <- ROI::ROI_solve(model, solver)
  
  # output
  model.slt <- ROI::solution(model.solver)
  G <- model.slt[seq.int(p)] |>
    matrix(nrow = n_b, ncol = n, byrow = FALSE)
  
  return(list(model = model, solver = model.solver, solution = model.slt,
              G = G, z = z))
}