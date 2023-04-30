#' Forecast reconciliation
#' 
#' Reconcile forecasts in accordance with its key structure.
#' 
#' @param base_forecasts A vector or a matrix contains the base forecasts.
#' @param S Summing or structural matrix.
#' @param method Reconciliation method to use.
#' @param residuals Historical residuals used to estimate the variance-covariance matrix of the h-step-ahead base forecast errors. It should not be `NULL` when method is one of `c("wls_var", "mint_cov", "mint_shrink")`.
#' @param fitted_values Fitted values in training set.
#' @param train_data Training data.
#' @param subset Logical. If true, mixed integer programming is implemented to achieve subset selection.
#' @param lasso Logical. If true, `l1 norm` is included to achieve shrinkage.
#' @param ridge Logical. If true, `l2 norm` is included to achieve shrinkage.
#' @param G_bench Target G matrix shrinks towards.
#' @param lambda_0 A user supplied `lambda_0` value.
#' @param lambda_1 A user supplied `lambda_1` value.
#' @param lambda_2 A user supplied `lambda_2` value.
#' @param nlambda_0 Number of candidate `lambda_0` values to choose from.
#' @param M The value of the Big M.
#' @param solver A character vector specifying the solver to use. If missing, then `gurobi` is used.
#' @param parallel Logical. If true, optimal `lambda_0` will be found in parallel.
#' @param workers Number of workers when `parallel = TRUE`.
#' 
#' @import ROI
#' @import future
#' 
#' @export
reconcile <- function(base_forecasts, S, 
                      method = c("bu", "ols", "wls_struct", "wls_var", "mint_cov", "mint_shrink"), 
                      residuals = NULL, 
                      fitted_values = NULL, train_data = NULL,
                      subset = FALSE, lasso = FALSE, ridge = FALSE, G_bench = c("Zero", "G"),
                      lambda_0 = NULL, lambda_1 = 0, lambda_2 = 0, 
                      nlambda_0 = 20, M = NULL, solver = "gurobi",
                      parallel = FALSE, workers = 2){
  # Dimension info
  n <- NROW(S); n_b <- NCOL(S)
  if (is.vector(base_forecasts)){
    if (length(base_forecasts) != n){
      stop("The dimensions of the base_forecasts do not match S")
    }
    base_forecasts <- as.matrix(base_forecasts, nrow = 1, ncol = n)
  } else{
    if(NCOL(base_forecasts) != n){
      stop("The dimensions of the base_forecasts do not match S")
    }
  }
  
  method <- match.arg(method)
  G_bench <- match.arg(method)
  if (method == "bu"){
    # Buttom-up
    subset <- FALSE
    G <- cbind(matrix(0, n_b, n-n_b), diag(n_b))
  } else{
    # W matrix
    if (method == "ols"){
      W <- diag(n)
    } else if (method == "wls_struct"){
      W <- diag(rowSums(S))
    } else{
      # Unbiased sample covariance estimator
      if (!is.null(residuals)){
        Tt <- NROW(residuals)
        covm <- crossprod(stats::na.omit(residuals)) / Tt
      } else{
        if (is.null(fitted_values) | is.null(train_data)){
          stop("Historical residuals are required to estimate the variance-covariance matrix of the h-step-ahead base forecast errors")
        } else{
          residuals <- train_data - fitted_values
          Tt <- NROW(residuals)
          covm <- crossprod(stats::na.omit(residuals)) / Tt
        }
      }
      
      if (method == "wls_var"){
        W <- diag(diag(covm))
      } else if (method == "mint_cov"){
        W <- covm
      } else if (method == "mint_shrink"){
        tar <- diag(apply(residuals, 2, compose(crossprod, stats::na.omit))/Tt)
        corm <- cov2cor(covm)
        xs <- scale(residuals, center = FALSE, scale = sqrt(diag(covm)))
        xs <- xs[stats::complete.cases(xs),]
        v <- (1/(Tt * (Tt - 1))) * (crossprod(xs^2) - 1/Tt * (crossprod(xs))^2)
        diag(v) <- 0
        corapn <- cov2cor(tar)
        d <- (corm - corapn)^2
        lambda <- sum(v)/sum(d)
        lambda <- max(min(lambda, 1), 0)
        W <- lambda * tar + (1 - lambda) * covm
      } else {
        stop("Unknown reconciliation method")
      }
    }
    
    # Check positive definiteness of weights
    eigenvalues <- eigen(W, only.values = TRUE)[["values"]]
    if (any(eigenvalues < 1e-8)) {
      stop("MinT methods need covariance matrix to be positive definite.", call. = FALSE)
    }
    
    # G matrix
    R <- t(S) %*% solve(W)
    G <- solve(R %*% S) %*% R
    
    # Subset selection
    if (subset){
      # One-step ahead base forecasts
      fc <- base_forecasts[1, ] |> as.vector()
      
      # Candidate lambda_0
      if (is.null(lambda_0)){
        lambda_0_max <- 0.5 * (t(fc) %*% solve(W) %*% fc)/n_b
        lambda_0 <- c(0, 
                      exp(seq(from = log(1e-04*lambda_0_max), 
                              to = log(lambda_0_max), 
                              by = log(1e04)/(nlambda_0 - 2)))
                      )
      }
      if (lambda_1 == 0L & lambda_2 == 0L){
        lambda_1 <- 1e-5
      }
      
      # Shrinkage matrix
      if (G_bench == "Zero"){
        G_bench <- NULL
      } else{
        G_bench <- G
      }
      
      # Find optimal lambda_0 by minimizing sum of squared reconciled residuals
      if (length(lambda_0) > 1){
        if (is.null(train_data) | is.null(fitted_values)){
          stop("Training data and fitted values are required to find the optimal lambda_0")
        }
        
        if (parallel){
          future::plan(multisession, workers = workers)
          map_fun <- furrr::future_map
        } else {
          map_fun <- purrr::map
        }
        
        fit.lambda <- lambda_0[-1] |>
          map_fun(\(l0) mip_l0(lambda_0 = l0,
                               fc = fc, S = S, W = W, G_bench = G_bench,
                               lambda_1 = lambda_1, lambda_2 = lambda_2,
                               M = M, solver = solver)$G)
        fit.lambda <- append(fit.lambda, list(G), after = 0) # lambda_0 = 0
        sse <- purrr::map_dbl(fit.lambda, \(x) sum(stats::na.omit(train_data - fitted_values %*% t(x) %*% t(S))^2)) |> round(2)
        sse_summary <- data.frame(lambda0 = lambda_0, sse = sse)
        lambda_0 <- lambda_0[which.min(sse)]
      }
      if (lambda_0 == 0L){
        G <- G
        z <- NA
      } else{
        fit.mip <- mip_l0(fc = fc, S = S, W = W, G_bench = G_bench, 
                          lambda_0 = lambda_0, lambda_1 = lambda_1, lambda_2 = lambda_2, 
                          M = M, solver = solver)
        
        G <- fit.mip$G
        z <- fit.mip$z
      }
    }
  }
  
  # Reconciliation
  y_tilde <- base_forecasts %*% t(G) %*% t(S)
  colnames(y_tilde) <- colnames(base_forecasts)
  
  # Output
  if (subset){
    z <- z
  } else{
    z <- NA
  }
  if (exists("sse_summary")){
    lambda0_report <- sse_summary
  } else{
    lambda0_report <- NA
  }
  list(y_tilde = y_tilde, G = G,
       z = z,
       lambda0_report = lambda0_report)
}
