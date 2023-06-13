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
#' @param lasso Logical. If true, Group lasso is included to achieve shrinkage.
#' @param ridge Logical. If true, `l2 norm` is included to achieve shrinkage.
#' @param G_bench Target G matrix shrinks towards.
#' @param ELasso Logical. If both lasso and ELasso are TRUE, Group lasso based on training set is implemented.
#' @param lambda_0 A user supplied `lambda_0` value for subset selection.
#' @param lambda_1 A user supplied `lambda_1` value for group lasso.
#' @param lambda_2 A user supplied `lambda_2` value for subset selection.
#' @param nlambda Number of candidate `lambda_0`(or `lambda_1`) values to choose from for subset selection (group lasso) problem.
#' @param M The value of the Big M.
#' @param solver A character vector specifying the solver to use. If missing, then `gurobi` is used.
#' @param parallel Logical. If true, optimal `lambda_0` will be found in parallel.
#' @param workers Number of workers when `parallel = TRUE`.
#' @param .progress Logical. If true, a progress bar will be displayed when searching optimal `lambda_0`.
#' 
#' @import ROI
#' @import future
#' @import gglasso
#' @import reticulate
#' 
#' @export
.reconcile <- function(base_forecasts, S, 
                      method = c("bu", "ols", "wls_struct", "wls_var", "mint_cov", "mint_shrink"), 
                      residuals = NULL, 
                      fitted_values = NULL, train_data = NULL,
                      lasso = FALSE, ELasso = FALSE, lambda_1 = NULL,
                      subset = FALSE, ridge = FALSE, lambda_0 = NULL, lambda_2 = 0, 
                      G_bench = c("Zero", "G"), unbiased = TRUE,
                      nlambda = 20, M = NULL, solver = "gurobi", 
                      pythonpath = NULL,
                      pyfunpath = NULL,
                      parallel = FALSE, workers = 2, .progress = FALSE){
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
        lambda_0_max <- (0.5 * (t(fc) %*% solve(W) %*% fc)/n_b) |> as.vector()
        lambda_0 <- c(0,
                      exp(seq(from = log(1e-04*lambda_0_max),
                              to = log(lambda_0_max),
                              by = log(1e04)/(nlambda - 2)))
                      )
      }
      # if (lambda_2 == 0L){
      #   lambda_2 <- 1e-5
      # }

      # Shrinkage matrix
      if (G_bench == "Zero"){
        G_bench <- NULL
      } else{
        G_bench <- G
      }

      # Find optimal lambda_0 by minimizing sum of squared reconciled residuals
      if (is.null(train_data) | is.null(fitted_values)){
        stop("Training data and fitted values are required to find the optimal lambda_0")
      }
      
      if (parallel){
        future::plan(multisession, workers = workers)
        map_fun <- furrr::future_map
      } else {
        map_fun <- purrr::map
      }
      
      fit.out <- lambda_0 |>
        map_fun(\(l0) mip_l0(lambda_0 = l0, lambda_2 = lambda_2,
                             fc = fc, S = S, W = W, G_bench = G_bench,
                             M = M, solver = solver),
                .progress = .progress)
      
      index0 <- which(lambda_0 == 0)
      fit.G <- lapply(fit.out, function(l) l$G)
      fit.G[[index0]] <- G
      fit.z <- lapply(fit.out, function(l) l$z)
      sse <- purrr::map_dbl(fit.G, \(x) sum(stats::na.omit(train_data - fitted_values %*% t(x) %*% t(S))^2)) |> round(2)
      sse_summary <- data.frame(lambda0 = lambda_0, sse = sse)
      lambda_0 <- lambda_0[which.min(sse)]
      z <- fit.z[[which.min(sse)]]
      G <- fit.G[[which.min(sse)]]
    }
    
    # Subset selection
    # if (subset){
    #   # One-step ahead base forecasts
    #   fc <- base_forecasts[1, ] |> as.vector()
    #   
    #   # Candidate lambda_0
    #   if (is.null(lambda_0)){
    #     lambda_0_max <- (0.5 * (t(fc) %*% solve(W) %*% fc)/n_b) |> as.vector()
    #     lambda_0 <- c(0, 
    #                   exp(seq(from = log(1e-04*lambda_0_max), 
    #                           to = log(lambda_0_max), 
    #                           by = log(1e04)/(nlambda - 2)))
    #                   )
    #   }
    #   
    #   # Call python function
    #   py_fun <- function(y, S, W, l0, M){
    #     reticulate::use_python(pythonpath, required = T)
    #     reticulate::source_python(pyfunpath)
    #     A_diag <- miqp(y = y, S = S, W = W, l0 = l0, M = M)
    #     C <- t(S) %*% solve(W) %*% diag(A_diag) %*% S
    #     ev <- eigen(C, only.values = TRUE)[["values"]]
    #     if (any(ev < 1e-8)) {
    #       G <- NULL
    #     } else{
    #       G <- solve(C) %*% t(S) %*% solve(W) %*% diag(A_diag)
    #     }
    #     return(G)
    #   }
    #   
    #   # Find optimal lambda_0 by minimizing sum of squared reconciled residuals
    #   if (is.null(train_data) | is.null(fitted_values)){
    #     stop("Training data and fitted values are required to find the optimal lambda_0")
    #   }
    #   
    #   if (parallel){
    #     future::plan(multisession, workers = workers)
    #     map_fun <- furrr::future_map
    #   } else {
    #     map_fun <- purrr::map
    #   }
    #   
    #   G_candidates <- lambda_0 |>
    #     map_fun(\(l0) py_fun(l0 = l0, y = fc, S = S, W = W, M = 2),
    #             .progress = .progress)
    #   index <- !sapply(G_candidates, is.null)
    #   G_candidates <- G_candidates[index]
    #   lambda_0 <- lambda_0[index]
    #   sse <- purrr::map_dbl(G_candidates, \(x) sum(stats::na.omit(train_data - fitted_values %*% t(x) %*% t(S))^2)) |> round(2)
    #   sse_summary <- data.frame(lambda0 = lambda_0, sse = sse)
    #   lambda_0 <- lambda_0[which.min(sse)]
    #   G <- G_candidates[which.min(sse)]
    #   z <- ifelse(colSums(G) == 0, 0, 1)
    # }
    
    # Group lasso
    if (lasso){
      if (ELasso){
        # Group lasso based on in-sample observations and fitted values
        if (is.null(fitted_values) | is.null(train_data)){
          stop("Historical residuals are required to implement empirical group lasso")
        }
        if (NROW(fitted_values) != NROW(train_data)){
          stop("The dimensions of the fitted_values do not match train_data")
        }
        
        N <- NROW(train_data)
        X <- kronecker(S, fitted_values)
        y <- as.vector(train_data)
        index <- matrix(seq.int(n_b*n), nrow = n_b, ncol = n, byrow = FALSE) |> 
          t() |> 
          as.vector()
        
        # Reorder data to make 'group' as a vector of consecutive integers
        X <- X[ , order(index)]
        group <- rep(seq.int(n), each = n_b)
        
        cvfit <- gglasso::cv.gglasso(x = X, y = y, group = group, intercept = FALSE,
                                     loss = "ls", pred.loss = "L1", 
                                     nfolds = 5, nlambda = nlambda, lambda.factor = 1e-05)
        vec_G <- coef(cvfit, s = "lambda.min")[-1] # remove 0 intercept
        G <- matrix(vec_G, nrow = n_b, ncol = n, byrow = FALSE)
        sse_summary <- data.frame(lambda1 = cvfit$lambda, sse = cvfit$cvm)
        z <- ifelse(round(colSums(G), 3) == 0, 0, 1)
      } else{
        # Group lasso based on one-step ahead forecasts using different W (WLS)
        y_hat <- base_forecasts[1, ] |> as.vector() # One-step ahead base forecasts
        X <- kronecker(t(y_hat), S)
        group <- rep(seq.int(n), each = n_b)
        
        fit <- gglasso::gglasso(x = X, y = y_hat, group = group, intercept = FALSE,
                                loss = "wls", weight = solve(W), 
                                nlambda = nlambda, lambda.factor = 1e-05, eps = 1e-04)
        sse <- purrr::map_dbl(1:nlambda, \(i){
          G <- fit$beta[, i] |> as.vector() |> 
            matrix(nrow = n_b, ncol = n, byrow = FALSE)
          sum(stats::na.omit(train_data - fitted_values %*% t(G) %*% t(S))^2)
        }) |> round(2)
        lambda_index <- which.min(sse)
        G <- fit$beta[, lambda_index] |> as.vector() |> 
          matrix(nrow = n_b, ncol = n, byrow = FALSE)
        sse_summary <- data.frame(lambda1 = fit$lambda, sse = sse)
        z <- ifelse(round(colSums(G), 3) == 0, 0, 1)
      }
    } 
  }
  
  # Reconciliation
  y_tilde <- base_forecasts %*% t(G) %*% t(S)
  colnames(y_tilde) <- colnames(base_forecasts)
  
  # Output
  if (subset | lasso){
    z <- z
  } else{
    z <- NA
  }
  if (exists("sse_summary")){
    lambda_report <- sse_summary
  } else{
    lambda_report <- NA
  }
  list(y_tilde = y_tilde, G = G,
       z = z,
       lambda_report = lambda_report)
}
