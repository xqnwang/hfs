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
reconcile <- function(base_forecasts, S, 
                      method = c("bu", "ols", "wls_struct", "wls_var", "mint_cov", "mint_shrink"), 
                      residuals = NULL, fitted_values = NULL, train_data = NULL,
                      subset = FALSE, ridge = FALSE,
                      lambda_0 = NULL, lambda_2 = NULL, 
                      m = NULL, M = NULL, MIPGap = NULL, WarmStart = 1, MIPFocus = 0, Cuts = -1,
                      TimeLimit = 600, MIPVerbose = FALSE, SearchVerbose = FALSE){
  # Dimension info
  n <- NROW(S); nb <- NCOL(S)
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
  if (method == "bu"){
    # Buttom-up
    subset <- FALSE
    G <- cbind(matrix(0, nb, n-nb), diag(nb))
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
    
    # Group best-subset selection with ridge regularization
    if (subset){
      # One-step ahead base forecasts
      fc <- base_forecasts[1, ] |> matrix(nrow = n, ncol = 1)
      fc_tilde_init <- S %*% G %*% fc
      obj_init <- 0.5 * t(fc - fc_tilde_init) %*% solve(W) %*% (fc - fc_tilde_init) |> 
        as.numeric()
      
      # Find optimal lambda_0 by minimizing sum of squared reconciled residuals
      if (is.null(train_data) | is.null(fitted_values)){
        stop("Training data and fitted values are required to find the optimal lambda_0")
      }
      
      # Candidate lambda_0
      if (is.null(lambda_0)){
        ndigits <- floor(log10(abs(obj_init/nb))) + 2
        lambda_0 <- c(0, 10^seq(from = ndigits - 4, to = ndigits, by = 1))
      }
      
      # Candidate lambda_2
      if (ridge){
        if (is.null(lambda_2)){
          lambda_2 <- c(0, 10^seq(from = -2, to = 2, by = 1))
        }
      } else{
        lambda_2 <- 0
      }
      
      # Search optimal l0 and l2
      if (MIPVerbose){
        LogToConsole = 1
        OutputFlag = 1
      } else{
        LogToConsole = 0
        OutputFlag = 0
      }
      
      lambda <- expand.grid(l0 = lambda_0, l2 = lambda_2)
      mip.out <- purrr::pmap(lambda, function(l0, l2){
        if (l0 == 0 & l2 == 0){
          fit <- list(l0 = l0, l2 = l2, G = G, Z = rep(1, n), obj = obj_init, gap = 0, opt = 1)
        } else{
          fit <- miqp(fc, S, W, l0, l2, 
                      m, M, MIPGap, TimeLimit, 
                      LogToConsole, OutputFlag, 
                      WarmStart, MIPFocus, Cuts)
          names(fit) <- c("G", "Z", "obj", "gap", "opt")
          fit <- append(list(l0 = l0, l2 = l2), fit)
        }
        sse <- sum(stats::na.omit(train_data - fitted_values %*% t(fit$G) %*% t(S))^2)
        fit$sse <- sse
        fit
      }, .progress = SearchVerbose)

      sse_summary <- sapply(mip.out, function(l) c(l$l0, l$l2, sum(as.vector(l$Z)), l$sse, l$obj, l$gap, l$opt)) |> 
        t() |> 
        data.frame()
      names(sse_summary) <- c("lambda0", "lambda2", "k", "sse", "obj", "gap", "opt")
      min.index <- purrr::map_dbl(mip.out, function(l) l$sse) |> round(2) |> which.min()
      z <- mip.out[[min.index]]$Z |> as.vector()
      G <- mip.out[[min.index]]$G
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
    lambda_report <- sse_summary
  } else{
    lambda_report <- NA
  }
  list(y_tilde = y_tilde, G = G,
       z = z,
       lambda_report = lambda_report)
}
