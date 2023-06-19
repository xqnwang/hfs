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
#' @param lasso Data to use. Valid option are `NULL` no lasso, `Lasso` for out-of-sample based lasso and `ELasso` for in-sample based lasso.
#' @param nlambda Number of candidate `lambda_1` values to choose from for group lasso problem.
#' @param nfolds Number of folds in cross-validation
#' @param SearchVerbose Logical. If true, a progress bar will be displayed when searching optimal combination of `lambda_0` and `lambda_2`.
#' 
#' @import future
#' @import gglasso
#' @import reticulate
#' 
#' @export
lasso.reconcile <- function(base_forecasts, S, 
                            method = c("bu", "ols", "wls_struct", "wls_var", "mint_cov", "mint_shrink"), 
                            residuals = NULL, fitted_values = NULL, train_data = NULL,
                            lasso = NULL, nlambda = 20, nfolds = 5,
                            TimeLimit = 0, SearchVerbose = FALSE){
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
  
  if (is.null(lasso)){
    subset <- FALSE
    method <- match.arg(method)
  } else if(lasso == "Lasso"){
    subset <- TRUE
    method <- match.arg(method)
  } else if(lasso == "ELasso"){
    subset <- TRUE
    method <- "ols"
  }
  
  if (method == "bu"){
    # Buttom-up
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
  }
  
  if (subset){
    # Group lasso with the unbiasedness constraint using different W (WLS)
    if (lasso == "Lasso"){
      # One-step ahead base forecasts
      fc <- base_forecasts[1, ] |> matrix(nrow = n, ncol = 1)
      fc_tilde_init <- S %*% G %*% fc
      obj_init <- 0.5 * t(fc - fc_tilde_init) %*% solve(W) %*% (fc - fc_tilde_init) |> 
      as.numeric()
      
      # Find optimal lambda_1 by minimizing sum of squared reconciled residuals
      if (is.null(train_data) | is.null(fitted_values)){
        stop("Training data and fitted values are required to find the optimal lambda_0")
      }
      
      # Candidate lambda_1
      ndigits <- floor(log10(abs(obj_init))) + 2
      lambda_max <- 10^ndigits
      lambda <- c(0,
                    exp(seq(from = log(1e-05*lambda_max),
                            to = log(lambda_max),
                            by = log(1e05)/(nlambda - 2)))
      )
      
      socp.out <- purrr::map(lambda, function(l1){
        if (l1 == 0){
          fit <- list(l1 = l1, G = G, Z = rep(1, n), obj = obj_init)
        } else{
          fit <- socp(y = fc, S = S, W = W, l1 = l1, weight = 1, unbiased = 1, TimeLimit = TimeLimit, LogToConsole = 0, OutputFlag = 0)
          names(fit) <- c("G", "Z", "obj")
          fit <- append(list(l1 = l1), fit)
        }
        sse <- sum(stats::na.omit(train_data - fitted_values %*% t(fit$G) %*% t(S))^2)
        fit$sse <- sse
        fit
      }, .progress = SearchVerbose)
      
      sse_summary <- sapply(socp.out, function(l) c(l$l1, sum(as.vector(l$Z)), l$sse, l$obj)) |> 
        t() |> 
        data.frame()
      names(sse_summary) <- c("lambda1", "k", "sse", "obj")
      min.index <- purrr::map_dbl(socp.out, function(l) l$sse) |> round(2) |> which.min()
      z <- socp.out[[min.index]]$Z |> as.vector()
      G <- socp.out[[min.index]]$G
      
      # # Directly using group lasso without unbiasedness constraint (can specify pf - penalty weights applied to each group)
      # y_hat <- base_forecasts[1, ] |> as.vector() # One-step ahead base forecasts
      # X <- kronecker(t(y_hat), S)
      # group <- rep(seq.int(n), each = nb)
      # 
      # fit <- gglasso::gglasso(x = X, y = y_hat, group = group, intercept = FALSE,
      #                         loss = "wls", weight = solve(W), 
      #                         nlambda = nlambda, lambda.factor = 1e-05, eps = 1e-04)
      # sse <- purrr::map_dbl(1:nlambda, \(i){
      #   G <- fit$beta[, i] |> as.vector() |> 
      #     matrix(nrow = nb, ncol = n, byrow = FALSE)
      #   sum(stats::na.omit(train_data - fitted_values %*% t(G) %*% t(S))^2)
      # }) |> round(2)
      # lambda_index <- which.min(sse)
      # G <- fit$beta[, lambda_index] |> as.vector() |> 
      #   matrix(nrow = nb, ncol = n, byrow = FALSE)
      # sse_summary <- data.frame(lambda1 = fit$lambda, sse = sse)
      # z <- ifelse(round(colSums(abs(G)), 5) == 0, 0, 1)
    }
    
    # Empirical group lasso
    if (lasso == "ELasso"){
      if (is.null(fitted_values) | is.null(train_data)){
        stop("Historical residuals are required to implement empirical group lasso")
      }
      if (NROW(fitted_values) != NROW(train_data)){
        stop("The dimensions of the fitted_values do not match train_data")
      }
      
      N <- NROW(train_data)
      X <- kronecker(S, fitted_values)
      y <- as.vector(train_data)
      index <- matrix(seq.int(nb*n), nrow = nb, ncol = n, byrow = FALSE) |> 
        t() |> 
        as.vector()
      
      # Separate penalty weights can be applied to each group
      pf <- 1/apply(G, 2, function(x) sqrt(sum(x^2))) |> as.vector()
      
      # Reorder data to make 'group' as a vector of consecutive integers
      X <- X[ , order(index)]
      group <- rep(seq.int(n), each = nb)
      
      # Cross-validation for gglasso
      if (N < 20) {
        stop("At least 20 observation should be used in ELasso")
      }
      eachfoldid <- sample(rep(seq(nfolds), length = N))
      foldid <- rep(eachfoldid, n)
      cvfit <- gglasso::cv.gglasso(x = X, y = y, group = group, intercept = FALSE, 
                                   pf = pf,
                                   loss = "ls", pred.loss = "L1", 
                                   foldid = foldid, nlambda = nlambda, lambda.factor = 1e-05,
                                   eps = 1e-04)
      vec_G <- coef(cvfit, s = "lambda.min")[-1] # remove 0 intercept
      G <- matrix(vec_G, nrow = nb, ncol = n, byrow = FALSE)
      sse_summary <- data.frame(lambda1 = cvfit$lambda, 
                                sse = cvfit$cvm)
      z <- ifelse(colSums(abs(G) > 1e-8) == 0, 0, 1)
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
