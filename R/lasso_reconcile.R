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
#' @param nlambda Number of candidate `lambda_1` values.
#' @param nfolds Number of folds in cross-validation
#' 
#' @import future
#' @import gglasso
#' @import reticulate
#' 
#' @export
lasso.reconcile <- function(base_forecasts, S, 
                            method = c("bu", "ols", "wls_struct", "wls_var", "mint_cov", "mint_shrink"), 
                            residuals = NULL, fitted_values = NULL, train_data = NULL, nvalid = NULL,
                            lasso = NULL, nlambda = 20,
                            TimeLimit = 0, MonARCH = FALSE, workers = 4){
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
      if (NROW(train_data) != NROW(fitted_values)){
        stop("The dimensions of the fitted_values do not match train_data")
      }
      nvalid <- ifelse(is.null(nvalid), NROW(train_data), nvalid)
      
      # Candidate lambda_1
      X <- kronecker(t(fc), S)
      w <- 1/apply(G, 2, function(x) sqrt(sum(x^2)))
      lambda_max <- sapply(1:n, function(j){
        w_j <- w[j]
        cj <- ((j-1)*nb + 1):(j*nb)
        delta_loss_j <- - t(X[, cj]) %*% solve(W) %*% fc
        return(sqrt(sum(delta_loss_j^2))/w_j)
      }) |> max()
      lambda_min <- 0.0001 * lambda_max
      lambda <- c(sapply(0:(nlambda-1), function(j) lambda_max*(lambda_min/lambda_max)^(j/(nlambda-1))), 0)
      
      cl <- parallel::makeCluster(workers)
      doParallel::registerDoParallel(cl)
      socp.out <- foreach::foreach(l1 = lambda) %dopar% {
        if (MonARCH){
          path <- "~/.local/share/r-miniconda/envs/r-reticulate/bin/python3.8"
          setwd(Sys.glob(file.path("~/wm15/", "*", "hfs")))
        } else{
          path <- "~/Library/r-miniconda-arm64/bin/python3.10"
        }
        reticulate::use_python(path, required = T)
        reticulate::source_python("Python/glasso.py")
        
        if (l1 == 0){
          fit <- list(l1 = l1, G = G, Z = rep(1, n), obj = obj_init)
        } else{
          fit <- glasso(y = fc, S = S, W = W, l1 = l1, weight = 1, unbiased = 1, TimeLimit = TimeLimit, LogToConsole = 0, OutputFlag = 0)
          names(fit) <- c("G", "Z", "obj")
          fit <- append(list(l1 = l1), fit)
        }
        sse <- sum(stats::na.omit(tail(train_data, nvalid) - tail(fitted_values, nvalid) %*% t(fit$G) %*% t(S))^2)
        fit$sse <- sse
        return(fit)
      }
      stopCluster(cl)
      
      # socp.out <- purrr::map(lambda, function(l1){
      #   if (l1 == 0){
      #     fit <- list(l1 = l1, G = G, Z = rep(1, n), obj = obj_init)
      #   } else{
      #     fit <- glasso(y = fc, S = S, W = W, l1 = l1, weight = 1, unbiased = 1, TimeLimit = TimeLimit, LogToConsole = 0, OutputFlag = 0)
      #     names(fit) <- c("G", "Z", "obj")
      #     fit <- append(list(l1 = l1), fit)
      #   }
      #   sse <- sum(stats::na.omit(train_data - fitted_values %*% t(fit$G) %*% t(S))^2)
      #   fit$sse <- sse
      #   fit
      # })
      
      sse_summary <- sapply(socp.out, function(l) c(l$l1, sum(as.vector(l$Z)), l$sse, l$obj)) |> 
        t() |> 
        data.frame()
      names(sse_summary) <- c("lambda1", "k", "sse", "obj")
      min.index <- purrr::map_dbl(socp.out, function(l) l$sse) |> round(2) |> which.min()
      lambda_opt <- socp.out[[min.index]]$l1
      names(lambda_opt) <- "l1"
      z <- socp.out[[min.index]]$Z |> as.vector()
      G <- socp.out[[min.index]]$G
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
      if (N < 20) {
        stop("At least 20 observation should be used in ELasso")
      }
      nvalid <- ifelse(is.null(nvalid), floor(N/10), nvalid)
      
      # Candidate lambda_1
      w <- 1/apply(G, 2, function(x) sqrt(sum(x^2)))
      lambda_max <- sapply(1:n, function(j){
        w_j <- w[j]
        cj <- (0:(nb-1))*n + j
        delta_loss_j <- - 1/N * t(kronecker(S, fitted_values)[, cj]) %*% as.vector(train_data)
        return(sqrt(sum(delta_loss_j^2))/w_j)
      }) |> max()
      lambda_min <- 0.0001 * lambda_max
      lambda <- c(sapply(0:(nlambda-1), function(j) lambda_max*(lambda_min/lambda_max)^(j/(nlambda-1))), 0)

      # Find the optimal lambda_1 by splitting training data
      Y <- head(train_data, -nvalid)
      Y_hat <- head(fitted_values, -nvalid)
      cl <- parallel::makeCluster(workers)
      doParallel::registerDoParallel(cl)
      socp.out <- foreach::foreach(l1 = lambda) %dopar% {
        if (MonARCH){
          path <- "~/.local/share/r-miniconda/envs/r-reticulate/bin/python3.8"
          setwd(Sys.glob(file.path("~/wm15/", "*", "hfs")))
        } else{
          path <- "~/Library/r-miniconda-arm64/bin/python3.10"
        }
        reticulate::use_python(path, required = T)
        reticulate::source_python("Python/eglasso.py")
        
        fit <- eglasso(Y = Y, Y_hat = Y_hat, S = S, l1 = l1, weight = 1, 
                       TimeLimit = TimeLimit, LogToConsole = 0, OutputFlag = 0)
        names(fit) <- c("G", "Z", "obj")
        fit <- append(list(l1 = l1), fit)
        sse <- sum(stats::na.omit(tail(train_data, nvalid) - tail(fitted_values, nvalid) %*% t(fit$G) %*% t(S))^2)
        fit$sse <- sse
        return(fit)
      }
      stopCluster(cl)
      
      # socp.out <- purrr::map(lambda, function(l1){
      #   fit <- eglasso(Y = Y, Y_hat = Y_hat, S = S, l1 = l1, weight = 1, 
      #                  TimeLimit = TimeLimit, LogToConsole = 0, OutputFlag = 0)
      #   names(fit) <- c("G", "Z", "obj")
      #   fit <- append(list(l1 = l1), fit)
      #   sse <- sum(stats::na.omit(tail(train_data, floor(N/10)) - tail(fitted_values, floor(N/10)) %*% t(fit$G) %*% t(S))^2)
      #   fit$sse <- sse
      #   fit
      # })
      
      sse_summary <- sapply(socp.out, function(l) c(l$l1, sum(as.vector(l$Z)), l$sse, l$obj)) |> 
        t() |> 
        data.frame()
      names(sse_summary) <- c("lambda1", "k", "sse", "obj")
      min.index <- purrr::map_dbl(socp.out, function(l) l$sse) |> round(2) |> which.min()
      lambda_opt <- socp.out[[min.index]]$l1
      names(lambda_opt) <- "l1"
      
      # Resolve model to get G matrix
      if (MonARCH){
        path <- "~/.local/share/r-miniconda/envs/r-reticulate/bin/python3.8"
        setwd(Sys.glob(file.path("~/wm15/", "*", "hfs")))
      } else{
        path <- "~/Library/r-miniconda-arm64/bin/python3.10"
      }
      reticulate::use_python(path, required = T)
      reticulate::source_python("Python/eglasso.py")
      fit_opt <- eglasso(Y = train_data, Y_hat = fitted_values, S = S, l1 = lambda_opt, weight = 1, 
                         TimeLimit = TimeLimit, LogToConsole = 0, OutputFlag = 0)
      names(fit_opt) <- c("G", "Z", "obj")
      z <- fit_opt$Z |> as.vector()
      G <- fit_opt$G
    }
  }
  
  # Reconciliation
  y_tilde <- base_forecasts %*% t(G) %*% t(S)
  colnames(y_tilde) <- colnames(base_forecasts)
  
  # Output
  if (subset){
    z <- z
    lambda_opt <- lambda_opt
  } else{
    z <- NA
    lambda_opt <- NA
  }
  if (exists("sse_summary")){
    lambda_report <- sse_summary
  } else{
    lambda_report <- NA
  }
  list(y_tilde = y_tilde, G = G,
       z = z,
       lambda_report = lambda_report,
       lambda_opt = lambda_opt)
}
