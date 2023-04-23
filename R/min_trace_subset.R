min_trace_subset <- function(models,
                             method = c("wls_var", "ols", "wls_struct", "mint_cov", "mint_shrink"),
                             subset = TRUE, lasso = FALSE, ridge = FALSE,
                             lambda_0 = NULL, lambda_1 = 0, lambda_2 = 0,
                             nlambda_0 = 20,
                             sparse = NULL){
  if(is.null(sparse)){
    sparse <- requireNamespace("Matrix", quietly = TRUE)
  }
  structure(models, class = c("lst_mintsst_mdl", "lst_mdl", "list"),
            method = match.arg(method),
            subset = subset, lasso = lasso, ridge = ridge,
            lambda_0 = lambda_0, lambda_1 = lambda_1, lambda_2 = lambda_2,
            nlambda_0 = nlambda_0,
            sparse = sparse)
}

forecast.lst_mintsst_mdl <- function(object, key_data, 
                                     new_data = NULL, h = NULL,
                                     point_forecast = list(.mean = mean), ...){
  method <- object%@%"method"
  subset <- object%@%"subset"
  lasso <- object%@%"lasso"
  ridge <- object%@%"ridge"
  lambda_0 <- object%@%"lambda_0"
  lambda_1 <- object%@%"lambda_1"
  lambda_2 <- object%@%"lambda_2"
  nlambda_0 <- object%@%"nlambda_0"
  sparse <- object%@%"sparse"
  
  if(sparse){
    fabletools:::require_package("Matrix")
    as.matrix <- Matrix::as.matrix
    t <- Matrix::t
    diag <- function(x) if(is.vector(x)) Matrix::Diagonal(x = x) else Matrix::diag(x)
    solve <- Matrix::solve
    cov2cor <- Matrix::cov2cor
  } else {
    cov2cor <- stats::cov2cor
  }
  
  if(!subset){
    lambda_0 <- 0L
    message("Value of lambda_0 is set to 0 to match the subset argument.")
  }
  if(!lasso){
    if(lambda_1 != 0L){
      lambda_1 <- 0L
      message("Value of lambda_1 is set to 0 to match the lasso argument.")
    }
  }
  if(!ridge){
    if(lambda_2 != 0L){
      lambda_2 <- 0L
      message("Value of lambda_2 is set to 0 to match the ridge argument.")
    }
  }
  # if(lambda_1 == 0L & lambda_2 == 0L){
  #   # Give a very small lasso/ridge penalty in order to return the MinT explicit solution when lambda = c(0L, 0L, 0L)
  #   lambda_1 <- 1e-4
  # }
  
  point_method <- point_forecast
  point_forecast <- list()
  # Get forecasts
  fc <- NextMethod()
  if(length(unique(map(fc, interval))) > 1){
    abort("Reconciliation of temporal hierarchies is not yet supported.")
  }
  
  # Extract historical values
  values <- map(object, function(x, ...) x$data)
  if(length(unique(map_dbl(values, nrow))) > 1){
    values <- unname(as.matrix(reduce(values, full_join, by = index_var(values[[1]]))[,-1]))
  } else {
    values <- matrix(invoke(c, map(values, `[[`, 2)), ncol = length(object))
  }
  
  # Extract fitted values
  fits <- map(object, function(x, ...) fitted(x, ...), type = "response")
  if(length(unique(map_dbl(fits, nrow))) > 1){
    fits <- unname(as.matrix(reduce(fits, full_join, by = index_var(fits[[1]]))[,-1]))
  } else {
    fits <- matrix(invoke(c, map(fits, `[[`, 2)), ncol = length(object))
  }
  
  # Compute weights (sample covariance)
  res <- map(object, function(x, ...) residuals(x, ...), type = "response")
  if(length(unique(map_dbl(res, nrow))) > 1){
    # Join residuals by index #199
    res <- unname(as.matrix(reduce(res, full_join, by = index_var(res[[1]]))[,-1]))
  } else {
    res <- matrix(invoke(c, map(res, `[[`, 2)), ncol = length(object))
  }
  
  # Construct S matrix
  agg_data <- fabletools:::build_key_data_smat(key_data)
  
  n <- nrow(res)
  covm <- crossprod(stats::na.omit(res)) / n
  if(method == "ols"){
    # OLS
    W <- diag(rep(1L, nrow(covm)))
  } else if(method == "wls_var"){
    # WLS variance scaling
    W <- diag(diag(covm))
  } else if (method == "wls_struct"){
    # WLS structural scaling
    W <- diag(vapply(agg_data$agg,length,integer(1L)))
  } else if (method == "mint_cov"){
    # min_trace covariance
    W <- covm
  } else if (method == "mint_shrink"){
    # min_trace shrink
    tar <- diag(apply(res, 2, compose(crossprod, stats::na.omit))/n)
    corm <- cov2cor(covm)
    xs <- scale(res, center = FALSE, scale = sqrt(diag(covm)))
    xs <- xs[stats::complete.cases(xs),]
    v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
    diag(v) <- 0
    corapn <- cov2cor(tar)
    d <- (corm - corapn)^2
    lambda <- sum(v)/sum(d)
    lambda <- max(min(lambda, 1), 0)
    W <- lambda * tar + (1 - lambda) * covm
  } else {
    abort("Unknown reconciliation method")
  }
  
  # Check positive definiteness of weights
  eigenvalues <- eigen(W, only.values = TRUE)[["values"]]
  if (any(eigenvalues < 1e-8)) {
    abort("min_trace needs covariance matrix to be positive definite.", call. = FALSE)
  }
  
  # Reconciliation matrices
  if(sparse){
    row_btm <- agg_data$leaf
    row_agg <- seq_len(nrow(key_data))[-row_btm]
    S <- Matrix::sparseMatrix(
      i = rep(seq_along(agg_data$agg), lengths(agg_data$agg)),
      j = vctrs::vec_c(!!!agg_data$agg),
      x = rep(1, sum(lengths(agg_data$agg))))
    J <- Matrix::sparseMatrix(i = S[row_btm,,drop = FALSE]@i+1, j = row_btm, x = 1L,
                              dims = rev(dim(S)))
    U <- cbind(
      Matrix::Diagonal(diff(dim(J))),
      -S[row_agg,,drop = FALSE]
    )
    U <- U[, order(c(row_agg, row_btm)), drop = FALSE]
    Ut <- t(U)
    WUt <- W %*% Ut
    P0 <- J - J %*% WUt %*% solve(U %*% WUt, U)
  }
  else {
    S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vctrs::vec_c(!!!agg_data$agg)))
    S[length(agg_data$agg)*(vctrs::vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L
    R <- t(S)%*%solve(W)
    P0 <- solve(R%*%S)%*%R
  }
  
  # Extract one-step ahead base forecasts
  fc_h1 <- map_vec(fc, function(f) pull(f, -1) |> mean() |> head(1)) # use the resulted P matrix for all forecast horizons
  
  # Generate candidate lambda_0
  if(is.null(lambda_0)){
    lambda_0_max <- (t(fc_h1) %*% solve(W) %*% fc_h1)/NCOL(S)
    lambda_0 <- c(0, exp(seq(from = log(1e-04*lambda_0_max), to = log(lambda_0_max), 
                             by = log(1e04)/(nlambda_0 - 2))))
  }
  
  find_lambda_0 <- function(lambda_0, lambda_1, lambda_2,
                            fc_h1, S, W, G0,
                            M, solver,
                            values, fits){
    fit.mip_i <- mip_l0(fc = fc_h1, S = S, W = W, G0 = G0,
                        lambda_0 = lambda_0, lambda_1 = lambda_1, lambda_2 = lambda_2,
                        M = M,
                        solver = solver)
    SSE <- sum((as.vector(values) - kronecker(S, fits) %*% as.vector(t(fit.mip_i$G)))^2, na.rm = TRUE)
    return(SSE)
  }
  SSE <- lambda_0 |>
    map_dbl(find_lambda_0, lambda_1 = lambda_1, lambda_2 = lambda_2,
            fc_h1, S, W, G0 = P0,
            M = NULL, solver = "gurobi",
            values, fits)
  # SSE <- NULL
  # for(i in seq.int(length(lambda_0))){
  #   fit.mip_i <- mip_l0(fc_h1, S, W, G0 = P0,
  #                       lambda_0 = lambda_0[i], lambda_1 = lambda_1, lambda_2 = lambda_2,
  #                       M = NULL,
  #                       solver = "gurobi")
  #   SSE[i] <- sum((as.vector(values) - kronecker(S, fits) %*% as.vector(t(fit.mip_i$G)))^2, na.rm = TRUE)
  # }
  
  # Estimate P matrix using the selected optimal lambda_0
  lambda_0_select <- lambda_0[which.min(SSE)]
  fit.mip <- mip_l0(fc_h1, S, W, G0 = P0,
                    lambda_0 = lambda_0_select, lambda_1 = lambda_1, lambda_2 = lambda_2,
                    M = NULL,
                    solver = "gurobi")
  P <- fit.mip$G
  
  fabletools:::reconcile_fbl_list(fc, S, P, W, point_forecast = point_method)
}
