mint_subset <- function(models,
                        method = c("wls_var", "ols", "wls_struct", "mint_cov", "mint_shrink"),
                        subset = FALSE,
                        lambda = 0L,
                        BigM = NULL,
                        sparse = NULL){
  if(is.null(sparse)){
    sparse <- requireNamespace("Matrix", quietly = TRUE)
  }
  structure(models, class = c("lst_mintadj_mdl", "lst_mdl", "list"),
            method = match.arg(method),
            subset = subset,
            lambda = lambda,
            BigM = BigM,
            sparse = sparse)
}

forecast.lst_mintadj_mdl <- function(object, key_data, 
                                     new_data = NULL, h = NULL,
                                     point_forecast = list(.mean = mean), ...){
  method <- object%@%"method"
  subset <- object%@%"subset"
  lambda_0 <- object%@%"lambda"
  M <- object%@%"BigM"
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
    if(lambda_0 != 0L){
      lambda_0 <- 0L
      message("Value of lambda_0 is set to 0 to match the subset argument.")
    }
  }
  
  point_method <- point_forecast
  point_forecast <- list()
  # Get forecasts
  fc <- NextMethod()
  if(length(unique(map(fc, interval))) > 1){
    abort("Reconciliation of temporal hierarchies is not yet supported.")
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
  
  fc_h1 <- map_vec(fc, function(f) pull(f, -1) |> mean() |> head(1)) # extract one-step ahead base forecasts and use the resulted P matrix for all forecast horizons
  fit.mip <- mint_adjust(fc_h1, S, W, G = P0,
                         lambda_0 = lambda_0, M,
                         solver = "gurobi")
  P <- fit.mip$G
  
  fabletools:::reconcile_fbl_list(fc, S, P, W, point_forecast = point_method)
}
