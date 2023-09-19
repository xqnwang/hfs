# Shrinking the sample covariance matrix towards a digonal matrix

# Functions required to get Shrinked covariance matrix
# Target matrix for shrinking towards a diagonal matrix
mytargetD <- function (x) 
{
  diag(colMeans(x^2))
}

# Shrinked covariance matrix - Schafer and strimmer approach
# Shrinkage intensity parameter is based on the correlation matrix
# Then DRD is used to get back the covariance matrix

myshrink.estim <- function (x, tar) 
{
  if (is.matrix(x) == TRUE && is.numeric(x) == FALSE) 
    stop("The data matrix must be numeric!")
  p <- ncol(x)
  n <- nrow(x)
  covm <- crossprod(x) / n
  corm <- cov2cor(covm)
  xs <- scale(x, center = FALSE, scale = sqrt(diag(covm)))
  v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
  diag(v) <- 0
  # corapn <- cov2cor(tar)
  corapn <- diag(p)
  d <- (corm - corapn)^2
  lambda <- sum(v)/sum(d)
  lambda <- max(min(lambda, 1), 0)
  shrink.cov <- lambda * tar + (1 - lambda) * covm
  return(list(shrink.cov, c("The shrinkage intensity lambda is:", 
                            round(lambda, digits = 4))))
}
