# Utility functions
cumsum_matrix <- function(x) apply(x, 2, cumsum)
brodcast_sum <- function(x, y){
  t(apply(x, 1, function(row) row + y))
}

# Simulate seasonal time series using the basic structural time series (STS) model
STSsim <- function(n, sigma2_varrho, sigma2_varsigma, sigma2_omega, eta_mat, m){
  nseries <- dim(sigma2_omega)[1]
  
  # level
  initial_level <- rnorm(nseries)
  varrho <- MASS::mvrnorm(n-1, rep(0, nseries), sigma2_varrho)
  level <- rbind(initial_level,
                 brodcast_sum(cumsum_matrix(varrho), initial_level))
  
  # trend
  initial_trend <- rnorm(nseries)
  varsigma <- MASS::mvrnorm(n-1, rep(0, nseries), sigma2_varsigma)
  trend <- rbind(initial_trend, 
                 brodcast_sum(cumsum_matrix(varsigma), initial_trend))
  
  # seasonal
  initial_seasonal <- MASS::mvrnorm(m-1, rep(0, nseries), diag(rep(1, nseries)))
  omega <- MASS::mvrnorm(n-m+1, rep(0, nseries), sigma2_omega)
  season <- initial_seasonal
  for (i in m:n) {
    season <- rbind(season, -apply(season[(i-m+1):(i-1),], 2, sum) + omega[i-m+1,])
  }
  
  # eta
  eta <- MTS::VARMAsim(n, sample(c(0, 1), 1, replace = TRUE), sample(c(0, 1), 1, replace = TRUE),
                       phi = diag(runif(4, 0.5, 0.7)), theta = diag(runif(4, 0.5, 0.7)),
                       sigma = eta_mat)$series
  
  bts <- level + trend + season + eta
  bts
}

# Exploring the Smoothing Effect of Aggregation - noisy time series in disaggregated levels
NOISYsim <- function(n, sigma2_varrho, sigma2_varsigma, sigma2_omega, eta_mat, error1, error2, m=1){
  nseries <- dim(sigma2_omega)[1]
  
  # level
  initial_level <- rnorm(nseries)
  varrho <- MASS::mvrnorm(n-1, rep(0, nseries), sigma2_varrho)
  level <- rbind(initial_level,
                 brodcast_sum(cumsum_matrix(varrho), initial_level))
  
  # trend
  initial_trend <- rnorm(nseries)
  varsigma <- MASS::mvrnorm(n-1, rep(0, nseries), sigma2_varsigma)
  trend <- rbind(initial_trend, 
                 brodcast_sum(cumsum_matrix(varsigma), initial_trend))
  
  # seasonal
  initial_seasonal <- MASS::mvrnorm(m-1, rep(0, nseries), diag(rep(1, nseries)))
  omega <- MASS::mvrnorm(n-m+1, rep(0, nseries), sigma2_omega)
  season <- initial_seasonal
  for (i in m:n) {
    season <- rbind(season, -apply(season[(i-m+1):(i-1),], 2, sum) + omega[i-m+1,])
  }
  
  # eta
  eta <- MTS::VARMAsim(n, sample(c(0, 1), 1, replace = TRUE), sample(c(0, 1), 1, replace = TRUE),
                       phi = diag(runif(4, 0.5, 0.7)), theta = diag(runif(4, 0.5, 0.7)),
                       sigma = eta_mat)$series
  
  bts <- level + trend + season + eta
  
  e1 <- rnorm(n, 0, sqrt(error1))
  e2 <- rnorm(n, 0, sqrt(error2))
  bts[, 1] = bts[, 1] - e1 - 0.5*e2
  bts[, 2] = bts[, 2] + e1 - 0.5*e2
  bts[, 3] = bts[, 3] - e1 + 0.5*e2
  bts[, 4] = bts[, 4] + e1 + 0.5*e2
  bts
}