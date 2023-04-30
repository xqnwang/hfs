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

#################################################
# Section 3.3-3.4 of Wickramasuriya et al. (2019)
#################################################
# Setup parameters
repeat_n = 500
n = 180
h = 16
freq = 4
sigma2_varrho = diag(rep(2, 4))
sigma2_varsigma = diag(rep(0.007, 4))
sigma2_varomega = diag(rep(7, 4))
eta_mat = diag(c(5, 4, 5, 4))
eta_mat[1, 2] = eta_mat[2, 1] = 3
eta_mat[1, 3] = eta_mat[3, 1] = 2
eta_mat[1, 4] = eta_mat[4, 1] = 1
eta_mat[2, 3] = eta_mat[3, 2] = 2
eta_mat[2, 4] = eta_mat[4, 2] = 1
eta_mat[3, 4] = eta_mat[4, 3] = 3
error1 = 10
error2 = 9

# Data simulation - STS
set.seed(123)
simulated_set = data.frame()
for (i in 1:repeat_n) {
  ts = STSsim(n, sigma2_varrho, sigma2_varsigma, sigma2_varomega, eta_mat, freq)
  ts = data.frame(ts, t=seq(as.Date("1978-01-01"), by="quarter", length.out = n),
                  index=i, row.names = NULL)
  simulated_set = rbind(simulated_set, ts)
}

colnames(simulated_set) <- c("AA", "AB", "BA", "BB", "Quarter", "Index")

write.csv(simulated_set, file = 'data/simulation_data.csv', row.names = FALSE)

# Data simulation - NOISY
set.seed(123)
simulated_set = data.frame()
for (i in 1:repeat_n) {
  ts = NOISYsim(n, sigma2_varrho, sigma2_varsigma, sigma2_varomega, eta_mat, error1, error2, freq)
  ts = data.frame(ts, t=seq(as.Date("1978-01-01"), by="quarter", length.out = n),
                  index=i, row.names = NULL)
  simulated_set = rbind(simulated_set, ts)
}

colnames(simulated_set) <- c("AA", "AB", "BA", "BB", "Quarter", "Index")

write.csv(simulated_set, file = 'data/simulation_data_noisy.csv', row.names = FALSE)

#################################################
# Negative error correlation
#################################################
# Setup parameters
repeat_n = 500
n = 180
h = 16
freq = 4
sigma2_varrho = diag(rep(2, 4))
sigma2_varsigma = diag(rep(0.007, 4))
sigma2_varomega = diag(rep(7, 4))
eta_mat = diag(rep(3, 4))
eta_mat[1, 2] = eta_mat[2, 1] = -2
eta_mat[3, 4] = eta_mat[4, 3] = -1

# Data simulation
set.seed(123)
simulated_set = data.frame()
for (i in 1:repeat_n) {
  ts = STSsim(n, sigma2_varrho, sigma2_varsigma, sigma2_varomega, eta_mat, freq)
  ts = data.frame(ts, t=seq(as.Date("1978-01-01"), by="quarter", length.out = n),
                  index=i, row.names = NULL)
  simulated_set = rbind(simulated_set, ts)
}

colnames(simulated_set) <- c("AA", "AB", "BA", "BB", "Quarter", "Index")

write.csv(simulated_set, file = 'data/simulation_data_negative.csv', row.names = FALSE)

