#################################################
# Section 3.3-3.4 of Wickramasuriya et al. (2019)
#################################################
source("R/STSsim.R")
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


#################################################
# Generating data for bottom level from VAR(1) processes
#################################################
source("R/cov-shrink.R")
library(MTS) # generating multivariate time series
library(forecast) # for fitting models and forecasting
library(hts) # For hierarchical forecasting
library(Matrix)

# Constructing A1 matrix
########################################
theta <- pi/3
r1 <- 0.6
z1 <- r1 * (cos(theta) + 1i * sin(theta))
z2 <- Conj(z1)

## Below A1 matrix is the same as putting real part of z1 along the diagonal and imaginary part along the off-diagonal with opposite signs i.e. [u v; -v w]
a1 <- d1 <- r1 * cos(theta)
b1 <- r1 * sin(theta)
c1 <- -r1 * sin(theta)

A1 <- matrix(c(a1, c1, b1, d1), 2, 2)

# Constructing A2 matrix
########################################
alpha <- pi/6
r2 <- 0.9
z3 <- r2 * (cos(alpha) + 1i * sin(alpha))
z4 <- Conj(z3)

## Below A2 matrix is the same as putting real part of z3 along the diagonal and imaginary part along the off-diagonal with opposite signs i.e. [u v; -v w]
a2 <- d2 <- r2 * cos(alpha)
b2 <- r2 * sin(alpha)
c2 <- -r2 * sin(alpha)

A2 <- matrix(c(a2, c2, b2, d2), 2, 2)

# Phi and Sigma
########################################
Phi <- as.matrix(bdiag(A1, A2))

sig2.AA <- sig2.BA <- 2
sig2.AB <- sig2.BB <- 3
rho <- seq(-0.8, 0.8, by = 0.2)
constant <- rep(1, 4)

# Simulating data
########################################
set.seed(123)
nsim <- 500
n <- 101
h <- 1
nfreq <- 1L

nodes <- list(2, c(2, 2))
gmat <- hts:::GmatrixH(nodes)
gmat <- apply(gmat, 1, table)
nt <- sum(unlist(nodes)) + 1
nb <- sum(nodes[[length(nodes)]])

for (k in 1:length(rho)){
  sigmaA <- matrix(c(sig2.AA, sqrt(sig2.AA * sig2.AB) * rho[k], sqrt(sig2.AA * sig2.AB) * rho[k], sig2.AB), 2, 2)
  sigma <- bdiag(sigmaA, sigmaA)
  
  j <- 1
  simulated_set = data.frame()
  while (j <= nsim) {
    allB <- VARMAsim(nobs = n, arlags = 1, phi = Phi, cnst = constant,
                     skip = 200, sigma = sigma)$series
    allB <- data.frame(allB, Time = 1:n, Index = j)
    simulated_set <- rbind(simulated_set, allB)
    j <- j + 1
  }
  colnames(simulated_set) <- c("AA", "AB", "BA", "BB", "Time", "Index")
  write.csv(simulated_set, file = paste0("data/corr_", k, "_data.csv"), row.names = FALSE)
}


