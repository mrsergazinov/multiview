
# d <- c()
# r1 <- 3
# r2 <- 3
# U1 <- svd(matrix(rnorm(100 * 50), nrow = 100, ncol = 50))$u[, 1:r1]
# V1 <-  svd(matrix(rnorm(100 * 50), nrow = 100, ncol = 50))$v[, 1:r1]
# X1 <- U1 %*% diag(rnorm(r1, mean=10)) %*% t(V1)
# V2 <-  svd(matrix(rnorm(100 * 60), nrow = 100, ncol = 60))$v[, 1:r2]
# X2 <- U1 %*% diag(rnorm(r2, mean=10)) %*% t(V2)
# 
# sing.vals <- c()
# for (i in 1:100) {
#   U1.hat <- svd(X1 + matrix(rnorm(100 * 50), 100, 50))$u[, 1:r1]
#   U2.hat <- svd(X2 + matrix(rnorm(100 * 60), 100, 60))$u[, 1:r2]
#   P <- U1.hat %*% t(U1.hat)
#   Q <- U2.hat %*% t(U2.hat)
#   sing.vals <- c(svd(P %*% Q)$d, sing.vals)
# }

library(RMTstat)

source("gen_data.R")
source("fd_control_joint.R")

# set simulation parameters
set.seed(235017)
rj <- 2
ri1 <- 2
ri2 <- 2
n <- 100
p1 <- 20
p2 <- 20
sigma1 <- 0.5
sigma2 <- 0.5
sim_iter <- 100
signal_strength <- 20
dj <- rnorm(rj, mean = signal_strength, sd = 3)
di1 <- rnorm(ri1, mean = signal_strength, sd = 3)
di2 <- rnorm(ri2, mean = signal_strength, sd = 3)

sing.vals <- c()
for (i in 1:sim_iter){
  data <- gen_data(n, p1, p2, rj, ri1, ri2, dj, di1, di2, sigma1, sigma2)
  X1 <- data[["X1"]]
  X2 <- data[["X2"]]
  Uj <- data[["joint"]]
  Ui1 <- data[["indiv1"]]
  Ui2 <- data[["indiv2"]]
  jointProj <- Uj %*% t(Uj) # true projection
  indiv1Proj <- Ui1 %*% t(Ui1)
  indiv2Proj <- Ui2 %*% t(Ui2)
  
  thresh.X1 <- thresh(X1, sigma=NA) # thresholding singular values
  thresh.X2 <- thresh(X2, sigma=NA)
  svd.X1 <- svd(X1)
  svd.X2 <- svd(X2)
  u1 <- svd.X1$u[, svd.X1$d > thresh.X1 + 1e-10] # thresholding using Gavish and Donoho 2014
  u2 <- svd.X2$u[, svd.X2$d > thresh.X2 + 1e-10]
  P1 <- (u1 %*% t(u1))
  P2 <- (u2 %*% t(u2))
  prod <- P1 %*% P2
  sing.vals <- svd(prod)$d
}



