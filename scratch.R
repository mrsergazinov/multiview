
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
source("bounds.R")

# set simulation parameters
set.seed(235017)
rj <- 3
ri1 <- 2
ri2 <- 2
n <- 100
p1 <- 20
p2 <- 20
sigma1 <- 1
sigma2 <- 1
sim_iter <- 100
signal_strength <- 40
dj <- rnorm(rj, mean = signal_strength, sd = 3)
di1 <- rnorm(ri1, mean = signal_strength, sd = 3)
di2 <- rnorm(ri2, mean = signal_strength, sd = 3)
# compute spectral bound
bound.val <- bound(c(dj, di1), c(dj, di2), rj, ri1, ri2, p1 / n)
print(paste("Spectral bound = ", bound.val))

thresh <- function(X, sigma = NA) {
  # Based on ...
  if (is.na(sigma)){
    sigma = median(svd(X)$d) / sqrt(qmp(0.5, nrow(X), ncol(X)))
  }
  return (sigma * (1 + sqrt(min(ncol(X), nrow(X)) / max(ncol(X), nrow(X)))))
}

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
  
  # sampling
  X1.sample <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
  X2.sample <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
  # svd of sample + thresholding
  svd.X1.sample <- svd(X1.sample)
  svd.X2.sample <- svd(X2.sample)
  thresh.X1 <- thresh(X1.sample, sigma=1) # thresholding singular values
  thresh.X2 <- thresh(X2.sample, sigma=1)
  u1.sample <- svd.X1.sample$u[, svd.X1.sample$d > thresh.X1]
  u2.sample <- svd.X2.sample$u[, svd.X2.sample$d > thresh.X2]
  # print(paste("dim after thresh of U1: ", dim(u1.sample)[2]))
  # print(paste("dim after thresh of U2: ", dim(u1.sample)[2]))
  # compute projection
  sample.P1 <- (u1.sample %*% t(u1.sample))
  sample.P2 <- (u2.sample %*% t(u2.sample))
  # compute product
  prod <- sample.P1 %*% sample.P2
  # print singular values
  svd.prod <- svd(prod)
  print(paste("max sing val of PQ:", max(svd.prod$d)))
  # save sing.vals of prod
  sing.vals <- c(sing.vals, svd.prod$d)
}

for (i in seq(0, 0.9, 0.1)) {
  print(paste("Range [", i, i+0.1, "]: ", sum((sing.vals >= i) & (sing.vals < i+0.1))))
}



