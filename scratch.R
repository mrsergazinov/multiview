library(tidyverse)
library(ajive)
library(RMTstat)
library(pracma)

source("gen_data.R")
source("fd_control_joint.R")
source("ajive_oracle.R")
source("ajive.R")
source("bounds.R")

set.seed(235017)
rj <- 2
ri1 <- 3
ri2 <- 2
n <- 20
phi.max <- 0.07
p1 <- 100
p2 <- 100
sigma1 <- 1
sigma2 <- 1
sim_iter <- 100
signal_strength <- 40
dj <- rnorm(rj, mean = signal_strength, sd = 3)
di1 <- rnorm(ri1, mean = signal_strength, sd = 3)
di2 <- rnorm(ri2, mean = signal_strength, sd = 3)

# generate data
U <- svd(matrix(rnorm(n * (p1+p2)), n, p1+p2))$u
# joint part
Uj <- U[, 1:rj]
Ujperp <- U[, (rj+1):ncol(U)]
O <- randortho(rj, type = 'orthonormal') # rotate Uj
Uj1 <- Uj %*% O
Uj2 <- Uj
# individual part
Ui1 <- U[, (rj+1):(rj+ri1)]
Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
O <- matrix(runif(ri1 * ri2, -phi.max, phi.max), ri1, ri2) # rotate
Ui2 <- Ui2 + Ui1 %*% O
Ui2 <- gramSchmidt(Ui2)$Q # orthonormalize
# loadings
V1 <- svd(matrix(rnorm((rj + ri1) * p1), rj + ri1, p1))$v
V2 <- svd(matrix(rnorm((rj + ri2) * p2), rj + ri2, p2))$v
# combine
X1 <- cbind(Uj1, Ui1) %*% diag(c(dj, di1)) %*% t(V1)
X2 <- cbind(Uj2, Ui2) %*% diag(c(dj, di2)) %*% t(V2)
# noise
X1 <- X1 + matrix(rnorm(n * p1), n, p1) * sigma1
X2 <- X2 + matrix(rnorm(n * p2), n, p2) * sigma2

# compute bound
bound.val <- bound.approx.angle(c(dj, di1), c(dj, di2), rj, ri1, ri2, n / p1, 0, phi.max)
print(paste("Spectral bound = ", bound.val))

# FD control -- subsampling
args = list("sigma1" = NA, "sigma2" = NA,
            "rj" = rj, "ri1" = ri1,
            "ri2" = ri2, "numSamples" = 100,
            "alpha" = 0.4, "boundJoint" = bound.val)
out <- fd_control_joint(X1, X2, args)
estimJointProj.subsample <- out$joint %*% t(out$joint)

# FD control -- no subsampling
X1 <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
X2 <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
svd.X1 <- svd(X1)
svd.X2 <- svd(X2)
thresh.X1 <- thresh(X1) # thresholding singular values
thresh.X2 <- thresh(X2)
u1 <- svd.X1$u[, svd.X1$d > thresh.X1]
u2 <- svd.X2$u[, svd.X2$d > thresh.X2]
sample.P1 <- (u1 %*% t(u1))
sample.P2 <- (u2 %*% t(u2))
prod <- sample.P1 %*% sample.P2
svd.prod <- svd(prod)
joint <- svd.prod$u[, svd.prod$d > bound.val, drop = FALSE]
estimJointProj.noSubsample <- joint %*% t(joint) # no subsampling

tr.nosubsample <- c()
tr.subsample <- c()
for (i in 1:ncol(Ujperp)){
    u <- Ujperp[, i, drop = FALSE]
    tr.nosubsample <- c(tr.nosubsample, sum(diag(estimJointProj.noSubsample %*% u %*% t(u))))
    tr.subsample <- c(tr.subsample, sum(diag(estimJointProj.subsample %*% u %*% t(u))))
}

