library(tidyverse)
library(ajive)
library(RMTstat)
library(pracma)
library(Ckmeans.1d.dp)

set.seed(235017)
rj <- 2
ri1 <- 3
ri2 <- 2
n <- 20
phi.max <- 0.8
p1 <- 100
p2 <- 120
sigma1 = 2
sigma2 = 3
signal_strength1 <- 30
signal_strength2 <- 40

sim_iter <- 100
results_fd <- list(joint_tp = rep(NA, sim_iter), joint_fd = rep(NA, sim_iter), 
                   indiv1_tp = rep(NA, sim_iter), indiv1_fd = rep(NA, sim_iter),
                   indiv2_tp = rep(NA, sim_iter), indiv2_fd = rep(NA, sim_iter))
results_ajive <- list(joint_tp = rep(NA, sim_iter), joint_fd = rep(NA, sim_iter), 
                      indiv1_tp = rep(NA, sim_iter), indiv1_fd = rep(NA, sim_iter),
                      indiv2_tp = rep(NA, sim_iter), indiv2_fd = rep(NA, sim_iter))
# results_fd_sub <- list(joint_tp = rep(NA, sim_iter), joint_fd = rep(NA, sim_iter), 
#                        indiv1_tp = rep(NA, sim_iter), indiv1_fd = rep(NA, sim_iter),
#                        indiv2_tp = rep(NA, sim_iter), indiv2_fd = rep(NA, sim_iter))
for (i in 1:sim_iter) {
  dj1 <- rnorm(rj, mean = signal_strength1, sd = 3)
  dj2 <- rnorm(rj, mean = signal_strength2, sd = 3)
  di1 <- rnorm(ri1, mean = signal_strength1, sd = 3)
  di2 <- rnorm(ri2, mean = signal_strength2, sd = 3)
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
  Y1 <- cbind(Uj1, Ui1) %*% diag(c(dj1, di1)) %*% t(V1)
  Y2 <- cbind(Uj2, Ui2) %*% diag(c(dj2, di2)) %*% t(V2)
  # noise
  X1 <- Y1 + matrix(rnorm(n * p1), n, p1) * sigma1
  X2 <- Y2 + matrix(rnorm(n * p2), n, p2) * sigma2

  # compute true
  jointProj <- Uj1 %*% t(Uj1) # true projection
  indiv1Proj <- Ui1 %*% t(Ui1)
  indiv2Proj <- Ui2 %*% t(Ui2)
  
  # signal rank
  # error1 <- sample(-1:2, 1)
  # error2 <- sample(-1:2, 1)
  rank1 <- rj + ri1
  rank2 <- rj + ri2

  # AJIVE
  out <- ajive(list(X1, X2), c(rank1, rank2),
              n_wedin_samples = 100, 
              n_rand_dir_samples = 100)
  joint <- out$joint_scores
  indiv1 <- out$block_decomps[[1]][['individual']][['u']]
  indiv2 <- out$block_decomps[[2]][['individual']][['u']]
  estimJointProj <- joint %*% t(joint) # compute projection matrix
  estimIndiv1Proj <- indiv1 %*% t(indiv1)
  estimIndiv2Proj <- indiv2 %*% t(indiv2)
  if (all(is.na(joint))){
    estimJointProj <- matrix(0, nrow=n, ncol=n)
  }
  if (all(is.na(indiv1))){
    estimIndiv1Proj <- matrix(0, nrow=n, ncol=n)
  }
  if (all(is.na(indiv2))){
    estimIndiv2Proj <- matrix(0, nrow=n, ncol=n)
  }
  results_ajive$joint_tp[i] <- sum(diag(jointProj %*% estimJointProj))
  results_ajive$joint_fd[i] <- sum(diag((diag(n)-jointProj) %*% estimJointProj))
  results_ajive$indiv1_tp[i] <- sum(diag(indiv1Proj %*% estimIndiv1Proj))
  results_ajive$indiv1_fd[i] <- sum(diag((diag(n)-indiv1Proj) %*% estimIndiv1Proj))
  results_ajive$indiv2_tp[i] <- sum(diag(indiv2Proj %*% estimIndiv2Proj))
  results_ajive$indiv2_fd[i] <- sum(diag((diag(n)-indiv2Proj) %*% estimIndiv2Proj))
  
  # Proposed method
  # compute estimated P, Q
  U1.hat <- svd(X1)$u[, 1:rank1, drop = FALSE]
  U2.hat <- svd(X2)$u[, 1:rank2, drop = FALSE]
  P.hat <- U1.hat %*% t(U1.hat)
  Q.hat <- U2.hat %*% t(U2.hat)

  prod <- (P.hat %*% Q.hat + Q.hat %*% P.hat) / 2
  svd.prod <- svd(prod)
  cluster <- Ckmedian.1d.dp(c(1, svd.prod$d), k=3)
  joint <- svd.prod$u[, cluster$cluster[-1] == 3, drop = FALSE]
  jointPerp <- diag(nrow(joint)) - joint %*% t(joint)

  P.hat <- jointPerp %*% P.hat
  svd.P.hat <- svd(P.hat)
  indiv1 <- svd.P.hat$u[, svd.P.hat$d > 0.5, drop = FALSE]

  Q.hat <- jointPerp %*% Q.hat
  svd.Q.hat <- svd(Q.hat)
  indiv2 <- svd.Q.hat$u[, svd.Q.hat$d > 0.5, drop = FALSE]

  estimJointProj <- joint %*% t(joint) # compute projection matrix
  estimIndiv1Proj <- indiv1 %*% t(indiv1)
  estimIndiv2Proj <- indiv2 %*% t(indiv2)
  if (all(is.na(joint))){
    estimJointProj <- matrix(0, nrow=n, ncol=n)
  }
  if (all(is.na(indiv1))){
    estimIndiv1Proj <- matrix(0, nrow=n, ncol=n)
  }
  if (all(is.na(indiv2))){
    estimIndiv2Proj <- matrix(0, nrow=n, ncol=n)
  }

  results_fd$joint_tp[i] <- sum(diag(jointProj %*% estimJointProj))
  results_fd$joint_fd[i] <- sum(diag((diag(n)-jointProj) %*% estimJointProj))
  results_fd$indiv1_tp[i] <- sum(diag(indiv1Proj %*% estimIndiv1Proj))
  results_fd$indiv1_fd[i] <- sum(diag((diag(n)-indiv1Proj) %*% estimIndiv1Proj))
  results_fd$indiv2_tp[i] <- sum(diag(indiv2Proj %*% estimIndiv2Proj))
  results_fd$indiv2_fd[i] <- sum(diag((diag(n)-indiv2Proj) %*% estimIndiv2Proj))
  
  # Proposed method with FD
  # thresh <- function(X, sigma = NA) {
  #   if (is.na(sigma)){
  #     sigma = median(svd(X)$d) / sqrt(qmp(0.5, ncol(X), nrow(X)))
  #   }
  #   beta = nrow(X) / ncol(X)
  #   lambda = 1 + sqrt(beta)
  #   return (sigma * lambda)
  # }
  # avg.P <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  # avg.P1 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  # avg.P2 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  # for (j in 1:100){
  #   X1.sample <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
  #   X2.sample <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
  #   svd.X1.sample <- svd(X1.sample)
  #   svd.X2.sample <- svd(X2.sample)
  #   thresh.X1 <- thresh(X1.sample) # thresholding singular values
  #   thresh.X2 <- thresh(X2.sample)
  #   u1.sample <- svd.X1.sample$u[, svd.X1.sample$d > thresh.X1]
  #   u2.sample <- svd.X2.sample$u[, svd.X2.sample$d > thresh.X2]
  #   sample.P1 <- (u1.sample %*% t(u1.sample))
  #   sample.P2 <- (u2.sample %*% t(u2.sample))
  #   avg.P1 <- avg.P1 + sample.P1
  #   avg.P2 <- avg.P2 + sample.P2
  #   prod <- sample.P1 %*% sample.P2
  #   avg.P <- avg.P + prod
  # }
  # avg.P1 <- avg.P1 / 100
  # avg.P2 <- avg.P2 / 100
  # avg.P <- avg.P / 100
  # 
  # svd.avg <- svd(avg.P)
  # cluster <- Ckmeans.1d.dp(svd.avg$d, k=3)
  # joint <- svd.avg$u[, 1:cluster$cluster == 3, drop = FALSE]
  # jointPerp <- diag(nrow(joint)) - joint %*% t(joint)
  # 
  # avg.P1 <- jointPerp %*% avg.P1
  # svd.avg1 <- svd(avg.P1)
  # indiv1 <- svd.avg1$u[, svd.avg1$d > 0.5, drop = FALSE]
  # 
  # avg.P2 <- jointPerp %*% avg.P2
  # svd.avg2 <- svd(avg.P2)
  # indiv2 <- svd.avg2$u[, svd.avg2$d > 0.5, drop = FALSE]
  # 
  # estimJointProj <- joint %*% t(joint) # compute projection matrix
  # estimIndiv1Proj <- indiv1 %*% t(indiv1)
  # estimIndiv2Proj <- indiv2 %*% t(indiv2)
  # if (all(is.na(joint))){
  #   estimJointProj <- matrix(0, nrow=n, ncol=n)
  # }
  # if (all(is.na(indiv1))){
  #   estimIndiv1Proj <- matrix(0, nrow=n, ncol=n)
  # }
  # if (all(is.na(indiv2))){
  #   estimIndiv2Proj <- matrix(0, nrow=n, ncol=n)
  # }
  # 
  # results_fd_sub$joint_tp[i] <- sum(diag(jointProj %*% estimJointProj))
  # results_fd_sub$joint_fd[i] <- sum(diag((diag(n)-jointProj) %*% estimJointProj))
  # results_fd_sub$indiv1_tp[i] <- sum(diag(indiv1Proj %*% estimIndiv1Proj))
  # results_fd_sub$indiv1_fd[i] <- sum(diag((diag(n)-indiv1Proj) %*% estimIndiv1Proj))
  # results_fd_sub$indiv2_tp[i] <- sum(diag(indiv2Proj %*% estimIndiv2Proj))
  # results_fd_sub$indiv2_fd[i] <- sum(diag((diag(n)-indiv2Proj) %*% estimIndiv2Proj))
}

# compute mean
results_ajive <- map(results_ajive, mean)
results_fd <- map(results_fd, mean)
# results_fd_sub <- map(results_fd_sub, mean)
# make table of AJIVE and FD results
results <- rbind(results_ajive, results_fd)
# , results_fd_sub)
rownames(results) <- c("AJIVE", "FD")
# , "FD_sub")
