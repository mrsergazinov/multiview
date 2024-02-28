library(tidyverse)
library(ajive)
library(r.jive)
library(SLIDE)
library(RMTstat)
library(pracma)
library(Ckmeans.1d.dp)

set.seed(235017)
rj <- 2
ri1 <- 3
ri2 <- 2
n <- 20
phi.max <- 0.5
p1 <- 100
p2 <- 120
sigma1 = 2
sigma2 = 3
signal_strength1 <- 30
signal_strength2 <- 40

# compute error
compute_fd <- function(P1, P2){
  sum(diag((diag(dim(P2)[1]) - P2) %*% P1))
}
compute_tp <- function(P1, P2){
  sum(diag(P1 %*% P2))
}

# models
check_na <- function(X, n){
  if (all(is.na(X))){
    X <- matrix(0, nrow=n, ncol=n)
  }
  return (X)
}
form_output <- function(joint, indiv1, indiv2){
  Pjoint <- check_na(joint %*% t(joint), nrow(X1))
  Pindiv1 <- check_na(indiv1 %*% t(indiv1), nrow(X1))
  Pindiv2 <- check_na(indiv2 %*% t(indiv2), nrow(X1))
  return (list("P1" = Pjoint + Pindiv1,
              "P2" = Pjoint + Pindiv2,
              "Pjoint" = Pjoint, 
              "Pindiv1" = Pindiv1, 
              "Pindiv2" = Pindiv2))
}
compute_ajive <- function(X1, X2, rank1, rank2){
  out <- ajive(list(X1, X2), c(rank1, rank2),
              n_wedin_samples = 100, 
              n_rand_dir_samples = 100)
  return (form_output(out$joint_scores, 
                     out$block_decomps[[1]][['individual']][['u']], 
                     out$block_decomps[[2]][['individual']][['u']]))
}
compute_jive <- function(X1, X2, rank1, rank2) {
  out <- jive(list(t(X1), t(X2)), rankA = c(rank1, rank2), 
              method='perm', showProgress=FALSE)
  joint <- svd(t(out$joint[[1]]))$u[, 1:out$rankJ]
  indiv1 <- svd(t(out$individual[[1]]))$u[, 1:out$rankA[1]]
  indiv2 <- svd(t(out$individual[[2]]))$u[, 1:out$rankA[2]]
  return (form_output(joint, indiv1, indiv2))
}
compute_slide <- function(X1, X2, rank1, rank2) {
  out <- slide(cbind(X1, X2), pvec = c(ncol(X1), ncol(X2)))
  joint <- out$model$U[, (out$S[1, ] == 1 & out$S[2, ] == 1)]
  indiv1 <- out$model$U[, (out$S[1, ] == 1 & out$S[2, ] == 0)]
  indiv2 <- out$model$U[, (out$S[1, ] == 0 & out$S[2, ] == 1)]
  return (form_output(joint, indiv1, indiv2))
}
compute_proposed <- function(X1, X2, rank1, rank2) {
  U1.hat <- svd(X1)$u[, 1:rank1, drop = FALSE]
  U2.hat <- svd(X2)$u[, 1:rank2, drop = FALSE]
  P.hat <- U1.hat %*% t(U1.hat)
  Q.hat <- U2.hat %*% t(U2.hat)

  prod <- (P.hat + Q.hat)
  svd.prod <- svd(prod)
  cluster <- Ckmedian.1d.dp(sqrt(svd.prod$d), k=3)
  joint <- svd.prod$u[, cluster$cluster == 3, drop = FALSE]
  jointPerp <- diag(nrow(joint)) - joint %*% t(joint)

  P.hat <- jointPerp %*% P.hat
  svd.P.hat <- svd(P.hat)
  indiv1 <- svd.P.hat$u[, svd.P.hat$d > 0.5, drop = FALSE]

  Q.hat <- jointPerp %*% Q.hat
  svd.Q.hat <- svd(Q.hat)
  indiv2 <- svd.Q.hat$u[, svd.Q.hat$d > 0.5, drop = FALSE]
  
  return (form_output(joint, indiv1, indiv2))
}
compute_proposed_subsampling <- function(X1, X2, rank1, rank2, numSamples) {
  avg.P <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  avg.P1 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  avg.P2 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  for (i in 1:numSamples){
    X1.sample <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
    X2.sample <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
    svd.X1.sample <- svd(X1.sample)
    svd.X2.sample <- svd(X2.sample)
    u1.sample <- svd.X1.sample$u[, 1:rank1]
    u2.sample <- svd.X2.sample$u[, 1:rank2]
    sample.P1 <- (u1.sample %*% t(u1.sample))
    sample.P2 <- (u2.sample %*% t(u2.sample))
    avg.P1 <- avg.P1 + sample.P1
    avg.P2 <- avg.P2 + sample.P2
    prod <- (sample.P1 %*% sample.P2 + sample.P2 %*% sample.P1) / 2
    avg.P <- avg.P + prod
  }
  avg.P1 <- avg.P1 / numSamples
  avg.P2 <- avg.P2 / numSamples
  avg.P <- avg.P / numSamples

  svd.avg <- svd(avg.P)
  cluster <- Ckmeans.1d.dp(svd.avg$d, k=3)
  joint <- svd.avg$u[, cluster$cluster == 3, drop = FALSE]
  jointPerp <- diag(nrow(joint)) - joint %*% t(joint)

  avg.P1 <- jointPerp %*% avg.P1
  svd.avg1 <- svd(avg.P1)
  indiv1 <- svd.avg1$u[, svd.avg1$d > 0.5, drop = FALSE]

  avg.P2 <- jointPerp %*% avg.P2
  svd.avg2 <- svd(avg.P2)
  indiv2 <- svd.avg2$u[, svd.avg2$d > 0.5, drop = FALSE]

  return (form_output(joint, indiv1, indiv2))
}

sim_iter <- 10
progress <- txtProgressBar(min=0, max=sim_iter, style=3, width = 100)
results <- list("ajive" = matrix(0, sim_iter, 10),
                "jive" = matrix(0, sim_iter, 10),
                "slide" = matrix(0, sim_iter, 10),
                "proposed" = matrix(0, sim_iter, 10),
                "proposed_subsampling" = matrix(0, sim_iter, 10))
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
  Pjoint <- Uj1 %*% t(Uj1) # true projection
  Pindiv1 <- Ui1 %*% t(Ui1)
  Pindiv2 <- Ui2 %*% t(Ui2)
  P1 <- Pjoint + Pindiv1
  P2 <- Pjoint + Pindiv2
  
  # signal rank
  error1 <- 0
    # sample(1:2, 1)
  error2 <- 0
    # sample(1:2, 1)
  rank1 <- rj + ri1 + error1
  rank2 <- rj + ri2 + error2

  # compute results
  out <- compute_ajive(X1, X2, rank1, rank2)
  results[["ajive"]][i,] <- c(compute_fd(out$P1, P1), compute_tp(out$P1, P1),
                              compute_fd(out$P2, P2), compute_tp(out$P2, P2),
                              compute_fd(out$Pjoint, Pjoint), compute_tp(out$Pjoint, Pjoint),
                              compute_fd(out$Pindiv1, Pindiv1), compute_tp(out$Pindiv1, Pindiv1),
                              compute_fd(out$Pindiv2, Pindiv2), compute_tp(out$Pindiv2, Pindiv2))
  out <- compute_jive(X1, X2, rank1, rank2)
  results[["jive"]][i,] <- c(compute_fd(out$P1, P1), compute_tp(out$P1, P1),
                              compute_fd(out$P2, P2), compute_tp(out$P2, P2),
                              compute_fd(out$Pjoint, Pjoint), compute_tp(out$Pjoint, Pjoint),
                              compute_fd(out$Pindiv1, Pindiv1), compute_tp(out$Pindiv1, Pindiv1),
                              compute_fd(out$Pindiv2, Pindiv2), compute_tp(out$Pindiv2, Pindiv2))
  out <- compute_slide(X1, X2, rank1, rank2)
  results[["slide"]][i,] <- c(compute_fd(out$P1, P1), compute_tp(out$P1, P1),
                              compute_fd(out$P2, P2), compute_tp(out$P2, P2),
                              compute_fd(out$Pjoint, Pjoint), compute_tp(out$Pjoint, Pjoint),
                              compute_fd(out$Pindiv1, Pindiv1), compute_tp(out$Pindiv1, Pindiv1),
                              compute_fd(out$Pindiv2, Pindiv2), compute_tp(out$Pindiv2, Pindiv2))
  out <- compute_proposed(X1, X2, rank1, rank2)
  results[["proposed"]][i,] <- c(compute_fd(out$P1, P1), compute_tp(out$P1, P1),
                              compute_fd(out$P2, P2), compute_tp(out$P2, P2),
                              compute_fd(out$Pjoint, Pjoint), compute_tp(out$Pjoint, Pjoint),
                              compute_fd(out$Pindiv1, Pindiv1), compute_tp(out$Pindiv1, Pindiv1),
                              compute_fd(out$Pindiv2, Pindiv2), compute_tp(out$Pindiv2, Pindiv2))
  out <- compute_proposed_subsampling(X1, X2, rank1, rank2, 100)
  results[["proposed_subsampling"]][i,] <- c(compute_fd(out$P1, P1), compute_tp(out$P1, P1),
                              compute_fd(out$P2, P2), compute_tp(out$P2, P2),
                              compute_fd(out$Pjoint, Pjoint), compute_tp(out$Pjoint, Pjoint),
                              compute_fd(out$Pindiv1, Pindiv1), compute_tp(out$Pindiv1, Pindiv1),
                              compute_fd(out$Pindiv2, Pindiv2), compute_tp(out$Pindiv2, Pindiv2))
  setTxtProgressBar(progress, i)
}
close(progress)
results <- lapply(results, function(x) colMeans(x))
results <- do.call(rbind, results)
results <- as.data.frame(results)
colnames(results) <- c("fd.P1", "tp.P1", 
                      "fd.P2", "tp.P2", 
                      "fd.Pjoint", "tp.Pjoint", 
                      "fd.Pindiv1", "tp.Pindiv1", 
                      "fd.Pindiv2", "tp.Pindiv2")
