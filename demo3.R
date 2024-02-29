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
ri3 <- 2
n <- 32
phi.max <- 0.0
p1 <- 100
p2 <- 120
p3 <- 80
sigma1 <- 2
sigma2 <- 3
sigma3 <- 2
signal_strength1 <- 30
signal_strength2 <- 40
signal_strength3 <- 30
rank.spec <- 'exact'

# compute error
compute_fd <- function(P1, P2){
  sum(diag((diag(dim(P2)[1]) - P2) %*% P1))
}
compute_tp <- function(P1, P2){
  sum(diag(P1 %*% P2))
}

# models
compute <- list()
form_output <- function(joint, indiv1, indiv2, indiv3){
  check_null <- function(X){
    if (is.null(X)){
      return(list("r" = Inf,
                  "P" = matrix(0, nrow=n, ncol=n)))
    }
    return (list("r" = ncol(X),
                 "P" = X %*% t(X)))
  }
  check_null2 <- function(X1, X2) {
    if (is.null(X1) && is.null(X2)){
      return (list("r" = Inf,
                   "P" = matrix(0, nrow=n, ncol=n)))
    } else if (is.null(X1)) {
      return (list("r" = ncol(X2),
                   "P" = X2 %*% t(X2)))
    } else if (is.null(X2)) {
      return (list("r" = ncol(X1),
                   "P" = X1 %*% t(X1)))
    } else{
      return (list("r" = ncol(X1) + ncol(X2),
                   "P" = X1 %*% t(X1) + X2 %*% t(X2)))
    }
  }
  Pjoint <- check_null(joint)
  Pindiv1 <- check_null(indiv1)
  Pindiv2 <- check_null(indiv2)
  Pindiv3 <- check_null(indiv3)
  P1 <- check_null2(joint, indiv1)
  P2 <- check_null2(joint, indiv2)
  P3 <- check_null2(joint, indiv3)
  return (list("P1" = P1$P,
               "P2" = P2$P,
               "P3" = P3$P,
               "Pjoint" = Pjoint$P,
               "Pindiv1" = Pindiv1$P,
               "Pindiv2" = Pindiv2$P,
               "Pindiv3" = Pindiv3$P,
               "r1" = P1$r,
               "r2" = P2$r,
               "r3" = P3$r,
               "rj" = Pjoint$r,
               "ri1" = Pindiv1$r,
               "ri2" = Pindiv2$r,
               "ri3" = Pindiv3$r))
}
compute[["ajive"]] <- function(X1, X2, X3, rank1, rank2, rank3){
  out <- ajive(list(X1, X2, X3), c(rank1, rank2, rank3),
               n_wedin_samples = 100, 
               n_rand_dir_samples = 100)
  check_null <- function(X) {
    if (all(X == 0)){
      return (NULL)
    }
    return (X)
  }
  joint <- check_null(out$joint_scores)
  indiv1 <- check_null(out$block_decomps[[1]][['individual']][['u']])
  indiv2 <- check_null(out$block_decomps[[2]][['individual']][['u']])
  indiv3 <- check_null(out$block_decomps[[3]][['individual']][['u']])
  return (form_output(joint, indiv1, indiv2, indiv3))
}
compute[["proposed"]] <- function(X1, X2, X3, rank1, rank2, rank3) {
  U1.hat <- svd(X1)$u[, 1:rank1, drop = FALSE]
  U2.hat <- svd(X2)$u[, 1:rank2, drop = FALSE]
  U3.hat <- svd(X3)$u[, 1:rank3, drop = FALSE]
  P.hat <- U1.hat %*% t(U1.hat)
  Q.hat <- U2.hat %*% t(U2.hat)
  R.hat <- U3.hat %*% t(U3.hat)
  
  prod <- (P.hat %*% Q.hat %*% R.hat + 
             Q.hat %*% R.hat %*% P.hat + 
             R.hat %*% P.hat %*% Q.hat) / 3
  svd.prod <- svd(prod)
  cluster <- Ckmedian.1d.dp(sqrt(svd.prod$d), k=3)
  joint <- svd.prod$u[, cluster$cluster == 3, drop = FALSE]
  jointPerp <- diag(nrow(joint)) - joint %*% t(joint)
  
  P.hat <- jointPerp %*% P.hat
  svd.P.hat <- svd(P.hat)
  cluster <- Ckmedian.1d.dp(sqrt(svd.P.hat$d), k=2)
  indiv1 <- svd.P.hat$u[, cluster$cluster == 2, drop = FALSE]
  
  Q.hat <- jointPerp %*% Q.hat
  svd.Q.hat <- svd(Q.hat)
  cluster <- Ckmedian.1d.dp(sqrt(svd.Q.hat$d), k=2)
  indiv2 <- svd.Q.hat$u[, cluster$cluster == 2, drop = FALSE]

  R.hat <- jointPerp %*% R.hat
  svd.R.hat <- svd(R.hat)
  cluster <- Ckmedian.1d.dp(sqrt(svd.R.hat$d), k=2)
  indiv3 <- svd.R.hat$u[, cluster$cluster == 2, drop = FALSE]
  
  return (form_output(joint, indiv1, indiv2, indiv3))
}
compute[["proposed_subsampling"]] <- function(X1, X2, X3, rank1, rank2, rank3, numSamples=100) {
  avg.P <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  avg.P1 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  avg.P2 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  avg.P3 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  for (i in 1:numSamples){
    X1.sample <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
    X2.sample <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
    X3.sample <- X3[, sample(1:ncol(X3), as.integer(ncol(X3)/2), replace=FALSE)]
    svd.X1.sample <- svd(X1.sample)
    svd.X2.sample <- svd(X2.sample)
    svd.X3.sample <- svd(X3.sample)
    u1.sample <- svd.X1.sample$u[, 1:rank1]
    u2.sample <- svd.X2.sample$u[, 1:rank2]
    u3.sample <- svd.X3.sample$u[, 1:rank3]
    sample.P1 <- (u1.sample %*% t(u1.sample))
    sample.P2 <- (u2.sample %*% t(u2.sample))
    sample.P3 <- (u3.sample %*% t(u3.sample))
    avg.P1 <- avg.P1 + sample.P1
    avg.P2 <- avg.P2 + sample.P2
    avg.P3 <- avg.P3 + sample.P3
    prod <- (sample.P1 %*% sample.P2 %*% sample.P3 +
               sample.P2 %*% sample.P3 %*% sample.P1 +
               sample.P3 %*% sample.P1 %*% sample.P2) / 3
    avg.P <- avg.P + prod
  }
  avg.P1 <- avg.P1 / numSamples
  avg.P2 <- avg.P2 / numSamples
  avg.P3 <- avg.P3 / numSamples
  avg.P <- avg.P / numSamples
  
  svd.avg <- svd(avg.P)
  cluster <- Ckmeans.1d.dp(svd.avg$d, k=3)
  joint <- svd.avg$u[, cluster$cluster == 3, drop = FALSE]
  jointPerp <- diag(nrow(joint)) - joint %*% t(joint)
  
  avg.P1 <- jointPerp %*% avg.P1
  svd.avg1 <- svd(avg.P1)
  cluster <- Ckmeans.1d.dp(svd.avg1$d, k=2)
  indiv1 <- svd.avg1$u[, cluster$cluster == 2, drop = FALSE]
  
  avg.P2 <- jointPerp %*% avg.P2
  svd.avg2 <- svd(avg.P2)
  cluster <- Ckmeans.1d.dp(svd.avg2$d, k=2)
  indiv2 <- svd.avg2$u[, cluster$cluster == 2, drop = FALSE]

  avg.P3 <- jointPerp %*% avg.P3
  svd.avg3 <- svd(avg.P3)
  cluster <- Ckmeans.1d.dp(svd.avg3$d, k=2)
  indiv3 <- svd.avg3$u[, cluster$cluster == 2, drop = FALSE]
  
  return (form_output(joint, indiv1, indiv2, indiv3))
}

sim_iter <- 100
progress <- txtProgressBar(min=0, max=sim_iter, style=3, width = 100)
results <- list()
for (model in c("ajive", "proposed", "proposed_subsampling")) {
  results[[model]] <- matrix(0, sim_iter, 14)
}
for (i in 1:sim_iter) {
  dj1 <- rnorm(rj, mean = signal_strength1, sd = 3)
  dj2 <- rnorm(rj, mean = signal_strength2, sd = 3)
  dj3 <- rnorm(rj, mean = signal_strength3, sd = 3)
  di1 <- rnorm(ri1, mean = signal_strength1, sd = 3)
  di2 <- rnorm(ri2, mean = signal_strength2, sd = 3)
  di3 <- rnorm(ri3, mean = signal_strength3, sd = 3)
  # generate data
  U <- svd(matrix(rnorm(n * (p1+p2+p3)), n, p1+p2+p3))$u
  # joint part
  Uj <- U[, 1:rj]
  Ujperp <- U[, (rj+1):ncol(U)]
  O1 <- randortho(rj, type = 'orthonormal') # rotate Uj
  O2 <- randortho(rj, type = 'orthonormal') # rotate Uj
  Uj1 <- Uj 
  Uj2 <- Uj %*% O1
  Uj3 <- Uj %*% O2
  # individual part
  Ui1 <- U[, (rj+1):(rj+ri1)]
  Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
  Ui3 <- U[, (rj+ri1+ri2+1):(rj+ri1+ri2+ri3)]
  O1 <- matrix(runif(ri1 * ri2, -phi.max, phi.max), ri1, ri2) # rotate
  Ui2 <- Ui2 + Ui1 %*% O1
  Ui2 <- gramSchmidt(Ui2)$Q # orthonormalize
  O2 <- matrix(runif(ri2 * ri3, -phi.max, phi.max), ri2, ri3) # rotate
  Ui3 <- Ui3 + Ui2 %*% O2
  Ui3 <- gramSchmidt(Ui3)$Q # orthonormalize
  # loadings
  V1 <- svd(matrix(rnorm((rj + ri1) * p1), rj + ri1, p1))$v
  V2 <- svd(matrix(rnorm((rj + ri2) * p2), rj + ri2, p2))$v
  V3 <- svd(matrix(rnorm((rj + ri3) * p3), rj + ri3, p3))$v
  # combine
  Y1 <- cbind(Uj1, Ui1) %*% diag(c(dj1, di1)) %*% t(V1)
  Y2 <- cbind(Uj2, Ui2) %*% diag(c(dj2, di2)) %*% t(V2)
  Y3 <- cbind(Uj3, Ui3) %*% diag(c(dj3, di3)) %*% t(V3)
  # noise
  X1 <- Y1 + matrix(rnorm(n * p1), n, p1) * sigma1
  X2 <- Y2 + matrix(rnorm(n * p2), n, p2) * sigma2
  X3 <- Y3 + matrix(rnorm(n * p3), n, p3) * sigma3
  
  # compute true
  Pjoint <- Uj1 %*% t(Uj1) # true projection
  Pindiv1 <- Ui1 %*% t(Ui1)
  Pindiv2 <- Ui2 %*% t(Ui2)
  Pindiv3 <- Ui3 %*% t(Ui3)
  P1 <- Pjoint + Pindiv1
  P2 <- Pjoint + Pindiv2
  P3 <- Pjoint + Pindiv3
  
  # signal rank
  error1 <- 0
  # sample(1:2, 1)
  error2 <- 0
  # sample(1:2, 1)
  error3 <- 0
  # sample(1:2, 1)
  rank1 <- rj + ri1 + error1
  rank2 <- rj + ri2 + error2
  rank3 <- rj + ri3 + error3
  
  # compute results
  for (model in c("ajive", "proposed", "proposed_subsampling")) {
    out <- compute[[model]](X1, X2, X3, rank1, rank2, rank3)
    res <- c(compute_fd(out$P1, P1) / out$r1, compute_tp(out$P1, P1) / rank1,
             compute_fd(out$P2, P2) / out$r2, compute_tp(out$P2, P2) / rank2,
             compute_fd(out$P3, P3) / out$r3, compute_tp(out$P3, P3) / rank3,
             compute_fd(out$Pjoint, Pjoint) / out$rj, compute_tp(out$Pjoint, Pjoint) / rj,
             compute_fd(out$Pindiv1, Pindiv1) / out$ri1, compute_tp(out$Pindiv1, Pindiv1) / ri1,
             compute_fd(out$Pindiv2, Pindiv2) / out$ri2, compute_tp(out$Pindiv2, Pindiv2) / ri2,
             compute_fd(out$Pindiv3, Pindiv3) / out$ri3, compute_tp(out$Pindiv3, Pindiv3) / ri3)
    # check if any NA
    if (any(is.na(res))){
      print(paste0("NA detected in ", model))
      break
    }
    results[[model]][i,] <- res
  }
  setTxtProgressBar(progress, i)
}
close(progress)
results <- lapply(results, function(x) colMeans(x))
results <- do.call(rbind, results)
results <- as.data.frame(results)
colnames(results) <- c("fdr.P1", "tpr.P1",
                       "fdr.P2", "tpr.P2",
                       "fdr.P3", "tpr.P3",
                       "fdr.Pjoint", "tpr.Pjoint", 
                       "fdr.Pindiv1", "tpr.Pindiv1", 
                       "fdr.Pindiv2", "tpr.Pindiv2",
                       "fdr.Pindiv3", "tpr.Pindiv3")
# compute precision = 1 - fdr
results <- results %>% mutate(precision.Pjoint = 1 - fdr.Pjoint, 
                              precision.Pindiv1 = 1 - fdr.Pindiv1, 
                              precision.Pindiv2 = 1 - fdr.Pindiv2,
                              precision.Pindiv3 = 1 - fdr.Pindiv3)
# avg precision and avg tpr
results$avg_precision <- (results$precision.Pjoint + 
                            results$precision.Pindiv1 + 
                            results$precision.Pindiv2 + 
                            results$precision.Pindiv3) / 4
results$avg_tpr <- (results$tpr.Pjoint + 
                      results$tpr.Pindiv1 + 
                      results$tpr.Pindiv2 +
                      results$tpr.Pindiv3) / 4
# compute F1 score
results$f1 <- 2 * (results$avg_precision * results$avg_tpr) / (results$avg_precision + results$avg_tpr)
# save results
results.save <- list()
results.save[["results"]] <- results
results.save[["sim_iter"]] <- sim_iter
results.save[["rj"]] <- rj
results.save[["ri1"]] <- ri1
results.save[["ri2"]] <- ri2
results.save[["ri3"]] <- ri3
results.save[["n"]] <- n
results.save[["phi.max"]] <- phi.max
results.save[["p1"]] <- p1
results.save[["p2"]] <- p2
results.save[["p3"]] <- p3
results.save[["sigma1"]] <- sigma1
results.save[["sigma2"]] <- sigma2
results.save[["sigma3"]] <- sigma3
results.save[["signal_strength1"]] <- signal_strength1
results.save[["signal_strength2"]] <- signal_strength2
results.save[["signal_strength3"]] <- signal_strength3
results.save[["rank_spec"]] <- rank.spec
save(results.save, file=paste0("demo3_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData"))

