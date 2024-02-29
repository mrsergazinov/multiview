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
form_output <- function(joint, indiv1, indiv2){
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
  P1 <- check_null2(joint, indiv1)
  P2 <- check_null2(joint, indiv2)
  return (list("P1" = P1$P,
              "P2" = P2$P,
              "Pjoint" = Pjoint$P,
              "Pindiv1" = Pindiv1$P,
              "Pindiv2" = Pindiv2$P,
              "r1" = P1$r,
              "r2" = P2$r,
              "rj" = Pjoint$r,
              "ri1" = Pindiv1$r,
              "ri2" = Pindiv2$r))
}
compute[["ajive"]] <- function(X1, X2, rank1, rank2){
  out <- ajive(list(X1, X2), c(rank1, rank2),
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
  return (form_output(joint, indiv1, indiv2))
}
compute[["jive"]] <- function(X1, X2, rank1, rank2) {
  out <- jive(list(t(X1), t(X2)), rankA = c(rank1, rank2),
              method='perm', showProgress=FALSE)
  check_null <- function(X, rank){
    if (rank == 0){ return (NULL) }
    return (X[, 1:rank, drop = FALSE])
  }
  joint <- check_null(svd(t(out$joint[[1]]))$u, out$rankJ)
  indiv1 <- check_null(svd(t(out$individual[[1]]))$u, out$rankA[1])
  indiv2 <- check_null(svd(t(out$individual[[2]]))$u, out$rankA[2])
  return (form_output(joint, indiv1, indiv2))
}
compute[["slide"]] <- function(X1, X2, rank1, rank2) {
  out <- slide(cbind(X1, X2), pvec = c(ncol(X1), ncol(X2)))
  check_null <- function(X, mask){
    if (any(mask)){return(X[, mask, drop = FALSE])}
    return (NULL)
  }
  joint <- check_null(out$model$U, (out$S[1, ] == 1 & out$S[2, ] == 1))
  indiv1 <- check_null(out$model$U, (out$S[1, ] == 1 & out$S[2, ] == 0))
  indiv2 <- check_null(out$model$U, (out$S[1, ] == 0 & out$S[2, ] == 1))
  return (form_output(joint, indiv1, indiv2))
}
compute[["proposed"]] <- function(X1, X2, rank1, rank2) {
  U1.hat <- svd(X1)$u[, 1:rank1, drop = FALSE]
  U2.hat <- svd(X2)$u[, 1:rank2, drop = FALSE]
  P.hat <- U1.hat %*% t(U1.hat)
  Q.hat <- U2.hat %*% t(U2.hat)

  prod <- (P.hat %*% Q.hat + Q.hat %*% P.hat) / 2
  svd.prod <- svd(prod)
  cluster <- Ckmedian.1d.dp(sqrt(svd.prod$d), k = 3)
  joint <- svd.prod$u[, cluster$cluster == 3, drop = FALSE]
  jointPerp <- diag(nrow(joint)) - joint %*% t(joint)

  P.hat <- jointPerp %*% P.hat
  svd.P.hat <- svd(P.hat)
  cluster <- Ckmedian.1d.dp(sqrt(svd.P.hat$d), k = 2)
  indiv1 <- svd.P.hat$u[, cluster$cluster == 2, drop = FALSE]

  Q.hat <- jointPerp %*% Q.hat
  svd.Q.hat <- svd(Q.hat)
  cluster <- Ckmedian.1d.dp(sqrt(svd.Q.hat$d), k = 2)
  indiv2 <- svd.Q.hat$u[, cluster$cluster == 2, drop = FALSE]
  
  return (form_output(joint, indiv1, indiv2))
}
compute[["proposed_subsampling"]] <- function(X1, X2, rank1, rank2, numSamples=100) {
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
  cluster <- Ckmeans.1d.dp(svd.avg$d, k = 3)
  joint <- svd.avg$u[, cluster$cluster == 3, drop = FALSE]
  jointPerp <- diag(nrow(joint)) - joint %*% t(joint)

  avg.P1 <- jointPerp %*% avg.P1
  svd.avg1 <- svd(avg.P1)
  cluster <- Ckmeans.1d.dp(svd.avg1$d, k = 2)
  indiv1 <- svd.avg1$u[, cluster$cluster == 2, drop = FALSE]
  
  avg.P2 <- jointPerp %*% avg.P2
  svd.avg2 <- svd(avg.P2)
  cluster <- Ckmeans.1d.dp(svd.avg2$d, k = 2)
  indiv2 <- svd.avg2$u[,  cluster$cluster == 2, drop = FALSE]

  return (form_output(joint, indiv1, indiv2))
}

sim_iter <- 50
progress <- txtProgressBar(min=0, max=sim_iter, style=3, width = 100)
results <- list()
for (model in c("ajive", "jive", "slide", "proposed", "proposed_subsampling")){
  results[[model]] <- matrix(0, sim_iter, 10)
}
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
  for (model in c("ajive", "jive", "slide", "proposed", "proposed_subsampling")) {
    out <- compute[[model]](X1, X2, rank1, rank2)
    res <- c(compute_fd(out$P1, P1) / out$r1, compute_tp(out$P1, P1) / rank1,
             compute_fd(out$P2, P2) / out$r2, compute_tp(out$P2, P2) / rank2,
             compute_fd(out$Pjoint, Pjoint) / out$rj, compute_tp(out$Pjoint, Pjoint) / rj,
             compute_fd(out$Pindiv1, Pindiv1) / out$ri1, compute_tp(out$Pindiv1, Pindiv1) / ri1,
             compute_fd(out$Pindiv2, Pindiv2) / out$ri2, compute_tp(out$Pindiv2, Pindiv2) / ri2)
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
                      "fdr.Pjoint", "tpr.Pjoint", 
                      "fdr.Pindiv1", "tpr.Pindiv1", 
                      "fdr.Pindiv2", "tpr.Pindiv2")
# compute precision = 1 - fdr
results <- results %>% mutate(precision.Pjoint = 1 - fdr.Pjoint, 
                              precision.Pindiv1 = 1 - fdr.Pindiv1, 
                              precision.Pindiv2 = 1 - fdr.Pindiv2)
# avg precision and avg tpr
results$avg_precision <- (results$precision.Pjoint + 
                            results$precision.Pindiv1 + 
                            results$precision.Pindiv2) / 3
results$avg_tpr <- (results$tpr.Pjoint + 
                      results$tpr.Pindiv1 + 
                      results$tpr.Pindiv2) / 3
# compute F1 score
results$f1 <- 2 * (results$avg_precision * results$avg_tpr) / (results$avg_precision + results$avg_tpr)
# save results
results.save <- list()
results.save[["results"]] <- results
results.save[["sim_iter"]] <- sim_iter
results.save[["rj"]] <- rj
results.save[["ri1"]] <- ri1
results.save[["ri2"]] <- ri2
results.save[["n"]] <- n
results.save[["phi.max"]] <- phi.max
results.save[["p1"]] <- p1
results.save[["p2"]] <- p2
results.save[["sigma1"]] <- sigma1
results.save[["sigma2"]] <- sigma2
results.save[["signal_strength1"]] <- signal_strength1
results.save[["signal_strength2"]] <- signal_strength2
results.save[["rank_spec"]] <- rank.spec
save(results.save, file=paste0("demo2_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData"))
