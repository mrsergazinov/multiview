library(tidyverse)
library(ajive)
library(r.jive)
library(SLIDE)
library(RMTstat)
library(pracma)
library(Ckmeans.1d.dp)

set.seed(235017)
rj <- 4
ri1 <- 3
ri2 <- 2
ri3 <- 4
m <- 50
phi_max <- 0.8
n1 <- 80
n2 <- 100
n3 <- 110
sigma1 <- 1
sigma2 <- 1
sigma3 <- 1
signal_strength1 <- 10
signal_strength2 <- 12
signal_strength3 <- 15
rank_spec <- 'exact'
no_joint <- FALSE
no_indiv <- FALSE
try(if (no_joint && no_indiv) stop("At least one of no_joint and no_indiv must be FALSE"))

# set args from command line
args <- commandArgs(trailingOnly = TRUE)
# Check if arguments are provided
if (length(args) > 0) {
  # Parse command-line arguments
  for (i in seq_along(args)) {
    print(args[[i]])
    eval(parse(text = args[[i]]))
  }
}

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
                  "P" = matrix(0, nrow=m, ncol=m)))
    }
    return (list("r" = ncol(X),
                 "P" = X %*% t(X)))
  }
  check_null2 <- function(Y1, Y2) {
    if (is.null(Y1) && is.null(Y2)){
      return (list("r" = Inf,
                   "P" = matrix(0, nrow=m, ncol=m)))
    } else if (is.null(Y1)) {
      return (list("r" = ncol(Y2),
                   "P" = Y2 %*% t(Y2)))
    } else if (is.null(Y2)) {
      return (list("r" = ncol(Y1),
                   "P" = Y1 %*% t(Y1)))
    } else{
      return (list("r" = ncol(Y1) + ncol(Y2),
                   "P" = Y1 %*% t(Y1) + Y2 %*% t(Y2)))
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
compute[["ajive"]] <- function(Y1, Y2, Y3, rank1, rank2, rank3){
  out <- ajive(list(Y1, Y2, Y3), c(rank1, rank2, rank3),
               n_wedin_samples = 100, 
               n_rand_dir_samples = 100)
  check_null <- function(X) {
    if (is.na(X)) {
      return (NULL)
    } else if (all(X == 0)){
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
compute[["proposed"]] <- function(Y1, Y2, Y3, rank1, rank2, rank3) {
  U1.hat <- svd(Y1)$u[, 1:rank1, drop = FALSE]
  U2.hat <- svd(Y2)$u[, 1:rank2, drop = FALSE]
  U3.hat <- svd(Y3)$u[, 1:rank3, drop = FALSE]
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
compute[["proposed_subsampling1"]] <- function(Y1, Y2, Y3, rank1, rank2, rank3, numSamples=100) {
  avg.P <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
  avg.P1 <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
  avg.P2 <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
  avg.P3 <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
  for (i in 1:numSamples){
    Y1.sample <- Y1[, sample(1:ncol(Y1), as.integer(ncol(Y1)/2), replace=FALSE)]
    Y2.sample <- Y2[, sample(1:ncol(Y2), as.integer(ncol(Y2)/2), replace=FALSE)]
    Y3.sample <- Y3[, sample(1:ncol(Y3), as.integer(ncol(Y3)/2), replace=FALSE)]
    svd.Y1.sample <- svd(Y1.sample)
    svd.Y2.sample <- svd(Y2.sample)
    svd.Y3.sample <- svd(Y3.sample)
    u1.sample <- svd.Y1.sample$u[, 1:rank1]
    u2.sample <- svd.Y2.sample$u[, 1:rank2]
    u3.sample <- svd.Y3.sample$u[, 1:rank3]
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
  cluster <- Ckmedian.1d.dp(svd.avg$d, k=2)
  joint <- svd.avg$u[, cluster$cluster == 2, drop = FALSE]
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
compute[["proposed_subsampling2"]] <- function(Y1, Y2, Y3, rank1, rank2, rank3, numSamples=100) {
  avg.P <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
  avg.P1 <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
  avg.P2 <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
  avg.P3 <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
  for (i in 1:numSamples){
    if (i %% 3 == 0) {
      Y1.sample <- Y1[, sample(1:ncol(Y1), as.integer(ncol(Y1)/2), replace=FALSE)]
      Y2.sample <- Y2
      Y3.sample <- Y3
    } else if (i %% 3 == 1) {
      Y1.sample <- Y1
      Y2.sample <- Y2[, sample(1:ncol(Y2), as.integer(ncol(Y2)/2), replace=FALSE)]
      Y3.sample <- Y3
    } else{
      Y1.sample <- Y1
      Y2.sample <- Y2
      Y3.sample <- Y3[, sample(1:ncol(Y3), as.integer(ncol(Y3)/2), replace=FALSE)]
    }
    svd.Y1.sample <- svd(Y1.sample)
    svd.Y2.sample <- svd(Y2.sample)
    svd.Y3.sample <- svd(Y3.sample)
    u1.sample <- svd.Y1.sample$u[, 1:rank1]
    u2.sample <- svd.Y2.sample$u[, 1:rank2]
    u3.sample <- svd.Y3.sample$u[, 1:rank3]
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
  cluster <- Ckmedian.1d.dp(svd.avg$d, k=2)
  joint <- svd.avg$u[, cluster$cluster == 2, drop = FALSE]
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
for (model in c("ajive", "proposed", "proposed_subsampling1", "proposed_subsampling2")) {
  results[[model]] <- matrix(0, sim_iter, 14)
}
for (i in 1:sim_iter) {
  d1 <- svd(matrix(signal_strength1 * rnorm(m * n1), m, n1))$d[1:(rj+ri1)]
  d2 <- svd(matrix(signal_strength2 * rnorm(m * n2), m, n2))$d[1:(rj+ri2)]
  d3 <- svd(matrix(signal_strength3 * rnorm(m * n3), m, n3))$d[1:(rj+ri3)]
  d1 <- sample(d1)
  d2 <- sample(d2)
  d3 <- sample(d3)
  dj1 <- d1[1:rj]
  dj2 <- d2[1:rj]
  dj3 <- d3[1:rj]
  di1 <- d1[(rj+1):(rj+ri1)]
  di2 <- d2[(rj+1):(rj+ri2)]
  di3 <- d3[(rj+1):(rj+ri3)]
  # generate data
  U <- svd(matrix(rnorm(m * (n1+n2+n3)), m, n1+n2+n3))$u
  # joint part
  Uj <- U[, 1:rj]
  O1 <- randortho(rj, type = 'orthonormal') # rotate Uj
  O2 <- randortho(rj, type = 'orthonormal') # rotate Uj
  Uj1 <- Uj 
  Uj2 <- Uj %*% O1
  Uj3 <- Uj %*% O2
  # individual part
  Ui1 <- U[, (rj+1):(rj+ri1)]
  Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
  Ui3 <- U[, (rj+ri1+ri2+1):(rj+ri1+ri2+ri3)]
  O1 <- matrix(runif(ri1 * ri2, -phi_max, phi_max), ri1, ri2)
  Ui2 <- Ui2 + Ui1 %*% O1
  Ui2 <- gramSchmidt(Ui2)$Q # orthonormalize
  O2 <- matrix(runif(ri2 * ri3, -phi_max, phi_max), ri2, ri3) # rotate
  Ui3 <- Ui3 + Ui2 %*% O2
  Ui3 <- gramSchmidt(Ui3)$Q # orthonormalize
  # loadings
  Vj1 <- svd(matrix(rnorm(m * n1), m, n1))$v[, 1:rj]
  Vj2 <- svd(matrix(rnorm(m * n2), m, n2))$v[, 1:rj]
  Vj3 <- svd(matrix(rnorm(m * n3), m, n3))$v[, 1:rj]
  Vi1 <- svd(matrix(rnorm(m * n1), m, n1))$v[, 1:ri1]
  Vi2 <- svd(matrix(rnorm(m * n2), m, n2))$v[, 1:ri2]
  Vi3 <- svd(matrix(rnorm(m * n3), m, n3))$v[, 1:ri3]
  # combine
  X1 <- Uj1 %*% diag(dj1) %*% t(Vj1) + Ui1 %*% diag(di1) %*% t(Vi1)
  X2 <- Uj2 %*% diag(dj2) %*% t(Vj2) + Ui2 %*% diag(di2) %*% t(Vi2)
  X3 <- Uj3 %*% diag(dj3) %*% t(Vj3) + Ui3 %*% diag(di3) %*% t(Vi3)
  # noise
  Z1 <- matrix(rnorm(m * n1), m, n1) * sigma1
  Z2 <- matrix(rnorm(m * n2), m, n2) * sigma2
  Z3 <- matrix(rnorm(m * n3), m, n3) * sigma3
  Y1 <- X1 + Z1
  Y2 <- X2 + Z2
  Y3 <- X3 + Z3
  
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
  error2 <- 0
  error3 <- 0
  if (rank_spec == 'over') {
    error1 <- sample(1:2, 1)
    error2 <- sample(1:2, 1)
    error2 <- sample(1:2, 1)
  } else if (rank_spec == 'under') {
    error1 <- (-1) * sample(1:2, 1)
    error2 <- (-1) * sample(1:2, 1)
    error3 <- (-1) * sample(1:2, 1)
  }
  rank1 <- rj + ri1 + error1
  rank2 <- rj + ri2 + error2
  rank3 <- rj + ri3 + error3
  
  # compute results
  for (model in c("ajive", "proposed", "proposed_subsampling1", "proposed_subsampling2")) {
    out <- compute[[model]](Y1, Y2, Y3, rank1, rank2, rank3)
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
results.save[["m"]] <- m
results.save[["phi_max"]] <- phi_max
results.save[["n1"]] <- n1
results.save[["n2"]] <- n2
results.save[["n3"]] <- n3
results.save[["SNR1"]] <- signal_strength1 / sigma1
results.save[["SNR2"]] <- signal_strength2 / sigma2 
results.save[["SNR3"]] <- signal_strength3 / sigma3
results.save[["sigma1"]] <- sigma1
results.save[["sigma2"]] <- sigma2
results.save[["sigma3"]] <- sigma3
results.save[["signal_strength1"]] <- signal_strength1
results.save[["signal_strength2"]] <- signal_strength2
results.save[["signal_strength3"]] <- signal_strength3
results.save[["rank_spec"]] <- rank_spec
results.save[["no_joint"]] <- no_joint
results.save[["no_indiv"]] <- no_indiv
save(results.save, file=paste0("results/demo3_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData"))

