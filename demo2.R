suppressPackageStartupMessages({
  library(tidyverse)
  library(ajive)
  library(r.jive)
  library(SLIDE)
  library(RMTstat)
  library(pracma)
  library(Ckmeans.1d.dp)
  library(foreach)
  library(doParallel)
})

# define number of cores and start parallel backend
set.seed(1234)
numCores <- detectCores() - 1  # Leave one core for system processes
cl <- makeCluster(numCores)
registerDoParallel(cl)

# define other params
rj <- 4
ri1 <- 3
ri2 <- 2
m <- 50
phi_max <- 0.0
n1 <- 80
n2 <- 100
snr1 <- 4
snr2 <- 4
signal_strength1 <- 10
signal_strength2 <- 12
sigma1 <- (signal_strength1 / snr1) / (sqrt(m) + sqrt(n1))
sigma2 <- (signal_strength2 / snr2) / (sqrt(m) + sqrt(n2))
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

sim_iter <- 50
models <- c("jive", "ajive", "dcca", "slide", "unifac", "proposed", "proposed_subsampling")
iters <- foreach(i = 1:sim_iter,
                   .packages=c("pracma", "r.jive", "ajive", "SLIDE", "Ckmeans.1d.dp")) %dopar% {
  # compute error
  compute_fd <- function(P1, P2){
    sum(diag((diag(dim(P2)[1]) - P2) %*% P1))
  }
  compute_tp <- function(P1, P2){
    sum(diag(P1 %*% P2))
  }
  # models
  compute <- list()
  form_output <- function(joint, indiv1, indiv2) {
    # Simplified function to handle null checks and calculations
    check_and_compute <- function(X, Y = NULL) {
      # If both X and Y are NULL, return infinity and a zero matrix
      if (is.null(X) && is.null(Y)) {
        return(list("r" = Inf, "P" = matrix(0, nrow = m, ncol = m)))
      }
      # Calculate for single or combined inputs
      compute_P <- function(Z) {
        if (is.null(Z)) matrix(0, nrow = m, ncol = m) else Z %*% t(Z)
      }
      P_X <- compute_P(X)
      P_Y <- compute_P(Y)
      return(list("r" = ifelse(is.null(X), 0, ncol(X)) + ifelse(is.null(Y), 0, ncol(Y)),
                  "P" = P_X + P_Y))
    }
    
    # Apply the function to each combination of inputs
    Pjoint <- check_and_compute(joint)
    Pindiv1 <- check_and_compute(indiv1)
    Pindiv2 <- check_and_compute(indiv2)
    P1 <- check_and_compute(joint, indiv1)
    P2 <- check_and_compute(joint, indiv2)
    
    # Return the results in a list
    return(list("P1" = P1$P, "P2" = P2$P,
                "Pjoint" = Pjoint$P, "Pindiv1" = Pindiv1$P, "Pindiv2" = Pindiv2$P,
                "r1" = P1$r, "r2" = P2$r,
                "rj" = Pjoint$r, "ri1" = Pindiv1$r, "ri2" = Pindiv2$r))
  }
  compute[["ajive"]] <- function(Y1, Y2, rank1, rank2){
    out <- ajive(list(Y1, Y2), c(rank1, rank2),
                n_wedin_samples = 100, 
                n_rand_dir_samples = 100)
    check_null <- function(X) {
      if (any(is.na(X))) {
        return (NULL)
      } else if (all(X == 0)){
        return (NULL)
      }
      return (X)
    }
    joint <- check_null(out$joint_scores)
    indiv1 <- check_null(out$block_decomps[[1]][['individual']][['u']])
    indiv2 <- check_null(out$block_decomps[[2]][['individual']][['u']])
    return (form_output(joint, indiv1, indiv2))
  }
  compute[["jive"]] <- function(Y1, Y2, rank1, rank2) {
    out <- jive(list(t(Y1), t(Y2)), rankA = c(rank1, rank2),
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
  compute[["dcca"]] <- function(Y1, Y2, rank1, rank2) {
    library(reticulate)
    use_condaenv("base", required = TRUE)
    source_python("src/dcca.py")
    res <- dCCA(t(Y1), t(Y2), r_1 = as.integer(rank1), r_2 = as.integer(rank2))
    
    #--------------------------------------------------------------
    # Below code is taken from Gaynanova et al.
    #--------------------------------------------------------------
    projection <- function(A){
      return(A %*% t(A))
    }
    angle_cal <- function(X, Y, tol = .Machine$double.eps^0.5){
      X = as.matrix(X)
      Y = as.matrix(Y)
      Xnorm = svd(X)$u
      Ynorm = svd(Y)$u
      M = crossprod(Xnorm,Ynorm)
      # Extreme case when both X and Y only contain one number
      if (dim(M)[1] == 1 && dim(M)[2] == 1){
        cos_angle = abs(M)
        principal_angle = NA
        if (cos_angle >= 1){principal_angle = 0}
        if (cos_angle <= 0){principal_angle = pi/2}
        if (cos_angle > 0 && cos_angle < 1){
          principal_angle = acos(cos_angle)
        }
        principal_mat1 = Xnorm 
        principal_mat2 = Ynorm 
        return(list("angle" = principal_angle, "cos_angle" = cos_angle,
                    "principal_vector1" = principal_mat1, 
                    "principal_vector2" = principal_mat2))
      }
      # Normal case when X and Y are matrices (data frames)
      else{
        svd_result = svd(M)
        cos_angle = svd_result$d
        l = length(cos_angle)
        principal_angle = rep(NA, l)
        for (i in 1:l){
          if (cos_angle[i] >= 1){principal_angle[i] = 0}
          if (cos_angle[i] <= 0){principal_angle[i] = pi/2}
          if (cos_angle[i] > 0 && cos_angle[i] < 1){
            principal_angle[i] = acos(cos_angle[i])
          }
        }
        principal_mat1 = Xnorm %*% svd_result$u
        principal_mat2 = Ynorm %*% svd_result$v
        return(list("angle" = principal_angle, "cos_angle" = cos_angle,
                    "principal_vector1" = principal_mat1, 
                    "principal_vector2" = principal_mat2))
      }
    }
    # First two elements of the lists are fitted thetas
    theta1hat = t(res[[1]])
    theta2hat = t(res[[2]])
    # Extract joint structure, 3 and 4 are Cs, 5 and 6 are Ds
    # Ranks are 7, 8, 9
    r0 = res[[9]]
    r1 = res[[7]] - r0
    r2 = res[[8]] - r0
    # DCCA joint from Dongbang
    col1_dcca = svd(theta1hat)$u[,1:(r0+r1)]
    col2_dcca = svd(theta2hat)$u[,1:(r0+r2)]
    angle_result_dcca = angle_cal(col1_dcca, col2_dcca)
    # Get DCCA joint canonical variables
    j1hat = angle_result_dcca$principal_vector1[,1:r0, drop = FALSE]
    j2hat = angle_result_dcca$principal_vector2[,1:r0, drop = FALSE]
    # Get DCCA individual canonical variables
    i1hat = svd(theta1hat - projection(j1hat) %*% theta1hat)$u[, 1:r1, drop=FALSE]
    i2hat = svd(theta2hat - projection(j2hat) %*% theta2hat)$u[, 1:r2, drop=FALSE]
    #--------------------------------------------------------------
    
    joint <- (j1hat + j2hat) / 2
    indiv1 <- i1hat
    indiv2 <- i2hat
    return (form_output(joint, indiv1, indiv2))
  }
  compute[["slide"]] <- function(Y1, Y2, rank1, rank2) {
    out <- slide(cbind(Y1, Y2), pvec = c(ncol(Y1), ncol(Y2)))
    check_null <- function(X, mask){
      if (any(mask)){return(X[, mask, drop = FALSE])}
      return (NULL)
    }
    joint <- check_null(out$model$U, (out$S[1, ] == 1 & out$S[2, ] == 1))
    indiv1 <- check_null(out$model$U, (out$S[1, ] == 1 & out$S[2, ] == 0))
    indiv2 <- check_null(out$model$U, (out$S[1, ] == 0 & out$S[2, ] == 1))
    return (form_output(joint, indiv1, indiv2))
  }
  compute[["unifac"]] <- function(Y1, Y2, rank1, rank2) {
    source('src/unifac.plus.given.R')
    Y <- cbind(Y1, Y2)
    p.ind <- list()
    p.ind[[1]] <- 1:ncol(Y1)
    p.ind[[2]] <- (ncol(Y1)+1):ncol(Y)
    p.ind.list <- list()
    p.ind.list[[1]] <- unlist(p.ind[1:2])
    p.ind.list[[2]] <- unlist(p.ind[c(1)])
    p.ind.list[[3]] <- unlist(p.ind[c(2)])
    res.g <- unifac.plus.given(t(Y), p.ind, p.ind.list, n=m, max.iter=1000)
    Jmat.ug <- t(res.g$S[[1]])
    I1.mat.ug <- t(res.g$S[[2]])
    I2.mat.ug <- t(res.g$S[[3]])
    
    svd.joint <- svd(Jmat.ug)
    svd.I1 <- svd(I1.mat.ug)
    svd.I2 <- svd(I2.mat.ug)
    joint <- svd.joint$u[, svd.joint$d > 1e-6, drop = FALSE]
    indiv1 <- svd.I1$u[, svd.I1$d > 1e-6, drop = FALSE]
    indiv2 <- svd.I2$u[, svd.I2$d > 1e-6, drop = FALSE]
    if (length(joint) == 0) joint <- NULL
    if (length(indiv1) == 0) indiv1 <- NULL
    if (length(indiv2) == 0) indiv2 <- NULL
    return (form_output(joint, indiv1, indiv2))
  }
  compute[["proposed"]] <- function(Y1, Y2, rank1, rank2) {
    U1.hat <- svd(Y1)$u[, 1:rank1, drop = FALSE]
    U2.hat <- svd(Y2)$u[, 1:rank2, drop = FALSE]
    P.hat <- U1.hat %*% t(U1.hat)
    Q.hat <- U2.hat %*% t(U2.hat)
  
    prod <- (P.hat %*% Q.hat + Q.hat %*% P.hat) / 2
    svd.prod <- svd(prod)
    if (svd.prod$d[1] < 0.5) {
      joint <- NULL
      jointPerp <- diag(nrow(prod))
    } else {
      cluster <- Ckmedian.1d.dp(sqrt(svd.prod$d), k = 3)
      joint <- svd.prod$u[, cluster$cluster == 3, drop = FALSE]
      jointPerp <- diag(nrow(joint)) - joint %*% t(joint)
    }
  
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
  compute[["proposed_subsampling"]] <- function(Y1, Y2, rank1, rank2, numSamples=200) {
    avg.P <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
    avg.P1 <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
    avg.P2 <- matrix(0, nrow=nrow(Y1), ncol=nrow(Y1))
    for (i in 1:numSamples){
      if (i %% 2 == 0) {
        Y1.sample <- Y1[, sample(1:ncol(Y1), as.integer(ncol(Y1)/2), replace=FALSE)]
        Y2.sample <- Y2
      } else {
        Y1.sample <- Y1
        Y2.sample <- Y2[, sample(1:ncol(Y2), as.integer(ncol(Y2)/2), replace=FALSE)]
      }
      svd.Y1.sample <- svd(Y1.sample)
      svd.Y2.sample <- svd(Y2.sample)
      u1.sample <- svd.Y1.sample$u[, 1:rank1]
      u2.sample <- svd.Y2.sample$u[, 1:rank2]
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
    if (svd.avg$d[1] < 0.5) {
      joint <- NULL
      jointPerp <- diag(nrow(avg.P))
    } else {
      cluster <- Ckmedian.1d.dp(svd.avg$d, k=3)
      joint <- svd.avg$u[, cluster$cluster == 3, drop = FALSE]
      jointPerp <- diag(nrow(joint)) - joint %*% t(joint)
    }
    
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
  # data generator
  d1 <- runif(rj+ri1, min=signal_strength1, max=2*signal_strength1)
  d2 <- runif(rj+ri2, min=signal_strength2, max=2*signal_strength2)
  d1 <- sample(d1)
  d2 <- sample(d2)
  dj1 <- d1[1:rj]
  dj2 <- d2[1:rj]
  di1 <- d1[(rj+1):(rj+ri1)]
  di2 <- d2[(rj+1):(rj+ri2)]
  # generate data
  U <- svd(matrix(rnorm(m * (n1+n2)), m, n1+n2))$u
  # joint parts
  Uj <- U[, 1:rj]
  O <- randortho(rj, type = 'orthonormal') # rotate Uj
  Uj1 <- Uj %*% O
  Uj2 <- Uj
  if (no_joint) {
    Uj1 <- matrix(0, m, rj)
    Uj2 <- matrix(0, m, rj)
  }
  # individual part
  Ui1 <- U[, (rj+1):(rj+ri1)]
  Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
  O <- matrix(runif(ri1 * ri2, -phi_max, phi_max), ri1, ri2) # rotate
  Ui2 <- Ui2 + Ui1 %*% O
  Ui2 <- gramSchmidt(Ui2)$Q # orthonormalize
  if (no_indiv) {
    Ui1 <- matrix(0, m, ri1)
    Ui2 <- matrix(0, m, ri2)
  }
  # loadings
  Vj1 <- svd(matrix(rnorm(m * n1), m, n1))$v[, 1:rj]
  Vj2 <- svd(matrix(rnorm(m * n2), m, n2))$v[, 1:rj]
  Vi1 <- svd(matrix(rnorm(m * n1), m, n1))$v[, 1:ri1]
  Vi2 <- svd(matrix(rnorm(m * n2), m, n2))$v[, 1:ri2]
  # combine
  X1 <- Uj1 %*% diag(dj1) %*% t(Vj1) + Ui1 %*% diag(di1) %*% t(Vi1)
  Z1 <- matrix(rnorm(m * n1), m, n1) * sigma1
  Y1 <- X1 + Z1
  X2 <- Uj2 %*% diag(dj2) %*% t(Vj2) + Ui2 %*% diag(di2) %*% t(Vi2)
  Z2 <- matrix(rnorm(m * n2), m, n2) * sigma2
  Y2 <- X2 + Z2

  # compute true
  Pjoint <- Uj1 %*% t(Uj1) # true projection
  Pindiv1 <- Ui1 %*% t(Ui1)
  Pindiv2 <- Ui2 %*% t(Ui2)
  P1 <- Pjoint + Pindiv1
  P2 <- Pjoint + Pindiv2
  
  # signal rank
  error1 <- 0
  error2 <- 0
  if (rank_spec == 'over') {
    error1 <- sample(1:2, 1)
    error2 <- sample(1:2, 1)
  } else if (rank_spec == 'under') {
    error1 <- (-1) * sample(1:2, 1)
    error2 <- (-1) * sample(1:2, 1)
  }
  if (no_joint) {
    rank1 <- ri1 + error1
    rank2 <- ri2 + error2
  } else if (no_indiv) {
    rank1 <- rj + error1
    rank2 <- rj + error2
  } else {
    rank1 <- rj + ri1 + error1
    rank2 <- rj + ri2 + error2
  }

  # compute results
  result <- list()
  for (model in models) {
    out <- compute[[model]](Y1, Y2, rank1, rank2)
    out <- c(compute_fd(out$P1, P1) / out$r1, compute_tp(out$P1, P1) / rank1,
             compute_fd(out$P2, P2) / out$r2, compute_tp(out$P2, P2) / rank2,
             compute_fd(out$Pjoint, Pjoint) / out$rj, compute_tp(out$Pjoint, Pjoint) / rj,
             compute_fd(out$Pindiv1, Pindiv1) / out$ri1, compute_tp(out$Pindiv1, Pindiv1) / ri1,
             compute_fd(out$Pindiv2, Pindiv2) / out$ri2, compute_tp(out$Pindiv2, Pindiv2) / ri2)
    # check if any NA
    if (any(is.na(out))){
      print(paste0("NA detected in ", model))
      break
    }
    result[[model]] <- matrix(out, nrow=1)
  }
  result
}
# shut down the parallel backend
stopCluster(cl)
# collect results
results <- list()
for (iter in iters) {
  for (model in models) {
    results[[model]] = rbind(results[[model]], iter[[model]])
  }
}
# process and save
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
results.save[["m"]] <- m
results.save[["phi_max"]] <- phi_max
results.save[["n1"]] <- n1
results.save[["n2"]] <- n2
results.save[["SNR1"]] <- snr1
results.save[["SNR2"]] <- snr2
results.save[["sigma1"]] <- sigma1
results.save[["sigma2"]] <- sigma2
results.save[["signal_strength1"]] <- signal_strength1
results.save[["signal_strength2"]] <- signal_strength2
results.save[["rank_spec"]] <- rank_spec
results.save[["no_joint"]] <- no_joint
results.save[["no_indiv"]] <- no_indiv
save(results.save, file=paste0("results/demo2_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData"))
