form_output <- function(joint, indiv1, indiv2, dims) {
  # Simplified function to handle null checks and calculations
  check_and_compute <- function(X, Y = NULL) {
    # If both X and Y are NULL, return infinity and a zero matrix
    if (is.null(X) && is.null(Y)) {
      return(list("r" = Inf, "P" = matrix(0, nrow = dims, ncol = dims)))
    }
    # Calculate for single or combined inputs
    compute_P <- function(Z) {
      if (is.null(Z)) matrix(0, nrow = dims, ncol = dims) else Z %*% t(Z)
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
ajive_func <- function(Y1, Y2, rank1, rank2, return_scores=FALSE){
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
  if (return_scores) {
    return (list("joint" = joint, "indiv1" = indiv1, "indiv2" = indiv2))
  }
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
jive_func <- function(Y1, Y2, rank1, rank2) {
  out <- jive(list(t(Y1), t(Y2)), rankA = c(rank1, rank2),
              method='perm', showProgress=FALSE)
  check_null <- function(X, rank){
    if (rank == 0){ return (NULL) }
    return (X[, 1:rank, drop = FALSE])
  }
  joint <- check_null(svd(t(out$joint[[1]]))$u, out$rankJ)
  indiv1 <- check_null(svd(t(out$individual[[1]]))$u, out$rankA[1])
  indiv2 <- check_null(svd(t(out$individual[[2]]))$u, out$rankA[2])
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
dcca_func <- function(Y1, Y2, rank1, rank2) {
  use_virtualenv("./multiview_pylibs", required = TRUE)
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
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
slide_func <- function(Y1, Y2, rank1, rank2) {
  out <- slide(cbind(Y1, Y2), pvec = c(ncol(Y1), ncol(Y2)))
  check_null <- function(X, mask){
    if (any(mask)){return(X[, mask, drop = FALSE])}
    return (NULL)
  }
  joint <- check_null(out$model$U, (out$S[1, ] == 1 & out$S[2, ] == 1))
  indiv1 <- check_null(out$model$U, (out$S[1, ] == 1 & out$S[2, ] == 0))
  indiv2 <- check_null(out$model$U, (out$S[1, ] == 0 & out$S[2, ] == 1))
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
unifac_func <- function(Y1, Y2, rank1, rank2) {
  source('src/unifac.plus.given.R')
  Y <- cbind(Y1, Y2)
  m <- nrow(Y)
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
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
global_null <- function(Y1, Y2, rank1, rank2, compute_prod = TRUE) {
  m <- nrow(Y1)
  q1 <- rank1 / m
  q2 <- rank2 / m
  q.plus <- q1 + q2 - 2*q1*q2 + 2*sqrt(q1*q2*(1-q1)*(1-q2))
  q.minus <- q1 + q2 - 2*q1*q2 - 2*sqrt(q1*q2*(1-q1)*(1-q2))
  if (q.minus == 0) {
    q.minus <- 1e-6 # offset to avoid comp error
  }
  A0 <- 1 - min(q1, q2) 
  pdf <- function(x) {
    sqrt((q.plus - x)*(x - q.minus)) /(2 * pi * x * (1-x))
  }
  cdf <- function(lam) {
    A0 + integrate(pdf, q.minus, lam)$value - 0.95
  }
  # find closest to 0.95
  lam <- 0
  if (A0 > 0.95) {
    lam <- 0
  } else {
    lam <- uniroot(cdf, c(q.minus, q.plus))$root
  }
  # compute product of projections
  if (compute_prod) {
    U1.hat <- svd(Y1)$u[, 1:rank1, drop = FALSE]
    U2.hat <- svd(Y2)$u[, 1:rank2, drop = FALSE]
    P.hat <- U1.hat %*% t(U1.hat)
    Q.hat <- U2.hat %*% t(U2.hat)
    prod <- P.hat %*% Q.hat
    prod.sym <- (prod + t(prod)) / 2
    svd.prod <- svd(prod)
    return (list("noJoint" = (svd.prod$d[1] < lam),
                 "P.hat" = P.hat,
                 "Q.hat" = Q.hat,
                 "prod.sym" = prod.sym,
                 "svd.prod" = svd.prod))
  } else {
    return (lam)
  }
  
}
proposed_func <- function(Y1, Y2, rank1, rank2) {
  out <- global_null(Y1, Y2, rank1, rank2)
  
  if (out$noJoint) {
    joint <- NULL
    jointPerp <- diag(nrow(prod))
  } else {
    cluster <- Ckmedian.1d.dp(sqrt(out$svd.prod$d), k = 3)
    joint <- svd(out$prod.sym)$u[, cluster$cluster == 3, drop = FALSE]
    jointPerp <- diag(nrow(joint)) - joint %*% t(joint)
  }
  
  P.hat <- jointPerp %*% out$P.hat
  svd.P.hat <- svd(P.hat)
  cluster <- Ckmedian.1d.dp(svd.P.hat$d, k = 2)
  indiv1 <- svd.P.hat$u[, cluster$cluster == 2, drop = FALSE]
  
  Q.hat <- jointPerp %*% out$Q.hat
  svd.Q.hat <- svd(Q.hat)
  cluster <- Ckmedian.1d.dp(svd.Q.hat$d, k = 2)
  indiv2 <- svd.Q.hat$u[, cluster$cluster == 2, drop = FALSE]
  
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
proposed_subsampling_func <- function(Y1, Y2, rank1, rank2, numSamples=300, return_scores=FALSE) {
  out <- global_null(Y1, Y2, rank1, rank2)
    
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
  if (out$noJoint) {
    joint <- NULL
    jointPerp <- diag(nrow(avg.P))
  } else {
    cluster <- Ckmedian.1d.dp(out$svd.prod$d, k=3)
    joint <- svd.avg$u[, cluster$cluster == 3, drop = FALSE]
    jointPerp <- diag(nrow(joint)) - joint %*% t(joint)
  }
  
  svd.avg1 <- svd(jointPerp %*% avg.P1)
  cluster <- Ckmedian.1d.dp(svd(jointPerp %*% out$P.hat)$d, k = 2)
  indiv1 <- svd.avg1$u[, cluster$cluster == 2, drop = FALSE]
  
  svd.avg2 <- svd(jointPerp %*% avg.P2)
  cluster <- Ckmedian.1d.dp(svd(jointPerp %*% out$Q.hat)$d, k = 2)
  indiv2 <- svd.avg2$u[,  cluster$cluster == 2, drop = FALSE]
  
  if (return_scores) {
    return (list("joint" = joint, 
                 "indiv1" = indiv1, 
                 "indiv2" = indiv2, 
                 "test" = out))
  }
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}