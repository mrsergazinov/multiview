source('src/utils.R')
form_output <- function(joint, indiv1, indiv2, dims) {
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
    return (
      list(
        "joint" = joint, 
        "indiv1" = indiv1, 
        "indiv2" = indiv2
        )
      )
  }
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
jive_func <- function(Y1, Y2, rank1, rank2, return_scores=FALSE) {
  out <- jive(list(t(Y1), t(Y2)),
              method='perm', showProgress=FALSE, orthIndiv=FALSE)
  check_null <- function(X, rank){
    if (rank == 0){ return (NULL) }
    return (X[, 1:rank, drop = FALSE])
  }
  joint <- check_null(svd(t(out$joint[[1]]))$u, out$rankJ)
  indiv1 <- check_null(svd(t(out$individual[[1]]))$u, out$rankA[1])
  indiv2 <- check_null(svd(t(out$individual[[2]]))$u, out$rankA[2])
  if (return_scores) {
    return (
      list(
        "joint" = joint, 
        "indiv1" = indiv1, 
        "indiv2" = indiv2
        )
      )
  }
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
dcca_func <- function(Y1, Y2, rank1, rank2, return_scores=FALSE) {
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
  if (return_scores) {
    return (
      list(
        "joint" = joint, 
        "indiv1" = indiv1, 
        "indiv2" = indiv2
        )
      )
  }
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
slide_func <- function(Y1, Y2, rank1, rank2, return_scores=FALSE) {
  out <- slide(cbind(Y1, Y2), pvec = c(ncol(Y1), ncol(Y2)))
  check_null <- function(X, mask){
    if (any(mask)){return(X[, mask, drop = FALSE])}
    return (NULL)
  }
  joint <- check_null(out$model$U, (out$S[1, ] == 1 & out$S[2, ] == 1))
  indiv1 <- check_null(out$model$U, (out$S[1, ] == 1 & out$S[2, ] == 0))
  indiv2 <- check_null(out$model$U, (out$S[1, ] == 0 & out$S[2, ] == 1))
  if (return_scores) {
    return (
      list(
        "joint" = joint, 
        "indiv1" = indiv1, 
        "indiv2" = indiv2
        )
      )
  }
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
unifac_func <- function(Y1, Y2, rank1, rank2, return_scores=FALSE) {
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
  if (return_scores) {
    return (
      list(
        "joint" = joint, 
        "indiv1" = indiv1, 
        "indiv2" = indiv2
        )
      )
  }
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}
global_null_2_views <- function(Y1, Y2, rank1, rank2, compute_prod = TRUE) {
  m <- nrow(Y1)
  q1 <- rank1 / m
  q2 <- rank2 / m
  q.plus <- q1 + q2 - 2*q1*q2 + 2*sqrt(q1*q2*(1-q1)*(1-q2))
  lam <- sqrt(q.plus)
  # compute product of projections
  if (compute_prod) {
    svd.Y1 <- svd(Y1)
    svd.Y2 <- svd(Y2)
    U1.hat <- svd.Y1$u[, 1:rank1, drop = FALSE]
    U2.hat <- svd.Y2$u[, 1:rank2, drop = FALSE]
    P.hat <- U1.hat %*% t(U1.hat)
    Q.hat <- U2.hat %*% t(U2.hat)
    prod <- P.hat %*% Q.hat
    prod.sym <- (prod + t(prod)) / 2
    svd.prod <- svd(prod)
    svd.prod.sym <- svd(prod.sym)
    return (list("noJoint" = (svd.prod$d[1] < lam),
                 "P.hat" = P.hat,
                 "Q.hat" = Q.hat,
                 "prod" = prod,
                 "prod.sym" = prod.sym,
                 "svd.prod" = svd.prod,
                 "svd.prod.sym" = svd.prod.sym,
                 "svd.Y1" = svd.Y1,
                 "svd.Y2" = svd.Y2,
                 "lam" = lam))
  } else {
    return (lam)
  }
}
est.sigma <- function(Y){ 
  # from Gavish and Donoho 2014
  sing.vals <- svd(Y)$d  
  med.sing.val <- median(sing.vals)
  # median from Marchenko-Pastur
  med.mp <- qmp(0.5, ndf = ncol(Y), pdim = nrow(Y))
  return (med.sing.val / sqrt(med.mp * ncol(Y)))
}
bootstrap.epsilon_2_views <- function(Y1, Y2, rank1, rank2, prod.spectrum, num_iter = 100) {
  print(rank2)
  # estimate rank
  shrink.Y1 <- optishrink(Y1)
  shrink.Y2 <- optishrink(Y2)
  
  # compute residual
  svd.Y1 <- svd(Y1)
  X1.hat <- svd.Y1$u[, 1:rank1, drop = FALSE] %*% diag(svd.Y1$d[1:rank1], nrow = rank1) %*% t(svd.Y1$v[, 1:rank1, drop = FALSE])
  svd.Y2 <- svd(Y2)
  print(rank2)
  print(dim(svd.Y2$u[, 1:rank2, drop = FALSE]))
  print(dim(diag(svd.Y2$d[1:rank2], nrow = rank2)))
  X2.hat <- svd.Y2$u[, 1:rank2, drop = FALSE] %*% diag(svd.Y2$d[1:rank2], nrow = rank2) %*% t(svd.Y2$v[, 1:rank2, drop = FALSE])
  
  E1.hat <- Y1 - X1.hat
  E2.hat <- Y2 - X2.hat
  
  # impute
  U1.noise <- shrink.Y1$low.rank$u[, 1:rank1, drop=FALSE]
  V1.noise <- shrink.Y1$low.rank$v[, 1:rank1, drop=FALSE]
  U2.noise <- shrink.Y2$low.rank$u[, 1:rank2, drop=FALSE]
  V2.noise <- shrink.Y2$low.rank$v[, 1:rank2, drop=FALSE]
  D1 <- rmp(rank1, svr = min(nrow(Y1), ncol(Y1)) / max(nrow(Y1), ncol(Y1)))
  D2 <- rmp(rank2, svr = min(nrow(Y2), ncol(Y2)) / max(nrow(Y2), ncol(Y2)))
  
  est.sigma.Y1 <- est.sigma(Y1)
  est.sigma.Y2 <- est.sigma(Y2)
  E1.hat <- E1.hat + est.sigma.Y1 * U1.noise %*% diag(D1, nrow=rank1) %*% t(V1.noise)
  E2.hat <- E2.hat + est.sigma.Y2 * U2.noise %*% diag(D2, nrow=rank2) %*% t(V2.noise)
  
  # resampling
  epsilon_1s <- c()
  epsilon_2s <- c()
  d <- acos(prod.spectrum)
  cos.d <- cos(d)
  sin.d <- sin(d)
  for (i in 1:num_iter) {
    # resample signal
    U <- svd(matrix(rnorm(nrow(Y1) * (rank1+rank2)), nrow = nrow(Y1), ncol = rank1+rank2))$u[, 1:(rank1+rank2)]
    V1 <- svd(matrix(rnorm(ncol(Y1) * rank1), nrow = ncol(Y1), ncol = rank1))$u[, 1:rank1]
    V2 <- svd(matrix(rnorm(ncol(Y2) * rank2), nrow = ncol(Y2), ncol = rank2))$u[, 1:rank2]
    # align basis
    U1 <- U[, 1:rank1, drop=FALSE]
    U2 <- U[, (rank1+1):(rank1+rank2), drop=FALSE]
    rank.rot <- min(rank1, rank2)
    if (rank1 >= rank2) {
      U2 <- svd(U1[, 1:rank.rot, drop=FALSE] %*% diag(cos.d, nrow = rank.rot) + U2 %*% diag(sin.d, nrow=rank.rot))$u
    } else {
      U1 <- svd(U2[, 1:rank.rot, drop=FALSE] %*% diag(cos.d, nrow = rank.rot) + U1 %*% diag(sin.d, nrow=rank.rot))$u
    }
    # form signal
    X1 <- U1 %*% diag(shrink.Y1$low.rank$d[1:rank1], nrow=rank1) %*% t(V1)
    X2 <- U2 %*% diag(shrink.Y2$low.rank$d[1:rank2], nrow=rank2) %*% t(V2)
    # form signal + noise
    Y1.resample <- X1 + E1.hat
    Y2.resample <- X2 + E2.hat
    
    # estimate column space
    U1.hat <- svd(Y1.resample)$u[, 1:rank1, drop=FALSE]
    U2.hat <- svd(Y2.resample)$u[, 1:rank2, drop=FALSE]
    
    # form projections
    P1 <- U1 %*% t(U1)
    P2 <- U2 %*% t(U2)
    P1.hat <- U1.hat %*% t(U1.hat)
    P2.hat <- U2.hat %*% t(U2.hat)
    
    # Adjustments
    P1E2 <-  P1 %*% (P2.hat - P2)
    E1P2 <- (P1.hat - P1) %*% P2
    E1E2 <- (P1.hat - P1) %*% (P2.hat - P2)
    P1E1P2 <- P1 %*% E1P2
    P1E2P2 <- P1E2 %*% P2
    P1E1E2P2 <- P1 %*% E1E2 %*% P2
    epsilon_1 <- svd(P1E1P2 + P1E2P2 + P1E1E2P2)$d[1]
    epsilon_2 <- svd(P1E2 + E1P2 + E1E2)$d[1]
    epsilon_1s <- c(epsilon_1s, epsilon_1)
    epsilon_2s <- c(epsilon_2s, epsilon_2)
  }
  return (list("epsilon1" = epsilon_1s, "epsilon2" = epsilon_2s))
}
proposed_func <- function(Y1, Y2, rank1, rank2, rank_joint = NULL, bootstrap_iters = 100, return_scores=FALSE) {
  print(rank2)
  out <- global_null_2_views(Y1, Y2, rank1, rank2)
  epsilon <- NULL
  if (!is.null(bootstrap_iters)) {
    bootstrap <- bootstrap.epsilon_2_views(Y1, Y2, 
                                           rank1, rank2, out$svd.prod$d[1:min(rank1, rank2)], 
                                           num_iter = bootstrap_iters)
    # bound <- max(1-mean(epsilon), out$lam)
    # bound <- 1-mean(bootstrap$epsilon1)
    bound <- max(1-mean(bootstrap$epsilon1), mean(bootstrap$epsilon2))
  }
  
  if (out$noJoint) {
    joint <- NULL
    jointPerp <- diag(nrow(Y1))
  } else {
    if (is.null(rank_joint)) {
      rank_joint <- sum(out$svd.prod$d > bound)
    }
    joint <- svd(out$prod.sym)$u[, 1:rank_joint, drop = FALSE]
    jointPerp <- diag(nrow(joint)) - joint %*% t(joint)
  }
  
  svd.P.hat <- svd(jointPerp %*% out$P.hat)
  cluster <- Ckmedian.1d.dp(svd.P.hat$d, k = 2)
  indiv1 <- svd.P.hat$u[, cluster$cluster == 2, drop = FALSE]
  
  svd.Q.hat <- svd(jointPerp %*% out$Q.hat)
  cluster <- Ckmedian.1d.dp(svd.Q.hat$d, k = 2)
  indiv2 <- svd.Q.hat$u[, cluster$cluster == 2, drop = FALSE]
  
  if (return_scores) {
    return (list("joint" = joint, 
                 "indiv1" = indiv1, 
                 "indiv2" = indiv2, 
                 "test" = out,
                 "bootstrap" = bootstrap,
                 "svd.indiv1" = svd.P.hat,
                 "svd.indiv2" = svd.Q.hat))
  }
  return (form_output(joint, indiv1, indiv2, nrow(Y1)))
}