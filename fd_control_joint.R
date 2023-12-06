thresh <- function(X, sigma = NA) {
  # Based on ...
  if (is.na(sigma)){
    sigma = median(svd(X)$d) / sqrt(qmp(0.5, nrow(X), ncol(X)) * ncol(X))
  }
  return (sigma * (1 + sqrt(min(ncol(X), nrow(X)) / max(ncol(X), nrow(X)))))
}

fd_control_joint <- function(X1, X2, args){
  avg.P <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  avg.P1 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  avg.P2 <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  for (i in 1:args$numSamples){
    X1.sample <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
    X2.sample <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
    svd.X1.sample <- svd(X1.sample)
    svd.X2.sample <- svd(X2.sample)
    thresh.X1 <- thresh(X1.sample, sigma=args$sigma1) # thresholding singular values
    thresh.X2 <- thresh(X2.sample, sigma=args$sigma2)
    u1.sample <- svd.X1.sample$u[, svd.X1.sample$d > thresh.X1 + 1e-10]
    u2.sample <- svd.X2.sample$u[, svd.X2.sample$d > thresh.X2 + 1e-10]
    sample.P1 <- (u1.sample %*% t(u1.sample))
    sample.P2 <- (u2.sample %*% t(u2.sample))
    avg.P1 <- avg.P1 + sample.P1
    avg.P2 <- avg.P2 + sample.P2
    avg.P <- avg.P + (sample.P1 %*% sample.P2)
    # avg.P <- avg.P + (sample.P1 %*% sample.P2 %*% sample.P1 + 
                        # sample.P2 %*% sample.P1 %*% sample.P2) / 2
    # avg.P <- avg.P + (sample.P1 + sample.P2) / 2
  }
  avg.P1 <- avg.P1 / args$numSamples
  avg.P2 <- avg.P2 / args$numSamples
  avg.P <- avg.P / args$numSamples
  
  # avg.P <- avg.P1 %*% avg.P2
  svd.avg <- svd(avg.P)
  joint <- svd.avg$u[, svd.avg$d > args$alpha, drop = FALSE]
  jointPerp <- diag(nrow(joint)) - joint %*% t(joint)

  avg.P1 <- jointPerp %*% avg.P1
  svd.avg1 <- svd(avg.P1)
  indiv1 <- svd.avg1$u[, svd.avg1$d > args$alpha, drop = FALSE]

  avg.P2 <- jointPerp %*% avg.P2
  svd.avg2 <- svd(avg.P2)
  indiv2 <- svd.avg2$u[, svd.avg2$d > args$alpha, drop = FALSE]
  
  return(list("joint" = joint, "indiv1" = indiv1, "indiv2" = indiv2, "avgPd" = svd.avg$d))
}


# avg.P1 <- avg.P1 / args$numSamples
# avg.P2 <- avg.P2 / args$numSamples
# svd.avg1 <- svd(avg.P1)
# svd.avg2 <- svd(avg.P2)
# signal1 <- svd.avg1$u[, svd.avg1$d > args$alpha, drop = FALSE]
# signal2 <- svd.avg2$u[, svd.avg2$d > args$alpha, drop = FALSE]
# svd.avg <- svd(signal1 %*% t(signal1) %*% signal2 %*% t(signal2))
# joint <- svd.avg$u[, svd.avg$d > args$alpha, drop = FALSE]
# jointPerp <- diag(nrow(joint)) - joint %*% t(joint)

# thresh <- function(X, sigma = NA){
#   # thresholding using Gavish and Donoho 2014
#   beta = nrow(X) / ncol(X)
#   if (is.na(sigma)){
#     sigma = median(svd(X)$d)
#     omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
#     return (omega * sigma)
# 
#   } else {
#     lambda = sqrt(2*(beta+1)+8*beta/((beta+1) + sqrt(beta^2+14*beta+1)))
#     return (lambda * sigma * sqrt(ncol(X)))
#   }
# }

# test different alphas
# results = list()
# for (alpha in seq(0.5, 0.95, 0.05)){
#   fds <- rep(0, sim_iter)
#   for (i in 1:sim_iter){
#     data <- gen_data(n, p1, p2, rj, ri1, ri2, dj, di1, di2, sigma1, sigma2)
#     X1 <- data[["X1"]]
#     X2 <- data[["X2"]]
#     Uj <- data[["joint"]]
#     
#     # apply oracle ajive
#     jointFdControlOut <- jointFdControl(X1, X2, alpha = alpha)
#     colProjJointFdControl <- jointFdControlOut %*% t(jointFdControlOut)
#     # compute true col space
#     trueProj <- Uj %*% t(Uj)
#     
#     # compute FD and metric
#     fds[i] <- sum(diag((diag(n)-trueProj) %*% colProjJointFdControl))
#   }
#   results[[as.character(alpha)]] <- mean(fds)
# }
