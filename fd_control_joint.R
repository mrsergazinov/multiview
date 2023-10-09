thresh <- function(X, sigma = NA){
  # thresholding using Gavish and Donoho 2014
  beta = nrow(X) / ncol(X)
  if (is.na(sigma)){
    sigma = median(svd(X)$d)
    omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
    return (omega * sigma)

  } else {
    lambda = sqrt(2*(beta+1)+8*beta/((beta+1) + sqrt(beta^2+14*beta+1)))
    return (lambda * sigma * sqrt(ncol(X)))
  }
}

jointFdControl <- function(X1, X2, args){
  avg.P <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  for (i in 1:args$numSamples){
    X1.sample <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
    X2.sample <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
    thresh.X1 <- thresh(X1.sample, sigma=args$sigma1) # thresholding singular values
    thresh.X2 <- thresh(X2.sample, sigma=args$sigma2)
    svd.X1.sample <- svd(X1.sample)
    svd.X2.sample <- svd(X2.sample)
    u1.sample <- svd.X1.sample$u[, svd.X1.sample$d > thresh.X1] # thresholding using Gavish and Donoho 2014
    u2.sample <- svd.X2.sample$u[, svd.X2.sample$d > thresh.X2]
    avg.step <- (u1.sample %*% t(u1.sample) + u2.sample %*% t(u2.sample)) / 2
    avg.P <- avg.P + avg.step
  }
  avg.P <- avg.P / args$numSamples
  svd.avg <- svd(avg.P)
  out <- svd.avg$u[, svd.avg$d > args$alpha]
  return(out)
}

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
