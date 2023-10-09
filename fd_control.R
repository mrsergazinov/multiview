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

fdControl <- function(X1, X2, args){
    X <- X1
    # X is n x p matrix
    n <- nrow(X)
    p <- ncol(X)
    avgP <- matrix(0, nrow = n, ncol = n)
    for (i in 1:args$numSamples){
        # resample cols
        Xsample <- X[, sample(1:p, as.integer(p/2), replace=FALSE)]
        thresh.X = thresh(Xsample, sigma=args$sigma) # thresholding singular valuess
        svd.Xsample = svd(Xsample)
        U <- svd.Xsample$u[, svd.Xsample$d > thresh.X] # decompose and extract left singular vectors
        P <- U %*% t(U) # projection matrix
        avgP <- avgP + P
    }
    avgP <- avgP / args$numSamples
    svdAvg <- svd(avgP) 
    return (svdAvg$u[, svdAvg$d > args$alpha]) # select left singular vectors for which singular values are > alpha
}