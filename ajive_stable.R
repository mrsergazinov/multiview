# stableAJIVE <- function(X1, X2, numSamples=100, alpha=0.7) {
#   avgP <-  matrix(0, nrow=nrow(X1), ncol=nrow(X1))
#   for (i in 1:numSamples){
#     X1.sample <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
#     X2.sample <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
#     ajiveOut <- ajive(list(X1.sample, X2.sample), c(10, 10), 
#                       n_wedin_samples = 100, n_rand_dir_samples = 100) # input -> true joint rank
#     P <- ajiveOut$joint_scores %*% t(ajiveOut$joint_scores)
#     if (is.na(ajiveOut$joint_scores)){
#       print(i)
#       print('null space')
#       P <- matrix(0, nrow=n, ncol=n)
#     }
#     avgP <- avgP + P
#   }
#   avgP <- avgP / numSamples
#   svdAvg <- svd(avgP)
#   return (svdAvg$u[, svdAvg$d > alpha])
# }