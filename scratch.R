

# fdControlConcat <- function(X1, X2, numSamples=100, alpha=0.7){
#   # X is n x p matrix
#   X = cbind(X1, X2)
#   n <- nrow(X)
#   p <- ncol(X)
#   thresh.X = thresh(X) # thresholding singular valuess
#   avgP <- matrix(0, nrow=n, ncol=n)
#   for (i in 1:numSamples){
#     # resample cols
#     Xsample <- X[, sample(1:p, as.integer(p/2), replace=FALSE)]
#     svd.Xsample = svd(Xsample)
#     U <- svd.Xsample$u[, svd.Xsample$d > thresh.X] # decompose and extract left singular vectors
#     P <- U %*% t(U) # projection matrix
#     avgP <- avgP + P
#   }
#   avgP <- avgP / numSamples
#   svdAvg <- svd(avgP) 
#   return (svdAvg$u[, svdAvg$d > alpha]) # select left singular vectors for which singular values are > alpha
# }

# ajive <- function(X1, X2, selectRankJoint){
#     # X1, X2 are n x p1, n x p2 matrices
#     thresh.X1 <- thresh(X1) # thresholding singular values
#     thresh.X2 <- thresh(X2)
#     svd.X1 <- svd(X1)
#     svd.X2 <- svd(X2)
#     u.X1 <- svd.X1$u[, svd.X1$d > thresh.X1] # decompose and extract left singular vectors
#     u.X2 <- svd.X2$u[, svd.X2$d > thresh.X2]
#     
#     return (svd(cbind(u.X1, u.X2))$u[, 1:selectRankJoint])
# }

# ajiveModFD1 <- function(X1, X2, selectRankJoint){
#     # select ranks using FD at individual level
#     colProjFD.X1 <- fdControl(X1)
#     colProjFD.X2 <- fdControl(X2)

#     svdInner <- svd(crossprod(colProjFD.X1, colProjFD.X2))
#     colProjFD.X1.rotated <- colProjFD.X1 %*% svdInner$u
#     colProjFD.X2.rotated <- colProjFD.X2 %*% svdInner$v

#     return (svd(cbind(colProjFD.X1.rotated, colProjFD.X2.rotated))$u[, 1:selectRankJoint])
# }

# ajiveModFD2 <- function(X1, X2){
#   # select both ranks using FD
#   colProjFD.X1 <- fdControl(X1)
#   colProjFD.X2 <- fdControl(X2)

#   # svdInner <- svd(crossprod(colProjFD.X1, colProjFD.X2))
#   # colProjFD.X1.rotated <- colProjFD.X1 %*% svdInner$u
#   # colProjFD.X2.rotated <- colProjFD.X2 %*% svdInner$v

#   return ((fdControl(cbind(colProjFD.X1, colProjFD.X2))))
# }

# plot results using ggplot: plot CI for each method, x-axis is SNR, y-axis is fd
# results.fd$lower <- results.fd$Mean - 1.96 * results.fd$SD
# results.fd$upper <- results.fd$Mean + 1.96 * results.fd$SD
# plt <- ggplot(results.fd, aes(x = SNR, y = Mean, ymin = lower, ymax = upper, group = Method, color = Method)) +
#   geom_line() +
#   geom_ribbon(alpha = 0.3) +
#   labs(x = "SNR", y = "FD") +
#   theme_minimal()
# write.csv(results.fd, "results_fd.csv")
# ggsave("results_fd.png", plt, width = 6, height = 4, dpi = 300)

# plot results for metric
# results.metric$lower <- results.metric$Mean - 1.96 * results.metric$SD
# results.metric$upper <- results.metric$Mean + 1.96 * results.metric$SD
# plt <- ggplot(results.metric, aes(x = SNR, y = Mean, ymin = lower, ymax = upper, group = Method, color = Method)) +
#   geom_line() +
#   geom_ribbon(alpha = 0.3) +
#   labs(x = "SNR", y = "Metric") +
#   theme_minimal()
# write.csv(results.metric, "results_metric.csv")
# ggsave("results_metric.png", plt, width = 6, height = 4, dpi = 300)

# scree plot of X1, X2
# scree.X1 <- prcomp(X1, scale = TRUE)
# scree.X2 <- prcomp(X2, scale = TRUE)
# var_explained.X1 = scree.X1$sdev^2 / sum(scree.X1$sdev^2)
# var_explained.X2 = scree.X2$sdev^2 / sum(scree.X2$sdev^2)