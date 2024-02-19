library(tidyverse)
library(ajive)
library(RMTstat)
library(pracma)

source("gen_data.R")
source("fd_control_joint.R")
source("ajive_oracle.R")
source("ajive.R")
source("bounds.R")

set.seed(235017)
# rj <- 2
# ri1 <- 3
# ri2 <- 2
# n <- 20
# phi.max <- 0.4
# p1 <- 100
# p2 <- 100
# sigma1 <- 1
# sigma2 <- 1
# sim_iter <- 100
# signal_strength <- 40
# dj <- rnorm(rj, mean = signal_strength, sd = 3)
# di1 <- rnorm(ri1, mean = signal_strength, sd = 3)
# di2 <- rnorm(ri2, mean = signal_strength, sd = 3)
# 
# # generate data
# U <- svd(matrix(rnorm(n * (p1+p2)), n, p1+p2))$u
# # joint part
# Uj <- U[, 1:rj]
# Ujperp <- U[, (rj+1):ncol(U)]
# O <- randortho(rj, type = 'orthonormal') # rotate Uj
# Uj1 <- Uj %*% O
# Uj2 <- Uj
# # individual part
# Ui1 <- U[, (rj+1):(rj+ri1)]
# Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
# O <- matrix(runif(ri1 * ri2, -phi.max, phi.max), ri1, ri2) # rotate
# Ui2 <- Ui2 + Ui1 %*% O
# Ui2 <- gramSchmidt(Ui2)$Q # orthonormalize
# # loadings
# V1 <- svd(matrix(rnorm((rj + ri1) * p1), rj + ri1, p1))$v
# V2 <- svd(matrix(rnorm((rj + ri2) * p2), rj + ri2, p2))$v
# # combine
# X1 <- cbind(Uj1, Ui1) %*% diag(c(dj, di1)) %*% t(V1)
# X2 <- cbind(Uj2, Ui2) %*% diag(c(dj, di2)) %*% t(V2)
# # noise
# X1 <- X1 + matrix(rnorm(n * p1), n, p1) * sigma1
# X2 <- X2 + matrix(rnorm(n * p2), n, p2) * sigma2
# 
# # compute bound
# bound.val <- 0.5
# # bound.approx.angle(c(dj, di1), c(dj, di2), rj, ri1, ri2, n / p1, 0, phi.max)
# print(paste("Spectral bound = ", bound.val))
# 
# # FD control -- subsampling
# args = list("sigma1" = NA, "sigma2" = NA,
#             "rj" = rj, "ri1" = ri1,
#             "ri2" = ri2, "numSamples" = 100,
#             "alpha" = 0.4, "boundJoint" = bound.val)
# out <- fd_control_joint(X1, X2, args)
# estimJointProj.subsample <- out$joint %*% t(out$joint)
# 
# tr.nosubsample <- c()
# tr.subsample <- c()
# for (i in 1:100) {
#   # FD control -- no subsampling
#   svd.X1 <- svd(X1)
#   svd.X2 <- svd(X2)
#   thresh.X1 <- thresh(X1) # thresholding singular values
#   thresh.X2 <- thresh(X2)
#   u1 <- svd.X1$u[, svd.X1$d > thresh.X1]
#   u2 <- svd.X2$u[, svd.X2$d > thresh.X2]
#   P1 <- (u1 %*% t(u1))
#   P2 <- (u2 %*% t(u2))
#   prod <- P1 %*% P2
#   svd.prod <- svd(prod)
#   joint <- svd.prod$u[, svd.prod$d > bound.val, drop = FALSE]
#   estimJointProj.noSubsample <- joint %*% t(joint) # no subsampling
#   
#   
#   for (i in 1:ncol(Ujperp)){
#       u <- Ujperp[, i, drop = FALSE]
#       tr.nosubsample <- c(tr.nosubsample, sum(diag(estimJointProj.noSubsample %*% u %*% t(u))))
#       tr.subsample <- c(tr.subsample, sum(diag(estimJointProj.subsample %*% u %*% t(u))))
#   }
# }
# 

# # data
# n = 20
# p1 <- 100
# p2 <- 100
# args = list("sigma1" = NA, "sigma2" = NA,
#             "rj" = 1, "ri1" = 1,
#             "ri2" = 2, "numSamples" = 100,
#             "alpha" = 0.4, "boundJoint" = NA)
# blocks <- sample_toy_data(n, p1, p2, 1, only_observations = FALSE)
# 
# # extract
# X1 <- blocks$obs[[1]]
# X2 <- blocks$obs[[2]]
# Uj <- svd(blocks$decomp[[1]]$joint$full)$u[, 1:args$rj]
# Ui1 <- svd(blocks$decomp[[1]]$individual$full)$u[, 1:args$ri1]
# Ui2 <- svd(blocks$decomp[[2]]$individual$full)$u[, 1:args$ri2]
# jointProj <- Uj %*% t(Uj) # true projection
# indiv1Proj <- Ui1 %*% t(Ui1)
# indiv2Proj <- Ui2 %*% t(Ui2)


# generate data
# n <- 20
# p1 <- 100
# p2 <- 100
# U <- svd(matrix(rnorm(n * (p1+p2)), n, p1+p2))$u
# # joint part
# Uj <- U[, 1]
# # individual part
# Ui1 <- U[, 2]
# Ui2 <- U[, 3:4]
# Ui2 <- Ui2 + Ui1 * 0.6
# Ui2 <- gramSchmidt(Ui2)$Q # orthonormalize
# # loadings
# V1 <- svd(matrix(rnorm(n * p1), n, p1))$v[, 1:2]
# V2 <- svd(matrix(rnorm(n * p2), n, p2))$v[, 1:3]
# # combine
# X1 <- cbind(Uj, Ui1) %*% diag(c(30, 40)) %*% t(V1)
# X2 <- cbind(Uj, Ui2) %*% diag(c(15, 20, 25)) %*% t(V2)
# # noise
# X1 <- X1 + matrix(rnorm(n * p1), n, p1) * 1
# X2 <- X2 + matrix(rnorm(n * p2), n, p2) * 1
# 
# args = list("sigma1" = NA, "sigma2" = NA,
#             "rj" = 1, "ri1" = 1,
#             "ri2" = 2, "numSamples" = 100,
#             "alpha" = 0.4, "boundJoint" = NA)
# jointProj <- Uj %*% t(Uj) # true projection
# indiv1Proj <- Ui1 %*% t(Ui1)
# indiv2Proj <- Ui2 %*% t(Ui2)
# 
# 
# # AJIVE
# out <- ajive_wrapper(X1, X2, args)
# estimJointProj <- out$joint %*% t(out$joint) # compute projection matrix
# estimIndiv1Proj <- out$indiv1 %*% t(out$indiv1)
# estimIndiv2Proj <- out$indiv2 %*% t(out$indiv2)
# if (all(is.na(out$joint))){
#   print(paste(c('joint is null space ->', model_name), collapse = " "))
#   estimJointProj <- matrix(0, nrow=n, ncol=n)
# }
# if (all(is.na(out$indiv1))){
#   print(paste(c('indiv1 is null space ->', model_name), collapse = " "))
#   estimIndiv1Proj <- matrix(0, nrow=n, ncol=n)
# }
# if (all(is.na(out$indiv2))){
#   print(paste(c('indiv2 is null space ->', model_name), collapse = " "))
#   estimIndiv2Proj <- matrix(0, nrow=n, ncol=n)
# }
# 
# results_ajive <- list()
# results_ajive$joint_tp <- sum(diag(jointProj %*% estimJointProj))
# results_ajive$joint_fd <- sum(diag((diag(n)-jointProj) %*% estimJointProj))
# results_ajive$indiv1_tp <- sum(diag(indiv1Proj %*% estimIndiv1Proj))
# results_ajive$indiv1_fd <- sum(diag((diag(n)-indiv1Proj) %*% estimIndiv1Proj))
# results_ajive$indiv2_tp <- sum(diag(indiv2Proj %*% estimIndiv2Proj))
# results_ajive$indiv2_fd <- sum(diag((diag(n)-indiv2Proj) %*% estimIndiv2Proj))
# 
# # fd control
# bound.approx <- function(min.X1, min.X2, r1, r2, beta1, beta2){
#   c1.min <- (min.X1^4 - beta1) / (min.X1 ^ 4 + beta1 * min.X1^2)
#   c2.min <- (min.X2^4 - beta2) / (min.X2 ^ 4 + beta2 * min.X2^2)
#   
#   ksi <- sqrt(2) * (sqrt(1 - sqrt(c1.min)) + 
#                       sqrt(1 - sqrt(c2.min)) + 
#                       sqrt(2) * sqrt(1 - sqrt(c1.min)) * sqrt(1 - sqrt(c2.min)))
#   bound <- (1/2) * (r1 + r2) * (2*ksi + ksi^2)
#   
#   eta <- sqrt(1 - c1.min) + sqrt(1 - c2.min) + sqrt(1 - c1.min) * sqrt(1 - c2.min)
#   bound <- bound + min(c(r1, r2)) * eta
#   
#   return(bound)
# }
# bound.approx.angle <- function(min.X1, min.X2, r1, r2, beta1, beta2, phi.min, phi.max){
#   c1.min <- (min.X1^4 - beta1) / (min.X1 ^ 4 + beta1 * min.X1^2)
#   c2.min <- (min.X2^4 - beta2) / (min.X2 ^ 4 + beta2 * min.X2^2)
#   
#   ksi <- sqrt(2) * (sqrt(1 - sqrt(c1.min)) + 
#                       sqrt(1 - sqrt(c2.min)) + 
#                       sqrt(2) * sqrt(1 - sqrt(c1.min)) * sqrt(1 - sqrt(c2.min)))
#   bound <- (2*ksi + ksi^2)
#   
#   eta <- phi.max * sqrt(c1.min) * sqrt(c2.min) + 
#     sqrt(1 - phi.min^2) * (sqrt(1 - c1.min) + sqrt(1 - c2.min)) + 
#     sqrt(1 - c1.min) * sqrt(1 - c2.min)
#   bound <- bound + eta 
#   
#   return(bound)
# }
# min.X1 <- 30 / (median(svd(X1)$d) / sqrt(qmp(0.5, ncol(X1), nrow(X1))))
#   # min(svd(blocks$decomp[[1]]$joint$full + blocks$decomp[[1]]$individual$full)$d[svd(blocks$decomp[[1]]$joint$full + blocks$decomp[[1]]$individual$full)$d > 0.1])
# min.X2 <- 15 / (median(svd(X2)$d) / sqrt(qmp(0.5, ncol(X2), nrow(X2))))
#   # min(svd(blocks$decomp[[2]]$joint$full + blocks$decomp[[2]]$individual$full)$d[svd(blocks$decomp[[2]]$joint$full + blocks$decomp[[2]]$individual$full)$d > 0.1])
# X1 <- X1 / (median(svd(X1)$d) / sqrt(qmp(0.5, ncol(X1), nrow(X1))))
# X2 <- X2 / (median(svd(X2)$d) / sqrt(qmp(0.5, ncol(X2), nrow(X2))))
# bound.approx.val <- bound.approx.angle(min.X1, min.X2, args$rj+args$ri1, args$rj+args$ri2, 0.2, 0.2, 0, 0.6)
# args$boundJoint <- bound.approx.val
# out <- fd_control_joint.nosubsample(X1, X2, args)
# estimJointProj <- out$joint %*% t(out$joint) # compute projection matrix
# estimIndiv1Proj <- out$indiv1 %*% t(out$indiv1)
# estimIndiv2Proj <- out$indiv2 %*% t(out$indiv2)
# if (all(is.na(out$joint))){
#   estimJointProj <- matrix(0, nrow=n, ncol=n)
# }
# if (all(is.na(out$indiv1))){
#   estimIndiv1Proj <- matrix(0, nrow=n, ncol=n)
# }
# if (all(is.na(out$indiv2))){
#   estimIndiv2Proj <- matrix(0, nrow=n, ncol=n)
# }
# 
# results_fd <- list()
# results_fd$joint_tp <- sum(diag(jointProj %*% estimJointProj))
# results_fd$joint_fd <- sum(diag((diag(n)-jointProj) %*% estimJointProj))
# results_fd$indiv1_tp <- sum(diag(indiv1Proj %*% estimIndiv1Proj))
# results_fd$indiv1_fd <- sum(diag((diag(n)-indiv1Proj) %*% estimIndiv1Proj))
# results_fd$indiv2_tp <- sum(diag(indiv2Proj %*% estimIndiv2Proj))
# results_fd$indiv2_fd <- sum(diag((diag(n)-indiv2Proj) %*% estimIndiv2Proj))
# 
# # make table of AJIVE and FD results
# results <- rbind(results_ajive, results_fd)
# rownames(results) <- c("AJIVE", "FD")
library(mgcv)

x <- sample(seq(0, 10, 0.05))
y <- sin(x) + rnorm(length(x), mean = 0, sd = 0.2)
model <- gam(y ~ s(x, k=20))
plot(x, y)
lines(x, model$fitted.values[order(x)])
plot(x, y)
lines(x, model$y)
plot(x, y)
lines(x[order(x)], model$y[order(x)])
plot(x, y)
lines(x[order(x)], model$fitted.values[order(x)])
plot(x, y)
lines(seq(0, 10, 0.05), predict(model, newdata = data.frame(x = seq(0, 10, 0.05))))
