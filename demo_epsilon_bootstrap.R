my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(reticulate)
library(denoiseR)
library(MASS)
library(ggplot2)

library(foreach)
library(doParallel)
source('src/generate_data_2_views.R')
source('src/utils.R')
source('src/models_2_views.R')
source('src/metrics.R')

# define other params
rj <- 4
ri1 <- 5
ri2 <- 4
m <- 50
phi_max <- 90
n1 <- 80
n2 <- 100
snr1 <- 10
snr2 <- 20
signal_strength1 <- 50
signal_strength2 <- 50
sigma1 <- (signal_strength1 / snr1) / (sqrt(m) + sqrt(n1))
sigma2 <- (signal_strength2 / snr2) / (sqrt(m) + sqrt(n2))
rank_spec <- 0
no_joint <- FALSE
no_indiv <- FALSE


# bootstrap
bootstrap.epsilon <- function(Y1, Y2, num_iter = 100) {
  # estimate rank
  shrink.Y1 <- optishrink(Y1)
  shrink.Y2 <- optishrink(Y2)
  rank1 <- shrink.Y1$nb.eigen
  rank2 <- shrink.Y2$nb.eigen
  
  # compute residual
  svd.Y1 <- svd(Y1)
  X1.hat <- svd.Y1$u[, 1:rank1] %*% diag(svd.Y1$d[1:rank1]) %*% t(svd.Y1$v[, 1:rank1])
  svd.Y2 <- svd(Y2)
  X2.hat <- svd.Y2$u[, 1:rank2] %*% diag(svd.Y2$d[1:rank2]) %*% t(svd.Y2$v[, 1:rank2])
  
  E1.hat <- Y1 - X1.hat
  E2.hat <- Y2 - X2.hat
  
  # impute
  U1.noise <- shrink.Y1$low.rank$u[, 1:rank1]
  V1.noise <- shrink.Y1$low.rank$v[, 1:rank1]
  U2.noise <- shrink.Y2$low.rank$u[, 1:rank2]
  V2.noise <- shrink.Y2$low.rank$v[, 1:rank2]
  D1 <- rmp(rank1, svr = min(nrow(Y1), ncol(Y1)) / max(nrow(Y1), ncol(Y1)))
  D2 <- rmp(rank2, svr = min(nrow(Y2), ncol(Y2)) / max(nrow(Y2), ncol(Y2)))

  E1.hat <- E1.hat + U1.noise %*% diag(D1) %*% t(V1.noise)
  E2.hat <- E2.hat + U2.noise %*% diag(D2) %*% t(V2.noise)
  
  # resampling
  out <- c()
  for (i in 1:num_iter) {
    # resample signal
    U1 <- svd(matrix(rnorm(nrow(Y1) * rank1), nrow = nrow(Y1), ncol = rank1))$u[, 1:rank1]
    V1 <- svd(matrix(rnorm(ncol(Y1) * rank1), nrow = ncol(Y1), ncol = rank1))$u[, 1:rank1]
    U2 <- svd(matrix(rnorm(nrow(Y2) * rank2), nrow = nrow(Y2), ncol = rank2))$u[, 1:rank2]
    V2 <- svd(matrix(rnorm(ncol(Y2) * rank2), nrow = ncol(Y2), ncol = rank2))$u[, 1:rank2]
    X1 <- U1 %*% diag(shrink.Y1$low.rank$d[1:rank1]) %*% t(V1)
    X2 <- U2 %*% diag(shrink.Y2$low.rank$d[1:rank2]) %*% t(V2)
    
    # form signal + noise
    Y1.resample <- X1 + E1.hat
    Y2.resample <- X2 + E2.hat
    
    # estimate column space
    U1.hat <- svd(Y1.resample)$u[, 1:rank1]
    U2.hat <- svd(Y2.resample)$u[, 1:rank2]
    
    # estimate spectral error in column space estimation
    # epsilon1 <- svd(U1.hat %*% t(U1.hat) - U1 %*% t(U1))$d[1]
    # epsilon2 <- svd(U2.hat %*% t(U2.hat) - U2 %*% t(U2))$d[1]
    # epsilon <- (1 + epsilon1) * (1 + epsilon2) - 1
    # out <- c(out, epsilon)
    # print(epsilon1)
    # print(epsilon2)
    # print(epsilon)
    
    # Adjustments
    P1epsilon1 <- svd(U1 %*% t(U1) %*% U1.hat %*% t(U1.hat) - U1 %*% t(U1))$d[1]
    P2epsilon2 <- svd(U2 %*% t(U2) %*% U2.hat %*% t(U2.hat) - U2 %*% t(U2))$d[1]
    E1E2 <- svd((U1.hat %*% t(U1.hat) - U1 %*% t(U1)) %*% (U2.hat %*% t(U2.hat) - U2 %*% t(U2)))$d[1]
    epsilon <- P1epsilon1 + P2epsilon2 + E1E2
    out <- c(out, epsilon)
    print(epsilon)
  }
  return (out)
}

# out <- c()
# for (snr in c(1, 2, 4, 8, 16, 32)) {
#   res <- c()
#   for (i in 1:10) {
#     sigma1 <- (signal_strength1 / snr) / (sqrt(m) + sqrt(n1))
#     sigma2 <- (signal_strength2 / snr) / (sqrt(m) + sqrt(n2))
#     data <- generate_data(m, n1, n2,
#                           rj, ri1, ri2, rank_spec,
#                           signal_strength1, signal_strength2,
#                           sigma1, sigma2,
#                           no_joint, no_indiv,
#                           phi_max)
#     epsilon <- mean(bootstrap.epsilon(data$Y1, data$Y2, num_iter = 100))
#     res <- c(res, epsilon)
#   }  
#   out <- c(out, mean(res))
# }

load("data/COADdata.rda")
# clean response data: delete 0 valued data and create class idx matrix
# all available
ID1 = (!is.na(COAD$subtype)) & (!apply(COAD$X1, 1, function(o) sum(is.na(o)))) &
  (!apply(COAD$X2, 1, function(o) sum(is.na(o))))
# check the numbers of obeservations for each case
sum(ID1 == T) # 167 - no missing

# extract subjects that have no missing information
rnaData = COAD$X1[ID1, ]
mRNAData = COAD$X2[ID1, ]
data = list(Y1 = rnaData,
            Y2 = mRNAData)

# bootstrap
epsilon <- bootstrap.epsilon(data$Y1, data$Y2, num_iter = 100)