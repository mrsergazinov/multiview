my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(latex2exp)
library(reticulate)
library(denoiseR)
library(MASS)
library(ggplot2)
library(pracma)
library(RMTstat)
library(Ckmeans.1d.dp)
library(gridExtra)

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

snrs <- c(1, 2, 4, 8, 16, 32)
phi_maxs <- c(90, 60, 30)
nsim <- 10
table <- matrix(0, nrow = length(snrs) * length(phi_maxs) * nsim, ncol=3)
table.sigma <- matrix(0, nrow = length(snrs) * length(phi_maxs) * nsim, ncol=3)
for (j in 1:length(snrs)) {
  snr <- snrs[j]
  for (k in 1:length(phi_maxs)){
    phi_max <- phi_maxs[k]
    for (i in 1:nsim) {
      sigma1 <- (signal_strength1 / snr) / (sqrt(m) + sqrt(n1))
      sigma2 <- (signal_strength2 / snr) / (sqrt(m) + sqrt(n2))
      data <- generate_data(m, n1, n2,
                            rj, ri1, ri2, rank_spec,
                            signal_strength1, signal_strength2,
                            sigma1, sigma2,
                            no_joint, no_indiv,
                            phi_max / 180 * pi)
      sigma1.hat <- est.sigma(data$Y1)
      sigma2.hat <- est.sigma(data$Y2)
      out <- proposed_func(data$Y1, data$Y2, data$rank1, data$rank2, return_scores = TRUE, bootstrap_iters = 200)
      
      P1.hat <- out$joint %*% t(out$joint) + out$indiv1 %*% t(out$indiv1)
      P2.hat <- out$indiv2 %*% t(out$indiv2) + out$joint %*% t(out$joint)
      E1 <- P1.hat - data$P1
      E2 <- P2.hat - data$P2
      P1E1P2 <- data$P1 %*% E1 %*% data$P2
      P1E2P2 <- data$P1 %*% E2 %*% data$P2
      P1E1E2P2 <- data$P1 %*% E1 %*% E2 %*% data$P2
      epsilon <- svd(P1E1P2 + P1E2P2 + P1E1E2P2)$d[1]
      
      table.sigma[(j-1)*nsim*length(phi_maxs) + (k-1)*nsim + i, ] <- c(phi_max, snr, mean(c(sigma1.hat - sigma1, sigma2.hat - sigma2)))
      table[(j-1)*nsim*length(phi_maxs) + (k-1)*nsim + i, ] <- c(phi_max, snr, mean(out$epsilon) - epsilon)
    }
  }
}
df <- data.frame(table)
df.sigma <- data.frame(table.sigma)
colnames(df) <- c("Angle", "SNR", "Bootstrap epsilon minus epsilon")
colnames(df.sigma) <- c("Angle", "SNR", "Sigma error")
plts <- list()
plts.sigma <- list()
for (phi_max in phi_maxs) {
  df1 <- df[df$Angle == phi_max, ]
  plt <- ggplot(data = df1, aes(x = SNR, y = `Bootstrap epsilon minus epsilon`)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_minimal() +
    xlab("SNR") +
    ylab(TeX("$\\hat{\\epsilon} - \\epsilon$")) + 
    ggtitle(paste("Angle = ", phi_max))
  plts[[length(plts) + 1]] <- plt
  
  df2 <- df.sigma[df.sigma$Angle == phi_max, ]
  plt <- ggplot(data = df2, aes(x = SNR, y = `Sigma error`)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_minimal() +
    xlab("SNR") +
    ylab(TeX("$\\hat{s} - s$"))
  plts.sigma[[length(plts.sigma) + 1]] <- plt
}
grid.arrange(grobs = c(plts, plts.sigma), ncol = 3)


# load("data/COADdata.rda")
# # clean response data: delete 0 valued data and create class idx matrix
# # all available
# ID1 = (!is.na(COAD$subtype)) & (!apply(COAD$X1, 1, function(o) sum(is.na(o)))) &
#   (!apply(COAD$X2, 1, function(o) sum(is.na(o))))
# # check the numbers of obeservations for each case
# sum(ID1 == T) # 167 - no missing
# 
# # extract subjects that have no missing information
# rnaData = COAD$X1[ID1, ]
# mRNAData = COAD$X2[ID1, ]
# data = list(Y1 = rnaData,
#             Y2 = mRNAData)
# 
# rank1 <- 20
# rank2 <- 16
# out <- proposed_func(data$Y1, data$Y2, rank1, rank2, return_scores = TRUE)
# d <- out$test$svd.prod$d[1:min(rank1, rank2)]
# epsilon <- bootstrap.epsilon(data$Y1, data$Y2, rank1, rank2, d, num_iter = 100)