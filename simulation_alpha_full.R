library(tidyverse)
library(ajive)
library(RMTstat)

source("gen_data.R")
source("fd_control_joint.R")

file_name <- "results_alpha_full.RData"
load(file_name)

# set simulation parameters
set.seed(235017)
rj <- 4
ri1 <- 10
ri2 <- 20
n <- 100
p1 <- 120
p2 <- 150
sigma1 <- 0.3
sigma2 <- 0.3
sim_iter <- 10
signal_strength <- 15
dj <- rnorm(rj, mean = signal_strength, sd = 3)
di1 <- rnorm(ri1, mean = signal_strength, sd = 3)
di2 <- rnorm(ri2, mean = signal_strength, sd = 3)
snr1hat1 <- (sum(dj ^ 2) + sum(di1 ^ 2)) / sum(n * p1 * sigma1 ^ 2)
snr2hat1 <- (sum(dj ^ 2) + sum(di2 ^ 2)) / sum(n * p2 * sigma2 ^ 2)
maxSigma1 <- sigma1 * sqrt(2*log(6)*(n+p1))
maxSigma2 <- sigma2 * sqrt(2*log(6)*(n+p2))
snr1hat2 <- min(c(di1, dj)) / maxSigma1
snr2hat2 <- min(c(di2, dj)) / maxSigma2

for (alpha in seq(0.7, 0.95, 0.01)) {
    args = list("sigma1" = NA, "sigma2" = NA,
            "rj" = rj, "ri1" = ri1, "ri2" = ri2,
            "sigma" = NA, "numSamples" = 100,
            "alpha" = alpha)

    fds_joint <- rep(0, sim_iter)
    tps_joint <- rep(0, sim_iter)
    fds_indiv1 <- rep(0, sim_iter)
    tps_indiv1 <- rep(0, sim_iter)
    fds_indiv2 <- rep(0, sim_iter)
    tps_indiv2 <- rep(0, sim_iter)
    for (i in 1:sim_iter){
        data <- gen_data(n, p1, p2, rj, ri1, ri2, dj, di1, di2, sigma1, sigma2)
        X1 <- data[["X1"]]
        X2 <- data[["X2"]]
        Uj <- data[["joint"]]
        Ui1 <- data[["indiv1"]]
        Ui2 <- data[["indiv2"]]
        jointProj <- Uj %*% t(Uj) # true projection
        indiv1Proj <- Ui1 %*% t(Ui1)
        indiv2Proj <- Ui2 %*% t(Ui2)

        # apply models from the list
        out <- fd_control_joint(X1, X2, args)

        # compute projection matrix
        estimJointProj <- out$joint %*% t(out$joint)
        estimIndiv1Proj <- out$indiv1 %*% t(out$indiv1)
        estimIndiv2Proj <- out$indiv2 %*% t(out$indiv2)
        if (all(is.na(out$joint))){
            print(paste(c('joint is null space ->', model_name), collapse = " "))
            estimJointProj <- matrix(0, nrow=n, ncol=n)
        }
        if (all(is.na(out$indiv1))){
            print(paste(c('indiv1 is null space ->', model_name), collapse = " "))
            estimIndiv1Proj <- matrix(0, nrow=n, ncol=n)
        }
        if (all(is.na(out$indiv2))){
            print(paste(c('indiv2 is null space ->', model_name), collapse = " "))
            estimIndiv2Proj <- matrix(0, nrow=n, ncol=n)
        }

        fds_joint[i] <- sum(diag((diag(n)-jointProj) %*% estimJointProj))
        tps_joint[i] <- sum(diag(jointProj %*% estimJointProj))
        fds_indiv1[i] <- sum(diag((diag(n)-indiv1Proj) %*% estimIndiv1Proj))
        tps_indiv1[i] <- sum(diag(indiv1Proj %*% estimIndiv1Proj))
        fds_indiv2[i] <- sum(diag((diag(n)-indiv2Proj) %*% estimIndiv2Proj))
        tps_indiv2[i] <- sum(diag(indiv2Proj %*% estimIndiv2Proj))
    }
    # add a row to results.fd
    df <- data.frame(
      Method = paste(c("alpha = ", alpha), collapse = " "),
      MeanJointFD = mean(fds_joint),
      MeanJointTP = mean(tps_joint),
      MeanIndiv1FD = mean(fds_indiv1),
      MeanIndiv1TP = mean(tps_indiv1),
      MeanIndiv2FD = mean(fds_indiv2),
      MeanIndiv2TP = mean(tps_indiv2),
      Description = paste(c('snr1hat1 =', round(snr1hat1, 2), 'snr2hat1 =', round(snr2hat1, 2),
                            'snr1hat2 =', round(snr1hat2, 2), 'snr2hat2 =', round(snr2hat2, 2),
                            'rj =', rj, 'ri1 =', ri1, 'ri2 =', ri2, 
                            'n = ', n, 'p1 = ', p1, 'p2 = ', p2), collapse = " "))
    results <- rbind(results, df)
}

# save results.fd and results.metric
save(results, file = file_name)

    
    