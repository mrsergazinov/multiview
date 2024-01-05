library(tidyverse)
library(ajive)
library(RMTstat)

source("gen_data.R")
source("fd_control_joint.R")
source("ajive_oracle.R")
source("ajive.R")
source("bounds.R")

file_name <- "results.RData"
load(file_name)

# models to test
model_list = list("fd_control_joint" = fd_control_joint,
                  "ajive" = ajive_wrapper,
                  "ajive_oracle" = ajive_oracle_wrapper)
model_name <- "fd_control_joint"
model <- model_list[[model_name]]

# set simulation parameters
set.seed(235017)
rj <- 2
ri1 <- 3
ri2 <- 2
n <- 20
p1 <- 100
p2 <- 100
sigma1 <- 1
sigma2 <- 1
sim_iter <- 100
signal_strength <- 40
dj <- rnorm(rj, mean = signal_strength, sd = 3)
di1 <- rnorm(ri1, mean = signal_strength, sd = 3)
di2 <- rnorm(ri2, mean = signal_strength, sd = 3)
# compute SNR -- 2 ways
snr1hat1 <- (sum(dj ^ 2) + sum(di1 ^ 2)) / sum(n * p1 * sigma1 ^ 2)
snr2hat1 <- (sum(dj ^ 2) + sum(di2 ^ 2)) / sum(n * p2 * sigma2 ^ 2)
maxSigma1 <- sigma1 * sqrt(2 * log(6) * (n + p1))
maxSigma2 <- sigma2 * sqrt(2 * log(6) * (n + p2))
snr1hat2 <- min(c(di1, dj)) / maxSigma1
snr2hat2 <- min(c(di2, dj)) / maxSigma2
# compute spectral bound
bound.val <- bound(c(dj, di1), c(dj, di2), rj, ri1, ri2, n / p1)
bound.approx.val <- bound.approx(c(dj, di1), c(dj, di2), rj, ri1, ri2, n / p1)
bound.approx.angle.val <- bound.approx.angle(c(dj, di1), c(dj, di2), rj, ri1, ri2, n / p1, 0.01, 0.05)
print(paste("Spectral bound = ", bound.val))
print(paste("Approx spectral bound = ", bound.approx.val))
print(paste("Approx angle spectral bound = ", bound.approx.angle.val))

# set args for models
args = list("sigma1" = NA, "sigma2" = NA,
            "rj" = rj, "ri1" = ri1,
            "ri2" = ri2, "numSamples" = 100,
            "alpha" = 0.4, "boundJoint" = bound.val)

# for (alpha in seq(0.7, 0.9, 0.05)) {
# args$alpha = alpha

sing.vals <- c()
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
    out <- model(X1, X2, args)
    sing.vals <- c(out$avgPd, sing.vals1)

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
    Method = model_name,
    MeanJointFD = mean(fds_joint),
    MeanJointTP = mean(tps_joint),
    MeanIndiv1FD = mean(fds_indiv1),
    MeanIndiv1TP = mean(tps_indiv1),
    MeanIndiv2FD = mean(fds_indiv2),
    MeanIndiv2TP = mean(tps_indiv2),
    SNR_X1 = snr1hat1,
    SNR_X2 = snr2hat1,
    SNR2_X1 = snr1hat2,
    SNR2_X2 = snr2hat2,
    RankJoint = rj,
    RankIndiv1 = ri1,
    RankIndiv2 = ri2,
    sigma1 = sigma1,
    sigma2 = sigma2,
    nrow = n,
    ncol1 = p1,
    ncol2 = p2,
    alpha = args$alpha,
    stabilityIter = args$numSamples,
    paramSigma1 = args$sigma1,
    paramSigma2 = args$sigma2
)
results <- rbind(results, df)
# }

hist_plot <- ggplot(data = data.frame(x = sing.vals), aes(x = x)) +
  geom_histogram(binwidth = 0.1, fill = "orange", color = "black", alpha = 0.7) +
  ggtitle(paste("Histogram of eigenvalues\nBound = ", bound.val, 
                "\nMean signal eigenvalue = ", signal_strength,
                "\nBeta = ", p1/n),
          ) +
  xlab("Values") +
  ylab("Frequency") +
  theme_minimal()

# save results.fd and results.metric
save(results, file = file_name)


