library(tidyverse)
library(ajive)
source("gen_data.R")
file_fd <- "results_diffNoise_fd.RData"
file_metric <- "results_diffNoise_metric.RData"
load(file_fd)
load(file_metric)

# set simulation parameters
set.seed(235017)
alpha <- 0.9
rj <- 4
ri1 <- 10
ri2 <- 20
n <- 100
p1 <- 100
p2 <- 150
sigma1 <- 1
sigma2 <- 10
sim_iter <- 100
signal_strength <- 15
dj <- rnorm(rj, mean = signal_strength, sd = 2)
di1 <- rnorm(ri1, mean = signal_strength, sd = 2)
# rep(0, ri1)
di2 <- rnorm(ri2, mean = signal_strength, sd = 2)
# rep(0, ri2)
snr1 <- (sum(dj ** 2) + sum(di1 ** 2)) / sum(n * p1)
snr2 <- (sum(dj ** 2) + sum(di2 ** 2)) / sum(n * p2)
snr <- (snr1 + snr2) / 2

thresh <- function(X1, sigma = 1){
  # thresholding using Gavish and Donoho 2017
  return(1 + sqrt(min(ncol(X1), nrow(X1)) / max(ncol(X1), nrow(X1)))) 
}

fdControl <- function(X, numSamples=100, alpha=0.7){
    # X is n x p matrix
    n <- nrow(X)
    p <- ncol(X)
    avgP <- matrix(0, nrow=n, ncol=n)
    for (i in 1:numSamples){
        # resample cols
        Xsample <- X[, sample(1:p, as.integer(p/2), replace=FALSE)]
        thresh.X = thresh(Xsample) # thresholding singular valuess
        svd.Xsample = svd(Xsample)
        U <- svd.Xsample$u[, svd.Xsample$d > thresh.X] # decompose and extract left singular vectors
        P <- U %*% t(U) # projection matrix
        avgP <- avgP + P
    }
    avgP <- avgP / numSamples
    svdAvg <- svd(avgP) 
    return (svdAvg$u[, svdAvg$d > alpha]) # select left singular vectors for which singular values are > alpha
}

# create empty table with colnames "SNR", "Method", "Mean", "SD"
fds <- rep(0, sim_iter)
metric <- rep(0, sim_iter)
for (i in 1:sim_iter){
  data <- gen_data(n, p1, p2, rj, ri1, ri2, dj, di1, di2, sigma1, sigma2)
  X1 <- data[["X1"]]
  X2 <- data[["X2"]]
  Uj <- data[["joint"]]
  
  # apply oracle ajive
  fdControlOut <- fdControl(X2, alpha = alpha)
  colProjfdControl <- fdControlOut %*% t(fdControlOut)
  # compute true col space
  trueProj <- Uj %*% t(Uj)
  
  # compute FD and metric
  fds[i] <- sum(diag((diag(n)-trueProj) %*% colProjfdControl))
  metric[i] <- sum((trueProj - colProjfdControl)^2)
}
# add a row to results.fd
df <- data.frame(
  SNR = snr,  # Creating an empty numeric vector for SNR
  Method = "fdControl X2",
  Mean = mean(fds),
  SD = sd(fds)
)
results.fd <- rbind(results.fd, df)

df <- data.frame(
  SNR = snr,  # Creating an empty numeric vector for SNR
  Method = "fdControl X2",
  Mean = mean(metric),
  SD = sd(metric)
)
results.metric <- rbind(results.metric, df)

# save results.fd and results.metric
save(results.fd, file = file_fd)
save(results.metric, file = file_metric)