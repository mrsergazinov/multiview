library(tidyverse)
library(ajive)
source("gen_data.R")
file_fd <- "results_fd.RData"
file_metric <- "results_metric.RData"
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
sigma2 <- 1
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

jointFdControl <- function(X1, X2, numSamples=100, alpha=0.7){
  avg.P <- matrix(0, nrow=nrow(X1), ncol=nrow(X1))
  thresh.X1 <- thresh(X1) # thresholding singular values
  thresh.X2 <- thresh(X2)
  for (i in 1:numSamples){
    X1.sample <- X1[, sample(1:ncol(X1), as.integer(ncol(X1)/2), replace=FALSE)]
    X2.sample <- X2[, sample(1:ncol(X2), as.integer(ncol(X2)/2), replace=FALSE)]
    svd.X1.sample <- svd(X1.sample)
    svd.X2.sample <- svd(X2.sample)
    u1.sample <- svd.X1.sample$u[, svd.X1.sample$d > thresh.X1] # thresholding using Gavish and Donoho 2014
    u2.sample <- svd.X2.sample$u[, svd.X2.sample$d > thresh.X2]
    avg.step <- (u1.sample %*% t(u1.sample) + u2.sample %*% t(u2.sample)) / 2
    avg.P <- avg.P + avg.step
  }
  avg.P <- avg.P / numSamples
  svd.avg <- svd(avg.P)
  out <- svd.avg$u[, svd.avg$d > alpha]
  return(out)
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
  jointFdControlOut <- jointFdControl(X1, X2, alpha = alpha)
  colProjJointFdControl <- jointFdControlOut %*% t(jointFdControlOut)
  # compute true col space
  trueProj <- Uj %*% t(Uj)
  
  # compute FD and metric
  fds[i] <- sum(diag((diag(n)-trueProj) %*% colProjJointFdControl))
  metric[i] <- sum((trueProj - colProjJointFdControl)^2)
}
# add a row to results.fd
df <- data.frame(
  SNR = snr,  # Creating an empty numeric vector for SNR
  Method = "fdControl joint",
  Mean = mean(fds),
  SD = sd(fds)
)
results.fd <- rbind(results.fd, df)

df <- data.frame(
  SNR = snr,  # Creating an empty numeric vector for SNR
  Method = "fdControl joint",
  Mean = mean(metric),
  SD = sd(metric)
)
results.metric <- rbind(results.metric, df)

# save results.fd and results.metric
save(results.fd, file = file_fd)
save(results.metric, file = file_metric)