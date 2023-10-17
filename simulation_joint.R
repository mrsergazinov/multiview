library(tidyverse)
library(ajive)
library(RMTstat)

source("gen_data.R")
source("fd_control_joint.R")
source("ajive_oracle.R")
source("ajive.R")
source("fd_control.R")

file_name <- "results.RData"
load(file_name)

# models to test
model_names <- c("fd_control_joint")
# ajive, "ajive_oracle", "fd_control_joint", "fd_control"
models <- c(fd_control_joint)
# ajive_wrapper, ajive_oracle_wrapper, fd_control_joint, fd_control

# set simulation parameters
set.seed(235017)
rj <- 4
ri1 <- 8
ri2 <- 8
n <- 100
p1 <- 120
p2 <- 150
sigma1 <- 0.3
sigma2 <- 0.3
sim_iter <- 100
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

# set args for models
args = list("sigma1" = NA, "sigma2" = NA,
            "rj" = rj, "ri1" = ri1, "ri2" = ri2,
            "sigma" = NA, "numSamples" = 100,
            "alpha" = 0.9)


for (i in 1:length(model_names)){
    model_name <- model_names[i]
    model <- models[[i]]

    fds <- rep(0, sim_iter)
    tps <- rep(0, sim_iter)
    tds <- rep(0, sim_iter)
    for (i in 1:sim_iter){
        data <- gen_data(n, p1, p2, rj, ri1, ri2, dj, di1, di2, sigma1, sigma2)
        X1 <- data[["X1"]]
        X2 <- data[["X2"]]
        Uj <- data[["joint"]]
        trueProj <- Uj %*% t(Uj) # true projection matrix

        # apply models from the list
        out <- model(X1, X2, args)

        # extract joint
        if (model_name == 'fd_control_joint' || 
            model_name == 'ajive' || 
            model_name == 'ajive_oracle'){
            out <- out$joint
        }

        # compute projection matrix
        td = ncol(out)
        colProj <- out %*% t(out)
        if (all(is.na(out))){
            print(paste(c('joint is null space ->', model_name), collapse = " "))
            colProj <- matrix(0, nrow=n, ncol=n)
            td = 0
        }

        # compute FD and TP
        fds[i] <- sum(diag((diag(n)-trueProj) %*% colProj))
        tps[i] <- sum(diag(trueProj %*% colProj))
        tds[i] <- td
    }
    # add a row to results.fd
    df <- data.frame(
      Method = model_name,
      MeanFD = mean(fds),
      SdFD = sd(fds),
      MeanTP = mean(tps),
      SdTP = sd(tps),
      MeanTD = mean(tds),
      SdTD = sd(tds),
      Description = paste(c('snr1hat1 =', round(snr1hat1, 2), 'snr2hat1 =', round(snr2hat1, 2),
                            'snr1hat2 =', round(snr1hat2, 2), 'snr2hat2 =', round(snr2hat2, 2),
                            'rj =', rj, 'ri1 =', ri1, 'ri2 =', ri2, 
                            'n = ', n, 'p1 = ', p1, 'p2 = ', p2), collapse = " "))
    results <- rbind(results, df)
}

# save results.fd and results.metric
save(results, file = file_name)