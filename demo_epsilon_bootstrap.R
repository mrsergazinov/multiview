my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(denoiseR)
library(MASS)
library(pracma)
library(RMTstat)
library(Ckmeans.1d.dp)

library(ggpubr)
library(ggplot2)
library(tidyr)

library(foreach)
library(doParallel)
source('src/generate_data_2_views.R')
source('src/utils.R')
source('src/models_2_views.R')

# define other params
rj <- 2
ri1 <- 2
ri2 <- 2
m <- 50
n1 <- 80
n2 <- 100
signal_strength1 <- 50
signal_strength2 <- 50
rank_spec <- 0
no_joint <- FALSE
no_indiv <- FALSE

snrs <- c(0.75, 0.8, 0.9, 1, 1.5, 2, 3)
phi_maxs <- c(60)
nsim <- 10
table <- matrix(0, nrow = length(snrs) * length(phi_maxs) * nsim, ncol=10)
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
      # bootstrap epsilon1
      out <- proposed_func(data$Y1, data$Y2, data$rank1, data$rank2, 
                           return_scores = TRUE, bootstrap_iters = 100, rotation_correction = TRUE)
      epsilon1.bootstrap.corrected <- mean(out$bootstrap$epsilon1)
      # bootstrap epsilon1
      out <- proposed_func(data$Y1, data$Y2, data$rank1, data$rank2, 
                           return_scores = TRUE, bootstrap_iters = 100, rotation_correction = FALSE)
      epsilon1.bootstrap <- mean(out$bootstrap$epsilon1)
      # true epsilon1
      if (is.null(out$joint)) {
        P1.hat <- out$indiv1 %*% t(out$indiv1)
        P2.hat <- out$indiv2 %*% t(out$indiv2)
      } else {
        P1.hat <- out$joint %*% t(out$joint) + out$indiv1 %*% t(out$indiv1)
        P2.hat <- out$joint %*% t(out$joint) + out$indiv2 %*% t(out$indiv2)  
      }
      E1 <- P1.hat - data$P1
      E2 <- P2.hat - data$P2
      P1E1P2 <- data$P1 %*% E1 %*% data$P2
      P1E2P2 <- data$P1 %*% E2 %*% data$P2
      P1E1E2P2 <- data$P1 %*% E1 %*% E2 %*% data$P2
      epsilon1.true <- svd(P1E1P2 + P1E2P2 + P1E1E2P2)$d[1]
      # lambda
      lambda <- out$test$lam
      # estimate variance
      sigma1.hat <- est.sigma(data$Y1)
      sigma2.hat <- est.sigma(data$Y2)
      # store results
      table[(j-1)*nsim*length(phi_maxs) + (k-1)*nsim + i, ] <- c(phi_max, snr, 
                                                                 sigma1.hat, sigma2.hat, sigma1, sigma2,
                                                                 epsilon1.true, epsilon1.bootstrap, 
                                                                 epsilon1.bootstrap.corrected, lambda)
    }
  }
}
df <- data.frame(table)
colnames(df) <- c("Angle", "SNR",
                  "Sigma1 hat", "Sigma2 hat", "Sigma1", "Sigma2",
                  "True epsilon", "Bootstrap epsilon (No Correction)", 
                  "Bootstrap epsilon", "Lambda")


df <- df[df$Angle == phi_maxs[1], c("SNR", "True epsilon", 
                                    "Bootstrap epsilon (No Correction)", "Bootstrap epsilon")]
df <- df %>% pivot_longer(
  cols = c("True epsilon", "Bootstrap epsilon (No Correction)", "Bootstrap epsilon"),
  names_to = "Variable", 
  values_to = "value"
)
plt <- ggplot(df, aes(x = SNR, y = value, color = Variable)) +
  stat_summary(fun = mean, geom = "line") +  # Plot the mean as a line
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", aes(fill = Variable), alpha = 0.2) +  # Add confidence intervals
  labs(title = paste0("Rank ratio = ", round((rj+ri1)/m*10, 2)*10, "%"), x = "SNR", y = "Epsilon") +
  theme_minimal()


# add back legend to all plots
ggarrange(plt4, plt1, plt5, plt3, plt2, plt6, ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")
