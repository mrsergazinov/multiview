my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(denoiseR)
library(MASS)
library(pracma)
library(RMTstat)
library(Ckmeans.1d.dp)

library(gridExtra)
library(ggplot2)
library(tidyr)

library(foreach)
library(doParallel)
source('src/generate_data_2_views.R')
source('src/utils.R')
source('src/models_2_views.R')

# define other params
rj <- 4
ri1 <- 5
ri2 <- 4
m <- 90
n1 <- 110
n2 <- 120
signal_strength1 <- 50
signal_strength2 <- 50
rank_spec <- 0
no_joint <- FALSE
no_indiv <- FALSE

snrs <- seq(0.5, 2, 0.1)
phi_maxs <- c(90, 60, 30)
nsim <- 10
table <- matrix(0, nrow = length(snrs) * length(phi_maxs) * nsim, ncol=9)
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
      # bootstrap epsilon2
      out <- proposed_func(data$Y1, data$Y2, data$rank1, data$rank2, return_scores = TRUE, bootstrap_iters = 200)
      epsilon2.bootstrap <- mean(out$bootstrap$epsilon2)
      # true epsilon2
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
      epsilon2.true <- svd(P1E1P2 + P1E2P2 + P1E1E2P2)$d[1]
      # lambda
      lambda <- out$test$lam
      # estimate variance
      sigma1.hat <- est.sigma(data$Y1)
      sigma2.hat <- est.sigma(data$Y2)
      # store results
      table[(j-1)*nsim*length(phi_maxs) + (k-1)*nsim + i, ] <- c(phi_max, snr, 
                                                                 sigma1.hat, sigma2.hat, sigma1, sigma2,
                                                                 epsilon2.bootstrap, epsilon2.true, lambda)
    }
  }
}
df <- data.frame(table)
colnames(df) <- c("Angle", "SNR",
                  "Sigma1 hat", "Sigma2 hat", "Sigma1", "Sigma2",
                  "Bootstrap epsilon", "True epsilon", "Lambda")


df <- df[df$Angle == 60, c("SNR", "Bootstrap epsilon", "True epsilon", "Lambda")]
df <- df %>% pivot_longer(
  cols = c("Bootstrap epsilon", "True epsilon", "Lambda"),
  names_to = "variable", 
  values_to = "value"
)
ggplot(df, aes(x = SNR, y = value, color = variable)) +
  stat_summary(fun = mean, geom = "line") +  # Plot the mean as a line
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", aes(fill = variable), alpha = 0.2) +  # Add confidence intervals
  labs(title = "Curves with Confidence Intervals", x = "SNR", y = "Y") +
  theme_minimal()

# plts <- list()
# plts.sigma <- list()
# for (phi_max in phi_maxs) {
#   df1 <- df[df$Angle == phi_max, ]
#   plt <- ggplot(data = df1, aes(x = SNR, y = `Bootstrap epsilon minus epsilon`)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE) +
#     theme_minimal() +
#     xlab("SNR") +
#     ylab(TeX("$\\hat{\\epsilon} - \\epsilon$")) + 
#     ggtitle(paste("Angle = ", phi_max))
#   plts[[length(plts) + 1]] <- plt
#   
#   df2 <- df.sigma[df.sigma$Angle == phi_max, ]
#   plt <- ggplot(data = df2, aes(x = SNR, y = `Sigma error`)) +
#     geom_point() +
#     geom_smooth(method = "lm", se = FALSE) +
#     theme_minimal() +
#     xlab("SNR") +
#     ylab(TeX("$\\hat{s} - s$"))
#   plts.sigma[[length(plts.sigma) + 1]] <- plt
# }
# grid.arrange(grobs = c(plts, plts.sigma), ncol = 3)