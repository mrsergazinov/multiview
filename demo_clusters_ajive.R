my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(MASS)
library(tidyverse)
library(Ckmeans.1d.dp)
library(ajive)
library(RMTstat)
library(denoiseR)
library(pracma)
library(gridExtra)
source('src/bounds_ajive.R')
source('src/models_2_views.R')
source('src/generate_data_2_views.R')
source('src/bounds.R')

set.seed(1)
# create list with params for simulation
num_exp <- 3
rj <- 4
ri1 <- 4
ri2 <- 5
rank1 <- rj + ri1
rank2 <- rj + ri2
m <- 50
n1 <- 60
n2 <- 70
sim_iter <- 100
signal_strength1 <- 100
signal_strength2 <- 100
params <- list(
  angles = matrix(c(50/180*pi, 45/180*pi, 40/180*pi,
                    50/180*pi, 45/180*pi, 40/180*pi,
                    50/180*pi, 45/180*pi, 40/180*pi), 
                  nrow = num_exp, byrow = TRUE),
  sigma = c(1.32, 2.2, 4.5)
  # sigma = c(1.32,1.32,1.32)
  # c(0.1, 0.15, 0.2, 0.3, 0.4, 0.6)
)

# create dataframe with columns singular values, product, and parameters
plt_ppd <- list()
plt_ajive <- list()
offset = 0.02
for (exp_id in 1:num_exp) {
  data <- generate_data(m, n1, n2,
                        rj, ri1, ri2, 0,
                        signal_strength1, signal_strength2, 
                        params$sigma[exp_id], params$sigma[exp_id],
                        FALSE, FALSE,
                        angles=params$angles[exp_id,])
  U1 <- svd(data$X1)$u[, 1:rank1, drop = FALSE]
  U2 <- svd(data$X2)$u[, 1:rank2, drop = FALSE]
  # compute PPD
  out <- proposed_func(data$Y1, data$Y2, rank1, rank2, return_scores = TRUE)
  spectrum_ppd <- out$test$svd.prod$d
  spectrum_ppd.true <- svd(t(U1) %*% U2)$d
  noise_bound_ppd <- out$test$lam
  joint_bound_ppd <- mean(out$bootstrap$epsilon1)
  # compute ajive
  spectrum_ajive <- compute_spectrum_ajive(data$Y1, data$Y2, rank1, rank2)
  spectrum_ajive.true <- svd(cbind(U1, U2))$d^2/2
  out <- compute_bounds_ajive(data$Y1, data$Y2, rank1, rank2, n_wedin_samples=1000, n_rand_dir_samples=1000)
  noise_bound_ajive <- out$random_directions_bound
  joint_bound_ajive <- out$wedin_bound
  
  # make title
  snr1 <- signal_strength1 / (params$sigma[exp_id] * (sqrt(m) + sqrt(n1)))
  snr2 <- signal_strength2 / (params$sigma[exp_id] * (sqrt(m) + sqrt(n2)))
  snr <- mean(c(snr1, snr2))
  angle_print = params$angles[exp_id, 3] / pi * 180
  title <- paste0("SNR = ", round(snr, 1))
  
  # plot histogram for ppd spectrum
  plt_ppd[[exp_id]] <- ggplot(data.frame(sing.vals = spectrum_ppd), aes(x = sing.vals)) +
    geom_histogram(closed = "right") +
    ggtitle(paste0("PPD, " , title)) +
    xlab('Singular values') +
    ylab('Frequency') + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  plt_ppd[[exp_id]] <- plt_ppd[[exp_id]] + geom_vline(xintercept = spectrum_ppd.true, color = 'red')
  max_height <- max(ggplot_build(plt_ppd[[exp_id]])$data[[1]]$y)
  plt_ppd[[exp_id]] <- plt_ppd[[exp_id]] + annotate("rect",
                                                    xmin = 0, xmax = noise_bound_ppd,
                                                    ymin = 0, ymax = max_height,
                                                    alpha = 0.3, fill = 'blue')
  plt_ppd[[exp_id]] <- plt_ppd[[exp_id]] + annotate("rect",
                                                    xmin = 1-joint_bound_ppd, xmax = 1,
                                                    ymin = 0, ymax = max_height,
                                                    alpha = 0.3, fill = 'green')
  
  
  # plot histogram for ajive spectrum
  plt_ajive[[exp_id]] <- ggplot(data.frame(sing.vals = spectrum_ajive), aes(x = sing.vals)) +
    geom_histogram(closed = "right") +
    ggtitle(paste0("AJIVE, " , title)) +
    xlab('Singular values') +
    ylab('Frequency') + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  plt_ajive[[exp_id]] <- plt_ajive[[exp_id]] + geom_vline(xintercept = spectrum_ajive.true, color = 'red')
  max_height <- max(ggplot_build(plt_ajive[[exp_id]])$data[[1]]$y)
  max_length <- max(ggplot_build(plt_ajive[[exp_id]])$data[[1]]$x)
  plt_ajive[[exp_id]] <- plt_ajive[[exp_id]] + annotate("rect",
                                                        xmin = 0, xmax = noise_bound_ajive^2/2,
                                                        ymin = 0, ymax = max_height,
                                                        alpha = 0.3, fill = 'blue')
  plt_ajive[[exp_id]] <- plt_ajive[[exp_id]] + annotate("rect",
                                                        xmin = joint_bound_ajive^2/2, xmax = 1,
                                                        ymin = 0, ymax = max_height,
                                                        alpha = 0.3, fill = 'green')
}

# arrrange ppd and ajive plot: ajive on top row and ppd on bottom row
grid.arrange(grobs = list(plt_ajive[[1]], plt_ajive[[2]], plt_ajive[[3]], 
                          plt_ppd[[1]], plt_ppd[[2]], plt_ppd[[3]]),
             ncol = 3, nrow = 2)





