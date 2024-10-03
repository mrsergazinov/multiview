library(tidyverse)
library(ajive)
library(RMTstat)
library(pracma)
library(gridExtra)
source('src/generate_data_2_views.R')
source('src/bounds.R')

set.seed(1)
# create list with params for simulation
num_exp <- 6
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
  angles = matrix(c(NA, NA, NA,
                    80/180*pi, 75/180*pi, 65/180*pi,
                    70/180*pi, 65/180*pi, 60/180*pi,
                    60/180*pi, 55/180*pi, 50/180*pi,
                    50/180*pi, 45/180*pi, 40/180*pi,
                    40/180*pi, 35/180*pi, 30/180*pi), 
                  nrow = num_exp, byrow = TRUE),
  sigma = c(0.1, 0.15, 0.2, 0.3, 0.4, 0.6)
)

# create dataframe with columns singular values, product, and parameters
plts <- list()
offset = 0.02
for (exp_id in 1:num_exp) {
  data <- generate_data(m, n1, n2,
                        rj, ri1, ri2, 0,
                        signal_strength1, signal_strength2, 
                        params$sigma[exp_id], params$sigma[exp_id],
                        FALSE, FALSE,
                        angles=params$angles[exp_id,])
  # compute estimated P, Q
  U1.hat <- svd(data$Y1)$u[, 1:rank1, drop = FALSE]
  U2.hat <- svd(data$Y2)$u[, 1:rank2, drop = FALSE]
  P.hat <- U1.hat %*% t(U1.hat)
  Q.hat <- U2.hat %*% t(U2.hat)
  sing.vals <- svd(P.hat %*% Q.hat)$d
  # compute bound
  R.upper <- bound.upper(data$Y1, data$Y2, data$X1, data$X2, data$Z1, data$Z2, rank1, rank2)
  R.lower <- bound.lower(data$Y1, data$Y2, data$X1, data$X2, data$Z1, data$Z2, rank1, rank2)
  # compute true sing.vals
  xticks <- c(0, 1)
  if (! all(is.na(params$angles[exp_id, ]))) {
    xticks <- c(xticks, cos(params$angles[exp_id, ]))
  }
  
  # plot histogram of singular values
  snr1 <- signal_strength1 / (params$sigma[exp_id] * (sqrt(m) + sqrt(n1)))
  snr2 <- signal_strength2 / (params$sigma[exp_id] * (sqrt(m) + sqrt(n2)))
  snr <- mean(c(snr1, snr2))
  angle_print = params$angles[exp_id, 3] / pi * 180
  # check if angle is na
  if (is.na(angle_print)) {
    angle_print = 90
  }
  plts[[exp_id]] <- ggplot(data.frame(sing.vals = sing.vals), aes(x = sing.vals)) +
    geom_histogram(closed = "right") +
    ggtitle(paste0("Angle = ", round(angle_print, 0), 
                   ", ",
                   "SNR = ", round(snr, 0))) +
    xlab('Singular values') +
    ylab('Frequency') + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  # max bin height
  max_height <- max(ggplot_build(plts[[exp_id]])$data[[1]]$y)
  # add bounds via annotation
  plts[[exp_id]] <- plts[[exp_id]] + annotate("rect",
                                              xmin = -offset, xmax = R.upper+offset,
                                              ymin = 0, ymax = max_height,
                                              alpha = 0.3, fill = 'red')
  plts[[exp_id]] <- plts[[exp_id]] + annotate("rect",
                                              xmin = 1-R.lower-offset, xmax = 1+offset,
                                              ymin = 0, ymax = max_height,
                                              alpha=0.3, fill='green')
  if (!all(is.na(params$angles[exp_id, ]))) {
    indiv.sing.vals <- cos(params$angles[exp_id, ])
    plts[[exp_id]] <- plts[[exp_id]] + annotate("rect",
                                                xmin = max(0, min(indiv.sing.vals)-R.lower-offset),
                                                xmax = min(1, max(indiv.sing.vals)+R.upper+offset),
                                                ymin = 0, ymax = max_height,
                                                alpha=0.3, fill='blue')
  }
  plts[[exp_id]] <- plts[[exp_id]] + geom_vline(xintercept = xticks, color = 'red')
}

# facet grid of 6 plots
grid.arrange(grobs = plts, ncol = 3)

