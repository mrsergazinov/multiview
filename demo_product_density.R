my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(ggplot2)
library(gridExtra)
source('src/bounds.R')


generate_data_plot <- function(q1, q2) {
  sing.vals <- c()
  for (i in 1:100) {
    rank1 <- as.integer(100 * q1)
    rank2 <- as.integer(100 * q2)
    U1 <- svd(matrix(rnorm(10000), nrow = 100))$u[, 1:rank1]
    U2 <- svd(matrix(rnorm(10000), nrow = 100))$u[, 1:rank2]
    sing.vals <- c(sing.vals, svd(t(U1) %*% U2)$d[1:max(rank1, rank2)]^2)
  }
  # plot histogram of singular values and overlay with the analytical density
  plt <- ggplot(data.frame(sing.vals = sing.vals), aes(x = sing.vals)) + 
    geom_histogram(aes(y = ..density..), bins = 30, color = 'black', fill = 'white') + 
    stat_function(fun = function(x) f_lambda(x, q1, q2), color = 'red') + 
    theme_minimal() +
    xlab('Singular values') + 
    ylab('Density') +
    ggtitle(paste0("Ratios q1 = q2 = ", q1)) +
    xlim(0.01, 1)  
  
  return( plt)
}

plt1 <- generate_data_plot(0.1, 0.1)
plt2 <- generate_data_plot(0.2, 0.2)
plt3 <- generate_data_plot(0.3, 0.3)
plt4 <- generate_data_plot(0.4, 0.4)
plt5 <- generate_data_plot(0.5, 0.5)
plt6 <- generate_data_plot(0.6, 0.6)

grid.arrange(plt1, plt2, plt3, plt4, plt5, plt6, ncol = 3)
