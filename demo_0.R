library(tidyverse)
library(ajive)
library(RMTstat)
library(pracma)
library(latex2exp)
library(lemon)

# create list with params for simulation
num_exp <- 3
params <- list(
  rj = rep(4, num_exp),
  ri1 = rep(3, num_exp),
  ri2 = rep(3, num_exp),
  m = rep(30, num_exp),
  phi_max = c(0, 0.3, 0.5),
  n1 = rep(30, num_exp),
  n2 = rep(40, num_exp),
  sim_iter = rep(100, num_exp),
  signal_strength1 = rep(30, num_exp),
  signal_strength2 = rep(40, num_exp),
  sigma1 = c(1 / sqrt(30), 1, 3),
  sigma2 = c(1 / sqrt(40), 1, 4),
  dj1 = list(rnorm(4, mean = 30, sd = 3), rnorm(4, mean = 30, sd = 3), rnorm(4, mean = 30, sd = 3)),
  dj2 = list(rnorm(4, mean = 40, sd = 3), rnorm(4, mean = 40, sd = 3), rnorm(4, mean = 40, sd = 3)),
  di1 = list(rnorm(3, mean = 30, sd = 3), rnorm(3, mean = 30, sd = 3), rnorm(3, mean = 30, sd = 3)),
  di2 = list(rnorm(3, mean = 40, sd = 3), rnorm(3, mean = 40, sd = 3), rnorm(3, mean = 40, sd = 3))
)
plts <- list()
for (exp_id in 1:num_exp) {
  rj <- params$rj[exp_id]
  ri1 <- params$ri1[exp_id]
  ri2 <- params$ri2[exp_id]
  m <- params$m[exp_id]
  phi_max <- params$phi_max[exp_id]
  n1 <- params$n1[exp_id]
  n2 <- params$n2[exp_id]
  sim_iter <- params$sim_iter[exp_id]
  signal_strength1 <- params$signal_strength1[exp_id]
  signal_strength2 <- params$signal_strength2[exp_id]
  sigma1 = params$sigma1[exp_id]
  sigma2 = params$sigma2[exp_id]
  dj1 <- params$dj1[[exp_id]]
  dj2 <- params$dj2[[exp_id]]
  di1 <- params$di1[[exp_id]]
  di2 <- params$di2[[exp_id]]
  # generate data
  U <- svd(matrix(rnorm(m * (n1+n2)), m, n1+n2))$u
  # joint part
  Uj <- U[, 1:rj]
  Ujperp <- U[, (rj+1):ncol(U)]
  O <- randortho(rj, type = 'orthonormal') # rotate Uj
  Uj1 <- Uj %*% O
  Uj2 <- Uj
  # individual part
  Ui1 <- U[, (rj+1):(rj+ri1)]
  Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
  O <-  matrix(runif(ri1 * ri2, -phi_max, phi_max), ri1, ri2) # rotate
  Ui2 <- Ui2 + Ui1 %*% O
  Ui2 <- gramSchmidt(Ui2)$Q # orthonormalize
  # loadings
  Vj1 <- svd(matrix(rnorm(m * n1), m, n1))$v[, 1:rj]
  Vj2 <- svd(matrix(rnorm(m * n2), m, n2))$v[, 1:rj]
  Vi1 <- svd(matrix(rnorm(m * n1), m, n1))$v[, 1:ri1]
  Vi2 <- svd(matrix(rnorm(m * n2), m, n2))$v[, 1:ri2]
  # combine
  X1 <- Uj1 %*% diag(dj1) %*% t(Vj1) + Ui1 %*% diag(di1) %*% t(Vi1)
  Y1 <- X1 + matrix(rnorm(m * n1), m, n1) * sigma1
  X2 <- Uj2 %*% diag(dj2) %*% t(Vj2) + Ui2 %*% diag(di2) %*% t(Vi2)
  Y2 <- X2 + matrix(rnorm(m * n2), m, n2) * sigma2

  # compute P, Q - true
  P <- cbind(Uj1, Ui1) %*% t(cbind(Uj1, Ui1))
  Q <- cbind(Uj2, Ui2) %*% t(cbind(Uj2, Ui2))
  # compute estimated P, Q
  rank1 <- rj + ri1
  rank2 <- rj + ri2
  U1.hat <- svd(Y1)$u[, 1:rank1, drop = FALSE]
  U2.hat <- svd(Y2)$u[, 1:rank2, drop = FALSE]
  P.hat <- U1.hat %*% t(U1.hat)
  Q.hat <- U2.hat %*% t(U2.hat)

  # use ggplot to plot singular values of P %*% Q and P.hat %*% Q.hat
  df <- data.frame(
    singular_value = c(svd(P %*% Q)$d, svd(P.hat %*% Q.hat)$d),
    Product = rep(c('PQ', 'P.hat Q.hat'), each = length(svd(P %*% Q)$d))
  )
  plts[[exp_id]] <- ggplot(df, aes(x = singular_value, fill = Product)) +
    geom_histogram(binwidth = 0.05, position = 'dodge') +
    theme_minimal() +
    xlab('Singular values') +
    ylab('Count') +
    scale_fill_discrete(labels = c(expression(hat(P)[ColY1]*hat(P)[ColY2]), expression(P[ColX1]*P[ColX2]))) +
    theme(text = element_text(size = 20)) +
    ylim(c(0, 27))
}
grid_arrange_shared_legend(plts[[1]], plts[[2]], plts[[3]], ncol = 3, position='right')