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
  signal_strength1 = rep(10, num_exp),
  signal_strength2 = rep(12, num_exp),
  sigma1 = c(1, 5, 10),
  sigma2 = c(1, 6, 10)
)

# create dataframe with columns singular values, product, and parameters
df <- data.frame(
  Singular_value = numeric(0),
  Product = character(0),
  Parameters = character(0)
)
param.title <- c()
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
  d1 <- svd(matrix(signal_strength1 * rnorm(m * n1), m, n1))$d[1:(rj+ri1)]
  d2 <- svd(matrix(signal_strength2 * rnorm(m * n2), m, n2))$d[1:(rj+ri2)]
  d1 <- sample(d1)
  d2 <- sample(d2)
  dj1 <- d1[1:rj]
  dj2 <- d2[1:rj]
  di1 <- d1[(rj+1):(rj+ri1)]
  di2 <- d2[(rj+1):(rj+ri2)]
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
  # angle
  angle <- min(acos(svd(t(Ui1) %*% Ui2)$d))

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
  df <- rbind(df, data.frame(
    Singular_value = c(svd(P %*% Q)$d, svd(P.hat %*% Q.hat)$d),
    Product = rep(c('PQ', 'P.hat Q.hat'), each = length(svd(P %*% Q)$d)),
    Parameters = rep(exp_id, 2 * length(svd(P %*% Q)$d))
  ))
  param.title <- c(param.title, paste(c('SNR = ', signal_strength1 / sigma1, 
                              ', ', 'Angle = ', round(angle, 2)), collapse = ''))
}

#facet grid of histograms
p <- ggplot(df, aes(x = Singular_value)) +
  geom_histogram(binwidth = 0.1, position = 'dodge') +
  facet_wrap(~factor(Product, c('PQ', 'P.hat Q.hat'), labels = c('True product', 'Estimated product')) + 
               factor(Parameters, labels = param.title), scales = 'free') +
  theme_minimal() +
  xlab('Singular values') +
  ylab('Count')
