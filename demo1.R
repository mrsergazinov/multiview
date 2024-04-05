library(tidyverse)
library(ajive)
library(RMTstat)
library(pracma)
library(latex2exp)

bound <- function(X1.d, X2.d, indiv.cos, beta1, beta2) {
  c1 <- (X1.d^4 - beta1) / (X1.d^4 + beta1*X1.d^2)
  c2 <- (X2.d^4 - beta2) / (X2.d^4 + beta2*X2.d^2)

  R1 <- 2 * sqrt(2 * sum(1 - sqrt(c1))) + 2 * sum(1 - sqrt(c1))
  R2 <- 2 * sqrt(2 * sum(1 - sqrt(c2))) + 2 * sum(1 - sqrt(c2))

  return (indiv.cos + R1 + R2 + R1 * R2)
}
bound.v <- function(d1, d2, sigma1, sigma2, m, n1, n2) {
  projection_bound <- function(d, sigma, m, n) {
    normE <- (sqrt(m) + sqrt(n)) * sigma
    normE.quad <- 2 * max(d) * normE + normE^2
    s <- min(d)^2
    kappa <- 2 * normE.quad
    return (normE.quad * kappa / ((s - kappa) * (s - 3 * kappa / 2)))
  }
  E1 <- projection_bound(d1, sigma1, m, n1)
  E2 <- projection_bound(d2, sigma2, m, n2)
  return (E1 + E2 + E1 * E2)
}
bound.cai <- function(Y1, Y2, X1, X2, Z1, Z2, rank1, rank2) {
  projection_bound <- function(Y, X, Z, rank) {
    svd.X <- svd(X)
    U <- svd.X$u[, 1:rank]
    U.orth <- svd.X$u[, (rank+1):ncol(svd.X$u)]
    V <- svd.X$v[, 1:rank]
    V.orth <- svd.X$v[, (rank+1):ncol(svd.X$v)]
    alpha = min(svd(t(U) %*% Y %*% V)$d)
    beta = max(svd(t(U.orth) %*% Y %*% V.orth)$d)
    Pu <- U %*% t(U)
    Pv <- V %*% t(V)
    Pu.orth <- U.orth %*% t(U.orth)
    Pv.orth <- V.orth %*% t(V.orth)
    Z21 <- Pu.orth %*% Z %*% Pv
    Z12 <- Pu %*% Z %*% Pv.orth
    z21 <- max(svd(Z21)$d)
    z12 <- max(svd(Z12)$d)
    
    min((alpha * z21 + beta * z12) / 
          (alpha^2 - beta^2 - min(z21^2, z12^2)), 
        1)
  }
  E1 <- projection_bound(Y1, X1, Z1, rank1)
  E2 <- projection_bound(Y2, X2, Z2, rank2)
  return (E1 + E2 + E1 * E2)
}

# create list with params for simulation
num_exp <- 6
params <- list(
  rj = rep(4, num_exp),
  ri1 = rep(3, num_exp),
  ri2 = rep(3, num_exp),
  m = rep(30, num_exp),
  phi_max = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
  n1 = rep(30, num_exp),
  n2 = rep(40, num_exp),
  sim_iter = rep(100, num_exp),
  signal_strength1 = rep(10, num_exp),
  signal_strength2 = rep(12, num_exp),
  sigma1 = c(0.1, 0.15, 0.2, 0.3, 0.4, 0.6),
  sigma2 = c(0.1, 0.15, 0.2, 0.3, 0.4, 0.6)
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
  Z1 <- matrix(rnorm(m * n1), m, n1) * sigma1
  Y1 <- X1 + Z1
  X2 <- Uj2 %*% diag(dj2) %*% t(Vj2) + Ui2 %*% diag(di2) %*% t(Vi2)
  Z2 <- matrix(rnorm(m * n2), m, n2) * sigma2
  Y2 <- X2 + Z2
  # angle
  d <- svd(t(Ui1) %*% Ui2)$d
  angle <- min(acos(d))
  angle.min <- min(d)
  angle.max <- max(d)
  
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

  # compute bounds
  R <- bound.cai(Y1, Y2, X1, X2, Z1, Z2, rank1, rank2)

  # use ggplot to plot singular values of P %*% Q and P.hat %*% Q.hat
  df <- rbind(df, data.frame(
    Singular_value = svd(P.hat %*% Q.hat)$d,
    Parameters = exp_id,
    Bound1 = R,
    Bound2 = max(angle.min - R, 0),
    Bound3 = max(angle.max + R, 0),
    Bound4 = 1 - R
  ))
  param.title <- c(param.title, paste(c('SNR = ', round(signal_strength1 / sigma1, 2), 
                                        ', ', 'Angle = ', round(angle, 2)), collapse = ''))

}

#facet grid of histograms
p <- ggplot(df, aes(x = Singular_value)) +
  geom_rect(aes(xmin = -0.05, xmax = Bound1 + 0.05, ymin = 0, ymax = Inf), fill = 'red', alpha = 0.01) +
  geom_rect(aes(xmin = Bound2 - 0.05, xmax = Bound3 + 0.05, ymin = 0, ymax = Inf), fill = 'blue', alpha = 0.01) +
  geom_rect(aes(xmin = Bound4 - 0.05, xmax = 1.05, ymin = 0, ymax = Inf), fill = 'green', alpha = 0.01) +
  geom_histogram(binwidth = 0.1, position = 'dodge') +
  facet_wrap(~factor(Parameters, labels = param.title), scales = 'free', ncol=3) +
  theme_minimal() +
  xlab('Singular values') +
  ylab('Count')
