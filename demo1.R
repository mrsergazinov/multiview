library(tidyverse)
library(ajive)
library(RMTstat)
library(pracma)
library(latex2exp)

# set.seed(235017)
rj <- 4
ri1 <- 3
ri2 <- 3
m <- 30
phi_max <- 0.1
n1 <- 100
n2 <- 70
sim_iter <- 100
signal_strength1 <- 5
signal_strength2 <- 5
sigma1 = 1
sigma2 = 1
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

bound <- function(X1.d, X2.d, indiv.cos, beta1, beta2){
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
# compute bound
bound.val <- bound.cai(Y1, Y2, X1, X2, Z1, Z2, rank1, rank2)

# use ggplot to plot singular values of P %*% Q and P.hat %*% Q.hat
df <- data.frame(
  singular_value = c(svd(P %*% Q)$d, svd(P.hat %*% Q.hat)$d),
  Product = rep(c('PQ', 'P.hat Q.hat'), each = length(svd(P %*% Q)$d))
)
plt <- ggplot(df, aes(x = singular_value, fill = Product)) +
  geom_histogram(binwidth = 0.05, position = 'dodge') +
  theme_minimal() +
  xlab('Singular values') +
  ylab('Count') +
  scale_fill_discrete(labels = c(expression(hat(P)[ColY1]*hat(P)[ColY2]), expression(P[ColX1]*P[ColX2]))) + 
  labs(title = paste('Singular values of PQ and PQ.hat',
                     '\nn = ', m, ', n1 = ', n1, ', n2 = ', n2,
                     '\nphi_max = ', phi_max, ', mean(X1.d) = ', round(mean(c(dj1, di1)), 2),
                     ', mean(X2.d) = ', round(mean(c(dj2, di2)), 2),
                     '\nsigma1 = ', round(sigma1, 2), ', sigma2 = ', round(sigma2, 2))) +
  geom_vline(xintercept = bound.val, linetype = 'dashed', color = 'red')
print(plt)
