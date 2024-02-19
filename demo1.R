library(tidyverse)
library(ajive)
library(RMTstat)
library(pracma)

set.seed(235017)
rj <- 2
ri1 <- 2
ri2 <- 3
n <- 20
phi.max <- 0
p1 <- 100
p2 <- 120
sim_iter <- 100
signal_strength1 <- 40
signal_strength2 <- 30
sigma1 = 3.5
sigma2 = 2.5
dj1 <- rnorm(rj, mean = signal_strength1, sd = 3)
dj2 <- rnorm(rj, mean = signal_strength2, sd = 3)
di1 <- rnorm(ri1, mean = signal_strength1, sd = 3)
di2 <- rnorm(ri2, mean = signal_strength2, sd = 3)
# generate data
U <- svd(matrix(rnorm(n * (p1+p2)), n, p1+p2))$u
# joint part
Uj <- U[, 1:rj]
Ujperp <- U[, (rj+1):ncol(U)]
O <- randortho(rj, type = 'orthonormal') # rotate Uj
Uj1 <- Uj %*% O
Uj2 <- Uj
# individual part
Ui1 <- U[, (rj+1):(rj+ri1)]
Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
O <- matrix(runif(ri1 * ri2, -phi.max, phi.max), ri1, ri2) # rotate
Ui2 <- Ui2 + Ui1 %*% O
Ui2 <- gramSchmidt(Ui2)$Q # orthonormalize
# loadings
V1 <- svd(matrix(rnorm((rj + ri1) * p1), rj + ri1, p1))$v
V2 <- svd(matrix(rnorm((rj + ri2) * p2), rj + ri2, p2))$v
# combine
Y1 <- cbind(Uj1, Ui1) %*% diag(c(dj1, di1)) %*% t(V1)
Y2 <- cbind(Uj2, Ui2) %*% diag(c(dj2, di2)) %*% t(V2)
# noise
X1 <- Y1 + matrix(rnorm(n * p1), n, p1) * sigma1
X2 <- Y2 + matrix(rnorm(n * p2), n, p2) * sigma2


bound.angle <- function(X1.d, X2.d, rj, ri1, ri2, beta1, beta2, psi.min, psi.max){
  c1 <- (X1.d^4 - beta1) / (X1.d^4 + beta1*X1.d^2)
  c2 <- (X2.d^4 - beta2) / (X2.d^4 + beta2*X2.d^2)
  delta1 <- sqrt(1 - c1)
  delta2 <- sqrt(1 - c2)
  sqrt.delta1 <- sqrt(1 - sqrt(c1))
  sqrt.delta2 <- sqrt(1 - sqrt(c2))
  bound <- 0
  diags <- list()
  for (i in 1:(rj + ri1)) {
    for (j in 1:(rj + ri2)) {
      if (i == j & i < rj){
        ksi_1 <- sqrt(2) * (min(sqrt.delta1[i]) + 
                            sqrt.delta2[i] + 
                            sqrt(2) * min(sqrt.delta1[i]) * sqrt.delta2[i])
        ksi_2 <- sqrt(2) * (sqrt.delta1[i] + 
                            min(sqrt.delta2[i]) + 
                            sqrt(2) * sqrt.delta1[i] * min(sqrt.delta2[i]))
        ksi <- min(ksi_1, ksi_2)
        bound <- bound + 2 * ksi + ksi^2
      } else if (i > rj & j > rj){
        eta <- psi.max * sqrt(c1[i]) * sqrt(c2[j]) + 
          sqrt(1 - psi.min^2) * delta1[i] + delta2[j] + 
          (1 + psi.max) * delta1[i] * delta2[j]
        diags[[paste(i - j)]] <- c(diags[[paste(i - j)]], eta)
      } else {
        eta <- delta1[i] + delta2[j] + delta1[i] * delta2[j]
        diags[[paste(i - j)]] <- c(diags[[paste(i - j)]], eta)
      }
    }
  }
  for (name in names(diags)) {
    bound <- bound + max(diags[[name]])
  }
  return(bound)
}

# compute P, Q - true
P <- cbind(Uj1, Ui1) %*% t(cbind(Uj1, Ui1))
Q <- cbind(Uj2, Ui2) %*% t(cbind(Uj2, Ui2))
# compute estimated P, Q
thresh <- function(X, sigma = NA) {
  if (is.na(sigma)){
    sigma = median(svd(X)$d) / sqrt(qmp(0.5, ncol(X), nrow(X)))
  }
  beta = nrow(X) / ncol(X)
  lambda = 1 + sqrt(beta)
  return (sigma * lambda)
}
X1.tresh <- thresh(X1)
X2.tresh <- thresh(X2)
U1.hat <- svd(X1)$u[, svd(X1)$d > X1.tresh, drop = FALSE]
U2.hat <- svd(X2)$u[, svd(X2)$d > X2.tresh, drop = FALSE]
P.hat <- U1.hat %*% t(U1.hat)
Q.hat <- U2.hat %*% t(U2.hat)
# compute bound
bound.true <- bound.angle(c(dj, di1), c(dj, di2), rj, ri1, ri2, n / p1, n / p2, 0, phi.max)

# use ggplot to plot singular values of P %*% Q and P.hat %*% Q.hat
df <- data.frame(
  singular_value = c(svd(P %*% Q)$d, svd(P.hat %*% Q.hat)$d),
  matrix = rep(c('PQ', 'PQ.hat'), each = length(svd(P %*% Q)$d))
)
plt <- ggplot(df, aes(x = singular_value, fill = matrix)) +
  geom_histogram(binwidth = 0.05, position = 'dodge') +
  theme_minimal() +
  labs(title = paste('Singular values of PQ and PQ.hat', 
                     '\nn = ', n, ', p1 = ', p1, ', p2 = ', p2,
                     '\nphi.max = ', phi.max, ', mean(X1.d) = ', round(mean(c(dj1, di1)), 2),
                     ', mean(X2.d) = ', round(mean(c(dj2, di2)), 2), 
                     '\nsigma1 = ', round(sigma1, 2), ', sigma2 = ', round(sigma2, 2))) +
  geom_vline(xintercept = bound.true, linetype = 'dashed', color = 'red')
print(plt)



