my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(pracma)

generate_data <- function(m, n1, n2, 
                          rj, ri1, ri2, rank_spec, 
                          signal_strength1, signal_strength2, 
                          sigma1, sigma2,
                          no_joint, no_indiv, 
                          phi_max) {
  # data generator
  d1 <- runif(rj+ri1, min=signal_strength1, max=2*signal_strength1)
  d2 <- runif(rj+ri2, min=signal_strength2, max=2*signal_strength2)
  d1 <- sample(d1)
  d2 <- sample(d2)
  dj1 <- d1[1:rj]
  dj2 <- d2[1:rj]
  di1 <- d1[(rj+1):(rj+ri1)]
  di2 <- d2[(rj+1):(rj+ri2)]
  # generate data
  U <- svd(matrix(rnorm(m * m), m, m))$u
  # joint parts
  Uj <- U[, 1:rj]
  O <- randortho(rj, type = 'orthonormal') # rotate Uj
  Uj1 <- Uj %*% O
  Uj2 <- Uj
  # individual part
  Ui1 <- U[, (rj+1):(rj+ri1)]
  Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
  # rotate 
  if (phi_max != 90) {
    phi_max = phi_max * pi / 180
    ra <- as.integer(min(ri1, ri2)  / 2)
    angles <- c(phi_max, phi_max + runif(ra-1, min=0,max=pi/2-phi_max))
    Ui2[, 1:ra] <- Ui1[, 1:ra] %*% diag(cos(angles), nrow=ra) + Ui2[, 1:ra] %*% diag(sin(angles), nrow=ra) 
  }
  # no-joint or no-indiv cases
  if (no_joint) {
    Uj1 <- matrix(0, m, rj)
    Uj2 <- matrix(0, m, rj)
  }
  if (no_indiv) {
    Ui1 <- matrix(0, m, ri1)
    Ui2 <- matrix(0, m, ri2)
  }
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
  
  # compute true
  Pjoint <- Uj1 %*% t(Uj1) # true projection
  Pindiv1 <- Ui1 %*% t(Ui1)
  Pindiv2 <- Ui2 %*% t(Ui2)
  P1 <- Pjoint + Pindiv1
  P2 <- Pjoint + Pindiv2
  
  # signal rank
  error1 <- rank_spec * sample(1:2, 1)
  error2 <- rank_spec * sample(1:2, 1)
  if (no_joint) {
    rank1 <- ri1 + error1
    rank2 <- ri2 + error2
  } else if (no_indiv) {
    rank1 <- rj + error1
    rank2 <- rj + error2
  } else {
    rank1 <- rj + ri1 + error1
    rank2 <- rj + ri2 + error2
  }
  
  return(list(Y1=Y1, Y2=Y2, 
              P1=P1, P2=P2, 
              Pjoint=Pjoint, Pindiv1=Pindiv1, Pindiv2=Pindiv2,
              rank1=rank1, rank2=rank2,
              rj=rj, ri1=ri1, ri2=ri2))
}

