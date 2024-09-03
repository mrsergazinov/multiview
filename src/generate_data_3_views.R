generate_data <- function(m, n1, n2, n3,
                          rj, ri1, ri2, ri3,
                          rank_spec, 
                          signal_strength1, signal_strength2, signal_strength3,
                          sigma1, sigma2, sigma3,
                          no_joint, no_indiv, 
                          phi_max=NA, angles=NA) {
  # data generator
  d1 <- runif(rj+ri1, min=signal_strength1, max=2*signal_strength1)
  d2 <- runif(rj+ri2, min=signal_strength2, max=2*signal_strength2)
  d3 <- runif(rj+ri3, min=signal_strength3, max=2*signal_strength3)
  d1 <- sample(d1)
  d2 <- sample(d2)
  d3 <- sample(d3)
  dj1 <- d1[1:rj]
  dj2 <- d2[1:rj]
  dj3 <- d3[1:rj]
  di1 <- d1[(rj+1):(rj+ri1)]
  di2 <- d2[(rj+1):(rj+ri2)]
  di3 <- d3[(rj+1):(rj+ri3)]
  # generate data
  U <- svd(matrix(rnorm(m * m), m, m))$u
  # joint parts
  Uj <- U[, 1:rj]
  O1 <- randortho(rj, type = 'orthonormal') # rotate Uj
  O2 <- randortho(rj, type = 'orthonormal') # rotate Uj
  Uj1 <- Uj
  Uj2 <- Uj %*% O1
  Uj3 <- Uj %*% O2
  # individual part
  Ui1 <- U[, (rj+1):(rj+ri1)]
  Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
  Ui3 <- U[, (rj+ri1+ri2+1):(rj+ri1+ri2+ri3)]
  # rotate 
  if (phi_max != 90) {
    phi_max = phi_max * pi / 180
    ra1 <- as.integer(min(ri1, ri2)  / 2)
    angles1 <- c(phi_max, phi_max + runif(ra1-1, min=0, max=pi/2-phi_max))
    Ui2[, 1:ra1] <- Ui1[, 1:ra1] %*% diag(cos(angles1), nrow=ra1) + Ui2[, 1:ra1] %*% diag(sin(angles1), nrow=ra1) 
    ra2 <- as.integer(min(ri1, ri3)  / 2)
    angles2 <- c(phi_max, phi_max + runif(ra2-1, min=0, max=pi/2-phi_max))
    Ui3[, 1:ra2] <- Ui1[, 1:ra2] %*% diag(cos(angles2), nrow=ra2) + Ui3[, 1:ra2] %*% diag(sin(angles2), nrow=ra2)
  } else if (!all(is.na(angles))) {
    # raise not implemented error
    stop("Angles for 3 views through `agnles` argument are not implemented")
  }
  # no-joint or no-indiv cases
  if (no_joint) {
    Uj1 <- matrix(0, m, rj)
    Uj2 <- matrix(0, m, rj)
    Uj3 <- matrix(0, m, rj)
  }
  if (no_indiv) {
    Ui1 <- matrix(0, m, ri1)
    Ui2 <- matrix(0, m, ri2)
    Ui3 <- matrix(0, m, ri3)
  }
  # loadings
  Vj1 <- svd(matrix(rnorm(m * n1), m, n1))$v[, 1:rj]
  Vj2 <- svd(matrix(rnorm(m * n2), m, n2))$v[, 1:rj]
  Vj3 <- svd(matrix(rnorm(m * n3), m, n3))$v[, 1:rj]
  Vi1 <- svd(matrix(rnorm(m * n1), m, n1))$v[, 1:ri1]
  Vi2 <- svd(matrix(rnorm(m * n2), m, n2))$v[, 1:ri2]
  Vi3 <- svd(matrix(rnorm(m * n3), m, n3))$v[, 1:ri3]
  # combine
  X1 <- Uj1 %*% diag(dj1) %*% t(Vj1) + Ui1 %*% diag(di1) %*% t(Vi1)
  Z1 <- matrix(rnorm(m * n1), m, n1) * sigma1
  Y1 <- X1 + Z1
  X2 <- Uj2 %*% diag(dj2) %*% t(Vj2) + Ui2 %*% diag(di2) %*% t(Vi2)
  Z2 <- matrix(rnorm(m * n2), m, n2) * sigma2
  Y2 <- X2 + Z2
  X3 <- Uj3 %*% diag(dj3) %*% t(Vj3) + Ui3 %*% diag(di3) %*% t(Vi3)
  Z3 <- matrix(rnorm(m * n3), m, n3) * sigma3
  Y3 <- X3 + Z3
  
  # compute true
  Pjoint <- Uj1 %*% t(Uj1) # true projection
  Pindiv1 <- Ui1 %*% t(Ui1)
  Pindiv2 <- Ui2 %*% t(Ui2)
  Pindiv3 <- Ui3 %*% t(Ui3)
  P1 <- Pjoint + Pindiv1
  P2 <- Pjoint + Pindiv2
  P3 <- Pjoint + Pindiv3
  
  # signal rank
  error1 <- rank_spec * sample(0:2, 1)
  error2 <- rank_spec * sample(0:2, 1)
  error3 <- rank_spec * sample(0:2, 1)
  if (no_joint) {
    rank1 <- ri1 + error1
    rank2 <- ri2 + error2
    rank3 <- ri3 + error3
  } else if (no_indiv) {
    rank1 <- rj + error1
    rank2 <- rj + error2
    rank3 <- rj + error3
  } else {
    rank1 <- rj + ri1 + error1
    rank2 <- rj + ri2 + error2
    rank3 <- rj + ri3 + error3
  }
  
  return(list(Y1=Y1, Y2=Y2, Y3=Y3,
              X1=X1, X2=X2, X3=X3,
              Z1=Z1, Z2=Z2, Z3=Z3,
              P1=P1, P2=P2, P3=P3,
              Pjoint=Pjoint, Pindiv1=Pindiv1, Pindiv2=Pindiv2, Pindiv3=Pindiv3,
              rank1=rank1, rank2=rank2, rank3=rank3,
              rj=rj, ri1=ri1, ri2=ri2, ri3=ri3))
}