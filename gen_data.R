gen_data <- function(n, p1, p2, rj, ri1, ri2, dj, di1, di2, sigma1, sigma2) {
  # Generate two views with dim n x p1 and n x p2 respectively
  # Arguments:
  # n: sample size
  # p1: dimension of view 1
  # p2: dimension of view 2
  # rj: joint rank
  # ri1: individual rank of view 1
  # ri2: individual rank of view 2
  # dj: joint eigenvalues
  # di1: individual eigenvalues of view 1
  # di2: individual eigenvalues of view 2
  # sigma1: noise level of view 1
  # sigma2: noise level of view 2
  # Return:
  # X1: view 1
  # X2: view 2
  # joint: joint part

  U <- svd(matrix(rnorm(n * (p1+p2)), n, p1+p2))$u
  # joint part
  Uj <- U[, 1:rj]
  V1 <- svd(matrix(rnorm(rj * p1), rj, p1))$v
  V2 <- svd(matrix(rnorm(rj * p2), rj, p2))$v
  X1.joint <- Uj %*% diag(dj) %*% t(V1)
  X2.joint <- Uj %*% diag(dj) %*% t(V2)
  # individual part
  X1.indiv <- matrix(0, n, p1)
  if (ri1 > 0){
    Ui1 <- U[, (rj+1):(rj+ri1)]
    V1 <- svd(matrix(rnorm(ri1 * p1), ri1, p1))$v
    X1.indiv <- Ui1 %*% diag(di1) %*% t(V1)
  }
  X2.indiv <- matrix(0, n, p2)
  if (ri2 > 0){
    Ui2 <- U[, (rj+ri1+1):(rj+ri1+ri2)]
    V2 <- svd(matrix(rnorm(ri2 * p2), ri2, p2))$v
    X2.indiv <- Ui2 %*% diag(di2) %*% t(V2)
  }
  # combine with noise
  X1 <- X1.joint + X1.indiv + matrix(rnorm(n * p1), n, p1) * sigma1
  X2 <- X2.joint + X2.indiv + matrix(rnorm(n * p2), n, p2) * sigma2
  
  return(list("X1" = X1, "X2" = X2, 
              "X1.joint" = X1.joint, "X1.indiv" = X1.indiv,
              "X2.joint" = X2.joint, "X1.indiv" = X2.indiv,
              "joint" = Uj, "jointPerp" = U[, (rj+1):ncol(U)], "indiv1" = Ui1, "indiv2" = Ui2))
}