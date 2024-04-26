bound.v <- function(d1, d2, sigma1, sigma2, m, n1, n2) {
  # bound using 
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
  # bound using Cai and Zhang
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