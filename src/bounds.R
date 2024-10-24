f_lambda <- function(lambda, q1, q2) {
  # density of the product of two free projections
  # Compute lambda_+ and lambda_-
  lambda_plus <- q1 + q2 - 2 * q1 * q2 + 2 * sqrt(q1 * q2 * (1 - q1) * (1 - q2))
  lambda_minus <- q1 + q2 - 2 * q1 * q2 - 2 * sqrt(q1 * q2 * (1 - q1) * (1 - q2))
  
  # Calculate f(lambda)
  numerator <- sqrt((lambda_plus - lambda) * (lambda - lambda_minus))
  denominator <- 2 * pi * lambda * (1 - lambda)
  
  # Ensure we handle cases where denominator might be zero
  result <- ifelse(denominator != 0, numerator / denominator, NA)
  
  return(result)
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
bound.lower <- function(Y1, Y2, X1, X2, Z1, Z2, rank1, rank2, rank1.hat, rank2.hat) {
  svd.X1 <- svd(X1)
  svd.X2 <- svd(X2)
  svd.Y1 <- svd(Y1)
  svd.Y2 <- svd(Y2)
  
  P1 <- svd.X1$u[, 1:rank1] %*% t(svd.X1$u[, 1:rank1])
  P2 <- svd.X2$u[, 1:rank2] %*% t(svd.X2$u[, 1:rank2])
  P1.hat <- svd.Y1$u[, 1:rank1.hat] %*% t(svd.Y1$u[, 1:rank1.hat])
  P2.hat <- svd.Y2$u[, 1:rank2.hat] %*% t(svd.Y2$u[, 1:rank2.hat])
  
  E1 <- P1.hat - P1
  E2 <- P2.hat - P2
  
  return(svd(P1 %*% E2 %*% P2 + P1 %*% E1 %*% P2 + P1 %*% E1 %*% E2 %*% P2)$d[1])
}
bound.upper <- function(Y1, Y2, X1, X2, Z1, Z2, rank1, rank2, rank1.hat, rank2.hat) {
  svd.X1 <- svd(X1)
  svd.X2 <- svd(X2)
  svd.Y1 <- svd(Y1)
  svd.Y2 <- svd(Y2)
  
  P1 <- svd.X1$u[, 1:rank1] %*% t(svd.X1$u[, 1:rank1])
  P2 <- svd.X2$u[, 1:rank2] %*% t(svd.X2$u[, 1:rank2])
  P1.hat <- svd.Y1$u[, 1:rank1.hat] %*% t(svd.Y1$u[, 1:rank1.hat])
  P2.hat <- svd.Y2$u[, 1:rank2.hat] %*% t(svd.Y2$u[, 1:rank2.hat])
  
  E1 <- P1.hat - P1
  E2 <- P2.hat - P2
  
  return(svd(P1 %*% E2 + E1 %*% P2 + E1 %*% E2)$d[1])
}