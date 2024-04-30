check_and_compute <- function(X, Y = NULL) {
  # If both X and Y are NULL, return infinity and a zero matrix
  if (is.null(X) && is.null(Y)) {
    return(list("r" = Inf, "P" = matrix(0, nrow = dims, ncol = dims)))
  }
  # Calculate for single or combined inputs
  compute_P <- function(Z) {
    if (is.null(Z)) matrix(0, nrow = dims, ncol = dims) else Z %*% t(Z)
  }
  P_X <- compute_P(X)
  P_Y <- compute_P(Y)
  return(list("r" = ifelse(is.null(X), 0, ncol(X)) + ifelse(is.null(Y), 0, ncol(Y)),
              "P" = P_X + P_Y))
}

global_null_2_views <- function(Y1, Y2, rank1, rank2, compute_prod = TRUE) {
  m <- nrow(Y1)
  q1 <- rank1 / m
  q2 <- rank2 / m
  q.plus <- q1 + q2 - 2*q1*q2 + 2*sqrt(q1*q2*(1-q1)*(1-q2))
  q.minus <- q1 + q2 - 2*q1*q2 - 2*sqrt(q1*q2*(1-q1)*(1-q2))
  if (q.minus == 0) {
    q.minus <- 1e-6 # offset to avoid comp error
  }
  A0 <- 1 - min(q1, q2) 
  pdf <- function(x) {
    sqrt((q.plus - x)*(x - q.minus)) /(2 * pi * x * (1-x))
  }
  cdf <- function(lam) {
    A0 + integrate(pdf, q.minus, lam)$value - 0.95
  }
  # find closest to 0.95
  lam <- 0
  if (A0 > 0.95) {
    lam <- 0
  } else {
    lam <- uniroot(cdf, c(q.minus, q.plus))$root
  }
  # compute product of projections
  if (compute_prod) {
    U1.hat <- svd(Y1)$u[, 1:rank1, drop = FALSE]
    U2.hat <- svd(Y2)$u[, 1:rank2, drop = FALSE]
    P.hat <- U1.hat %*% t(U1.hat)
    Q.hat <- U2.hat %*% t(U2.hat)
    prod <- P.hat %*% Q.hat
    prod.sym <- (prod + t(prod)) / 2
    svd.prod <- svd(prod)
    return (list("noJoint" = (svd.prod$d[1] < lam),
                 "P.hat" = P.hat,
                 "Q.hat" = Q.hat,
                 "prod.sym" = prod.sym,
                 "svd.prod" = svd.prod))
  } else {
    return (lam)
  }
  
}