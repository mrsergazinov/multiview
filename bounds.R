bound <- function(X1.d, X2.d, rj, ri1, ri2, beta){
  c1 <- (X1.d^4 - beta) / (X1.d^4 + beta*X1.d^2)
  c2 <- (X2.d^4 - beta) / (X2.d^4 + beta*X2.d^2)
  delta1 <- sqrt(1 - c1)
  delta2 <- sqrt(1 - c2)
  sqrt.delta1 <- sqrt(1 - sqrt(c1))
  sqrt.delta2 <- sqrt(1 - sqrt(c2))
  bound <- 0
  diags <- list()
  for (i in 1:(rj + ri1)) {
    for (j in 1:(rj + ri2)) {
      if (i == j & i < rj){
        ksi <- sqrt(2) * (sqrt.delta1[i] + 
                            sqrt.delta2[i] + 
                            sqrt(2) * sqrt.delta1[i] * sqrt.delta2[i])
        bound = bound + 2 * ksi + ksi^2
      } else{
        eta = delta1[i] + delta2[j] + delta1[i] * delta2[j]
        diags[[paste(i - j)]] = c(diags[[paste(i - j)]], eta)
      }
    }
  }
  for (name in names(diags)) {
    bound = bound + max(diags[[name]])
  }
  
  return(bound)
}

bound.approx <- function(min.X1, min.X2, rj, ri1, ri2, beta){
  c1.min <- (min.X1 ^4 - beta) / (min.X1 ^ 4 + beta * min.X1^2)
  c2.min <- (min.X2 ^4 - beta) / (min.X2 ^ 4 + beta * min.X2^2)

  ksi <- sqrt(2) * (sqrt(1 - sqrt(c1.min)) + 
                      sqrt(1 - sqrt(c2.min)) + 
                      sqrt(2) * sqrt(1 - sqrt(c1.min)) * sqrt(1 - sqrt(c2.min)))
  bound <- rj * (2*ksi + ksi^2)

  eta <- sqrt(1 - c1.min) + sqrt(1 - c2.min) + sqrt(1 - c1.min) * sqrt(1 - c2.min)
  bound <- bound + (2*rj + ri1 + ri2 - 1) * eta
  
  return(bound)
}

bound.approx.angle <- function(X1.d, X2.d, rj, ri1, ri2, beta, phi.min, phi.max){
  c1.min <- min((X1.d^4 - beta) / (X1.d^4 + beta*X1.d^2))
  c2.min <- min((X2.d^4 - beta) / (X2.d^4 + beta*X2.d^2))
  
  ksi <- sqrt(2) * (sqrt(1 - sqrt(c1.min)) + 
                      sqrt(1 - sqrt(c2.min)) + 
                      sqrt(2) * sqrt(1 - sqrt(c1.min)) * sqrt(1 - sqrt(c2.min)))
  bound <- rj * (2*ksi + ksi^2)
  
  eta <- phi.max * sqrt(c1.min) * sqrt(c2.min) + 
    sqrt(1 - phi.min^2) * (sqrt(1 - c1.min) + sqrt(1 - c2.min)) + 
    sqrt(1 - c1.min) * sqrt(1 - c2.min)
  bound <- bound + (2*rj + ri1 + ri2 - 1) * eta
  
  return(bound)
}

