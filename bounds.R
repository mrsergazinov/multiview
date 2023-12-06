bound <-function(X1.d, X2.d, rj, ri1, ri2, beta){
  c1 <- (X1.d^4 - beta) / (X1.d^4 + beta*X1.d^2)
  c2 <- (X2.d^4 - beta) / (X2.d^4 + beta*X2.d^2)
  delta1 <- sqrt(1 - c1)
  delta2 <- sqrt(1 - c2)
  bound <- c()
  for (i in 1:(rj + ri1)) {
    for (j in 1:(rj + ri2)) {
      if (i != j) {
        bound = c(bound, (delta1[i] + delta2[j] + delta1[i] * delta2[j]))
      }
    }
  }
  bound <- 2 * max(bound)
  for (i in 1:(rj + ri1)) {
    for (j in 1:(rj + ri2)) {
      if (i == j & i > rj) {
        bound = bound + (delta1[i] + delta2[j] + delta1[i] * delta2[j])
      }
    }
  }
  delta1 <- sqrt(1 - sqrt(c1))
  delta2 <- sqrt(1 - sqrt(c2))
  for (i in 1:rj) {
    ksi <- sqrt(2) * (delta1[i] + delta2[i] + sqrt(2) * delta1[i] * delta2[i])
    bound <- bound + 2 * ksi + ksi^2
  }
  
  return(bound)
}