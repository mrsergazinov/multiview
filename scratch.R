
d <- c()
r1 <- 3
r2 <- 3
u <- rep(1, 100)
v <- rep(1, 40)
L <- u %*% (10*diag(1)) %*% t(v)
for (j in 1:10){
  avgP <- 0
  for (i in 1:100) {
    U <- svd(L + matrix(rnorm(400), 100, 40))$u
    V <- svd(L + matrix(rnorm(400), 100, 40))$u
    P = U[,1:r1]%*%t(U[,1:r1])
    Q = V[,1:r2]%*%t(V[,1:r2])
    avgP = avgP + (P%*%Q + Q%*%P) / 2
  }
  avgP <- avgP / 100
  d <- c(d, svd(avgP)$d)
}
