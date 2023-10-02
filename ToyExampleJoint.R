# Dimensions fixed
# r = 3
rj = 4
n = 100
p1 = 100
p2 = 150

# Joint + noise, no individual
set.seed(235017)
# Parameters
Uj = svd(matrix(rnorm(n * rj), n, rj))$u
V1 = svd(matrix(rnorm(rj * p1), rj, p1))$v
V2 = svd(matrix(rnorm(rj * p2), rj, p2))$v
dj = c(10, 5, 7, 8) # signal strengths
# Add noise
X1 = Uj %*% diag(dj) %*% t(V1) + matrix(rnorm(n * p1), n, p1)
X2 = Uj %*% diag(dj) %*% t(V2) + matrix(rnorm(n * p2), n, p2)

# Extract joint given the rank
# Full version - not useful
svd1 = svd(X1)
svd2 = svd(X2)

# Reduced version
r1 = 30
r2 = 30

# Align bases by similarity and take the first rcand
svdInner = svd(crossprod(svd1$u[, 1:r1], svd2$u[, 1:r2])) # U1^t U2 = PLR^t, then (P^t U1^t) (U2 R) = L
# Get the rotated bases
U1rotated = svd1$u[, 1:r1] %*% svdInner$u
U2rotated = svd2$u[, 1:r2] %*% svdInner$v
round(crossprod(U1rotated, U2rotated), 2)
# Try candidate joint ranks
rcand = 4
Uj_hat = svd(cbind(U1rotated[, 1:rcand], U2rotated[,1:rcand]))$u[,1:rcand]

hatProj <- Uj_hat %*% t(Uj_hat)
trueProj <- Uj %*% t(Uj)
# compute FD and metric
fd <- sum(diag((diag(n)-trueProj) %*% hatProj))
print(fd)


