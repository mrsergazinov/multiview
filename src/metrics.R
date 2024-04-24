compute_fd <- function(P1, P2){
  sum(diag((diag(dim(P2)[1]) - P2) %*% P1))
}
compute_tp <- function(P1, P2){
  sum(diag(P1 %*% P2))
}