unifac.plus.given <- function(X0,p.ind,p.ind.list,n,num.comp=20,max.iter=500,conv.thresh=0.001){
  library(PRIMME)
  
  n.source <- length(p.ind)
  X0.resid <- X0
  S <- list()
  pen <- c()
  max.comb <- length(p.ind.list)
  for(i in 1:max.comb){
    S[[i]]=0
    pen[i]=0}
  obj.vec <- c(sum(X0^2))
  #obj.vec <- c()
  temp.fac <- svd(X0)$d[1]/(sum(sqrt(dim(X0))))-1
  obj.prev <- 10^10
  for(jj in 1:max.iter){
    # print(jj)
    if(jj < 100){
      lambda <- 1+(99-jj)/100*temp.fac
    }
    for(k in 1:max.comb){
      X0.temp <- X0.resid+S[[k]]
      X0.r <- crossprod(X0.temp[p.ind.list[[k]],])
      # print(paste('*',k))
      a <- sqrt(svd(X0.r)$d[1:num.comp])
      nc <- sum(a>(lambda*(sqrt(length(p.ind.list[[k]]))+sqrt(n))))
      S[[k]] <- array(rep(0,prod(dim(X0))),dim=dim(X0))
      pen[k] = 0
      #     Est =
      if(nc>0){
      SVD <- svd(X0.temp[p.ind.list[[k]],], nu=nc,nv=nc)
      s.vals <- pmax(SVD$d[1:nc]-lambda*(sqrt(length(p.ind.list[[k]]))+sqrt(n)),0)
      Diag <- diag(s.vals,nrow=nc,ncol=nc)
      Est <- SVD$u%*%Diag%*%t(SVD$v)
      S[[k]] <- array(rep(0,prod(dim(X0))),dim=dim(X0))
      S[[k]][p.ind.list[[k]],] <- Est
      pen[k] <- lambda*(sqrt(n)+sqrt(length(p.ind.list[[k]])))*(sum(s.vals))
      }
      X0.resid = X0-Reduce('+',S)
      obj.vec <- c(obj.vec,sum(X0.resid^2)+2*sum(pen))
      #     if(abs(obj.vec[jj+1]-obj.vec[jj])<0.001) break
    }
    obj.cur <- sum(X0.resid^2)+2*sum(pen)
    if(abs(obj.cur-obj.prev)<conv.thresh) break
    obj.prev <- sum(X0.resid^2)+2*sum(pen)
  }
  
  Sums <- matrix(nrow=max.comb,ncol=n.source)
  for(kk in 1:max.comb){for(j in 1:n.source){
    Sums[kk,j] = sum(norm(S[[kk]][p.ind[[j]],])^2)
  }}
  return(list(S=S,p.ind.list=p.ind.list,Sums=Sums,obj.vec=obj.vec))
}
