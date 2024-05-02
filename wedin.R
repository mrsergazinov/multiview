my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(foreach)
library(doParallel)
library(pracma)
source('src/generate_data_2_views.R')
source('src/utils.R')
source('src/models_2_views.R')
source('src/metrics.R')

rj <- 4
ri1 <- 5
ri2 <- 4
m <- 50
phi_max <- 30
n1 <- 80
n2 <- 100
snr1 <- 2
snr2 <- 2
signal_strength1 <- 10
signal_strength2 <- 12
sigma1 <- (signal_strength1 / snr1) / (sqrt(m) + sqrt(n1))
sigma2 <- (signal_strength2 / snr2) / (sqrt(m) + sqrt(n2))
rank_spec <- 0
no_joint <- FALSE
no_indiv <- FALSE

data <- generate_data(m, n1, n2, 
                      rj, ri1, ri2, rank_spec, 
                      signal_strength1, signal_strength2, 
                      sigma1, sigma2,
                      no_joint, no_indiv, 
                      phi_max)
Y1 <- data$Y1
Y2 <- data$Y2
rank1 <- data$rank1
rank2 <- data$rank2

get_wedin_bound_samples <- function(X, SVD, signal_rank, num_samples=1000){
  
  # resample for U and V
  U_perp <- SVD[['u']][ , -(1:signal_rank)]
  U_sampled_norms <- wedin_bound_resampling(X=X,
                                            perp_basis=U_perp,
                                            right_vectors=FALSE,
                                            num_samples=num_samples)
  
  V_perp <- SVD[['v']][ , -(1:signal_rank)]
  V_sampled_norms <- wedin_bound_resampling(X=X,
                                            perp_basis=V_perp,
                                            right_vectors=TRUE,
                                            num_samples=num_samples)
  
  sigma_min <- SVD[['d']][signal_rank]
  wedin_bound_samples <- mapply(function(u, v)  min(max(u, v)/sigma_min, 1), U_sampled_norms, V_sampled_norms)
  
  wedin_bound_samples
}
wedin_bound_resampling <- function(X, perp_basis, right_vectors, num_samples=100){
  
  rank <- dim(perp_basis)[2]
  resampled_norms <- rep(0, num_samples)
  
  for(s in 1:num_samples){
    
    sampled_col_index <- sample.int(n=dim(perp_basis)[2],
                                    size=rank,
                                    replace=TRUE)
    
    
    perp_resampled <- perp_basis[ , sampled_col_index]
    
    if(right_vectors){
      resampled_projection <- X %*% perp_resampled
    } else{
      resampled_projection <- t(perp_resampled) %*% X
    }
    
    # operator L2 norm
    resampled_norms[s] <- norm(resampled_projection,
                               type='2')
  }
  
  resampled_norms
}
get_svd <- function(X, rank=NULL){
  # SVD <- get_svd(X, rank=2)
  
  if(is.null(rank)){
    svd(X)
  } else if(rank == 0){
    # TODO: what to do
    decomposition <- list()
    decomposition[['u']] <- matrix(0, ncol=1, nrow=dim(X)[1])
    decomposition[['d']] <- 0
    decomposition[['v']] <- matrix(0, ncol=1, nrow=dim(X)[2])
    decomposition
    
  } else{
    decomposition <- svd(X, nu=rank, nv=rank)
    decomposition[['d']] <- decomposition[['d']][1:rank]
    decomposition
  }
  
}
get_random_direction_bound <- function(n_obs, dims, num_samples=1000){
  
  n_blocks <- length(dims)
  rand_dir_samples <- rep(0, num_samples)
  for(s in 1:num_samples){
    rand_subspaces <- list()
    for(b in 1:n_blocks){
      X <- matrix(rnorm(n_obs * dims[b], mean=0,sd=1), n_obs, dims[b])
      U <- get_svd(X)[['u']]
      rand_subspaces[[b]] <- U
      
    }
    M <- do.call(cbind, rand_subspaces)
    M_svd <- get_svd(M, rank=min(dims))
    
    rand_dir_samples[s] <- M_svd[['d']][1]^2
    
  }
  
  rand_dir_samples
}
get_random_direction_bound.analytical <- function(n_obs, dims) {
  m <- n_obs
  rank1 <- dims[1]
  rank2 <- dims[2]
  q1 <- rank1 / m
  q2 <- rank2 / m
  q.plus <- q1 + q2 - 2*q1*q2 + 2*sqrt(q1*q2*(1-q1)*(1-q2))
  return(q.plus)
}

svd.Y1 <- svd(Y1)
svd.Y2 <- svd(Y2)
U1.hat <- svd.Y1$u[, 1:rank1, drop = FALSE]
U2.hat <- svd.Y2$u[, 1:rank2, drop = FALSE]
svd.prod <- svd(U1.hat %*% t(U1.hat) %*% U2.hat %*% t(U2.hat))$d
print(paste0("prod svd: "))
print(svd.prod[1:10])
print(paste0("rand bound:", sqrt(rand.bound.analytical)))

out1 <- get_wedin_bound_samples(X=Y1, SVD=svd(Y1), signal_rank=rank1)
out2 <- get_wedin_bound_samples(X=Y2, SVD=svd(Y2), signal_rank=rank2)
sin1 <- quantile(out1, probs=0.05)
sin2 <- quantile(out2, probs=0.05)
print(paste0("Wedin-based bound:", 2 - sin1^2 - sin2^2))

rand.bound.sampling <- quantile(get_random_direction_bound(m, c(rank1, rank2)), probs=0.95)
print(paste0("Estimated rand dir bound: ", rand.bound.sampling))

rand.bound.analytical <- get_random_direction_bound.analytical(m, c(rank1, rank2))
print(paste0("Analytical rand dir bound: ", sqrt(rand.bound.analytical) + 1))

print("Sq. singular values from AJIVE: ")
print(svd(cbind(U1.hat, U2.hat))$d^2)


