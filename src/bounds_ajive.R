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
  wedin_bound_samples <- mapply(function(u, v)  min(max(u, v)/sigma_min, 1)^2, U_sampled_norms, V_sampled_norms)
  
  wedin_bound_samples
}


wedin_bound_resampling <- function(X, perp_basis, right_vectors, num_samples=1000){
  
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

compute_spectrum_ajive <- function(Y1, Y2, rank1, rank2) {
  U1 <- svd(Y1)$u[, 1:rank1, drop = FALSE]
  U2 <- svd(Y2)$u[, 1:rank2, drop = FALSE]
  M <- cbind(U1, U2)
  M_svd <- get_svd(M)
  return (M_svd[['d']])
}

compute_bounds_ajive <- function(Y1, Y2, rank1, rank2, n_wedin_samples=100, n_rand_dir_samples=100){
  n_obs <- nrow(Y1)
  blocks <- list(Y1, Y2)
  block_svd <- lapply(blocks, get_svd)
  initial_signal_ranks <- c(rank1, rank2)
  K <- 2
  
  # compute Wedin bound
  block_wedin_samples <- matrix(NA, K, n_wedin_samples)  
  for(k in 1:K){
    block_wedin_samples[k, ] <- get_wedin_bound_samples(X=blocks[[k]],
                                                        SVD=block_svd[[k]],
                                                        signal_rank=initial_signal_ranks[k],
                                                        num_samples=n_wedin_samples)
  }
  wedin_samples <-  K - colSums(block_wedin_samples)
  wedin_svsq_threshold <- quantile(wedin_samples, .05)
  
  # compute random directin bound
  rand_dir_samples <- get_random_direction_bound(n_obs=n_obs, 
                                                 dims=initial_signal_ranks, 
                                                 num_samples=n_rand_dir_samples)
  rand_dir_svsq_threshold <- quantile(rand_dir_samples, .95)
  
  return (list(
    "random_directions_bound" = sqrt(rand_dir_svsq_threshold),
    "wedin_bound" = sqrt(wedin_svsq_threshold)
  ))  
} 