library(tidyverse)
library(ajive)

ajive_wrapper <- function(X1, X2, args){
  r1 <- args$rj + args$ri1
  r2 <- args$rj + args$ri2
  error1 <- sample(-1:1, 1)
  error2 <- sample(-1:1, 1)
  ajive_out <- ajive(list(X1, X2), c(r1 + error1, r2 + error2),
                    n_wedin_samples = 100,
                    n_rand_dir_samples = 100)
  return(list("joint" = ajive_out$joint_scores,
              "indiv1" = ajive_out$block_decomps[[1]][['individual']][['u']],
              "indiv2" = ajive_out$block_decomps[[2]][['individual']][['u']]))
}