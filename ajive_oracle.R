library(tidyverse)
library(ajive)

ajive_oracle_wrapper <- function(X1, X2, args){
  ajive_out <- ajive(list(X1, X2), c(args$rj + args$ri1, args$rj + args$ri2),
                    n_wedin_samples = 100, 
                    n_rand_dir_samples = 100,
                    joint_rank = args$rj)
  return(ajive_out$joint_scores)
}