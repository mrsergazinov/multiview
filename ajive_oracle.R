library(tidyverse)
library(ajive)

ajive_oracle_wrapper <- function(X1, X2, rj, ri1, ri2){
  ajive_out <- ajive(list(X1, X2), c(rj + ri1, rj+ ri2),
                    n_wedin_samples = 100, 
                    n_rand_dir_samples = 100,
                    joint_rank = rj)
  return(ajive_out$joint_scores)
}