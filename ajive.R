library(tidyverse)
library(ajive)

ajive_wrapper <- function(X1, X2, args){
  ajive_out <- ajive(list(X1, X2), c(args$rj + args$ri1 + 10, args$rj + args$ri2 + 10),
                    n_wedin_samples = 100,
                    n_rand_dir_samples = 100)
  return(ajive_out$joint_scores)
}