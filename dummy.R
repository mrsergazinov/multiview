suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

# define number of cores and start parallel backend
set.seed(1234)
numCores <- detectCores() - 1  # Leave one core for system processes
cl <- makeCluster(numCores)
registerDoParallel(cl)

# run parallel
iters <- foreach(i = 1:sim_iter) %dopar% {
  source('src/generate_data_2_views.R')
  source('src/models_2_views.R')
  c(0,0)
}
# shut down the parallel backend
stopCluster(cl)
# collect results
results <- list()