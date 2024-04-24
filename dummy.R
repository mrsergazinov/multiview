my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)
suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})
source('src/generate_data_2_views.R')
source('src/models_2_views.R')
# define number of cores and start parallel backend
set.seed(1234)
numCores <- 48  # Leave one core for system processes
cl <- makeCluster(numCores)
registerDoParallel(cl)
clusterEvalQ(cl, my_lib_path)

# run parallel
packages <- c('reticulate', 'ajive', 'r.jive', 'SLIDE', 'Ckmeans.1d.dp', 'pracma', 'PRIMME')
iters <- foreach(i = 1:10, .packages=packages) %dopar% {
  data <- generate_data(50, 80, 70,
                       2, 2, 2, 0,
                       10, 10,
                       1, 1,
                       FALSE, FALSE,
                       30)
  out <- slide(cbind(data$Y1, data$Y2), c(80, 70))
  out
                                             
}
# shut down the parallel backend
stopCluster(cl)
print(iters)
