my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(foreach)
library(doParallel)
source('src/generate_data_2_views.R')
source('src/models_2_views.R')
source('src/metrics.R')

# file_paths <- list.files(path = './results', pattern = "^demo2.*\\.RData$", full.names = TRUE)
# 
# for (file in file_paths) {
#   load(file)
#   rj <- results.save$rj
#   ri1 <- results.save$ri1
#   ri2 <- results.save$ri2
#   m <- results.save$m
#   phi_max <- results.save$phi_max
#   n1 <- results.save$n1
#   n2 <- results.save$n2
#   snr1 <- results.save$SNR1
#   snr2 <- results.save$SNR2
#   signal_strength1 <- results.save$signal_strength1
#   signal_strength2 <- results.save$signal_strength2
#   sigma1 <- results.save$sigma1
#   sigma2 <- results.save$sigma2
#   rank_spec <- results.save$rank_spec
#   no_joint <- results.save$no_joint
#   no_indiv <- results.save$no_indiv
#   
#   print(paste0("snr1: ", snr1))
#   print(paste0("snr2: ", snr2))
#   print(paste0("phi_max: ", phi_max))
#   print(paste0("no joint: ", no_joint))
#   print(paste0("no indiv: ", no_indiv))
#   print(paste0("rank_spec: ", rank_spec))
  
# define other params
rj <- 4
ri1 <- 5
ri2 <- 4
m <- 50
phi_max <- 60
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
try(if (no_joint && no_indiv) stop("At least one of no_joint and no_indiv must be FALSE"))

# set args from command line
args <- commandArgs(trailingOnly = TRUE)
# Check if arguments are provided
if (length(args) > 0) {
  # Parse command-line arguments
  for (i in seq_along(args)) {
    print(args[[i]])
    eval(parse(text = args[[i]]))
  }
}
save_file <- paste0("results/demo2_", rank_spec, "_", no_joint, "_", no_indiv)

# define number of cores and start parallel backend
set.seed(1234)
numCores <- detectCores()-1  # Leave one core for system processes
cl <- makeCluster(numCores)
clusterEvalQ(cl, .libPaths("./multiview_rlibs"))
registerDoParallel(cl)

sim_iter <- 50
models <- c("ajive", "proposed_subsampling")
  # c("slide", "jive", "ajive", "dcca", "unifac", "proposed", "proposed_subsampling")
packages <- c('reticulate', 'ajive', 'r.jive', 'SLIDE', 'Ckmeans.1d.dp', 'pracma', 'PRIMME')
iters <- foreach(i = 1:sim_iter, .packages=packages) %dopar% {
# library(ajive)
# library(Ckmeans.1d.dp)
# library(pracma)
# for (i in 1:sim_iter) {
  compute <- list(slide = slide_func, 
                  jive = jive_func, 
                  ajive = ajive_func, 
                  dcca = dcca_func, 
                  unifac = unifac_func, 
                  proposed = proposed_func, 
                  proposed_subsampling = proposed_subsampling_func)
  
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
  # compute results
  result <- list()
  for (model in models) {
    out <- compute[[model]](data$Y1, data$Y2, data$rank1, data$rank2)
    out <- c(compute_fd(out$P1, data$P1) / out$r1, compute_tp(out$P1, data$P1) / data$rank1,
             compute_fd(out$P2, data$P2) / out$r2, compute_tp(out$P2, data$P2) / data$rank2,
             compute_fd(out$Pjoint, data$Pjoint) / out$rj, compute_tp(out$Pjoint, data$Pjoint) / data$rj,
             compute_fd(out$Pindiv1, data$Pindiv1) / out$ri1, compute_tp(out$Pindiv1, data$Pindiv1) / data$ri1,
             compute_fd(out$Pindiv2, data$Pindiv2) / out$ri2, compute_tp(out$Pindiv2, data$Pindiv2) / data$ri2)
    # check if any NA
    if (any(is.na(out))){
      print(paste0("NA detected in ", model))
      break
    }
    result[[model]] <- matrix(out, nrow=1)
  }
  result
}
# shut down the parallel backend
stopCluster(cl)
# collect results
results <- list()
for (iter in iters) {
  for (model in models) {
    results[[model]] = rbind(results[[model]], iter[[model]])
  }
}
for (model in models) {
  results[[model]] <- as.data.frame(results[[model]])
  colnames(results[[model]]) <- c("fdr.P1", "tpr.P1",
                                  "fdr.P2", "tpr.P2",
                                  "fdr.Pjoint", "tpr.Pjoint",
                                  "fdr.Pindiv1", "tpr.Pindiv1",
                                  "fdr.Pindiv2", "tpr.Pindiv2")
  
}
# save results
results.save <- list()
results.save[["results"]] <- results
results.save[["sim_iter"]] <- sim_iter
results.save[["rj"]] <- rj
results.save[["ri1"]] <- ri1
results.save[["ri2"]] <- ri2
results.save[["m"]] <- m
results.save[["phi_max"]] <- phi_max
results.save[["n1"]] <- n1
results.save[["n2"]] <- n2
results.save[["SNR1"]] <- snr1
results.save[["SNR2"]] <- snr2
results.save[["sigma1"]] <- sigma1
results.save[["sigma2"]] <- sigma2
results.save[["signal_strength1"]] <- signal_strength1
results.save[["signal_strength2"]] <- signal_strength2
results.save[["rank_spec"]] <- rank_spec
results.save[["no_joint"]] <- no_joint
results.save[["no_indiv"]] <- no_indiv
save(results.save, file=paste0(save_file, format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData"))

# results.save$results[['proposed_subsampling_new']] <- results[['proposed_subsampling']]
# save(results.save, file=file)
# }
