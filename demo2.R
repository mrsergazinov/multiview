suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

# define number of cores and start parallel backend
set.seed(1234)
numCores <- detectCores() - 1  # Leave one core for system processes
cl <- makeCluster(numCores)
registerDoParallel(cl)

# define other params
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
rank_spec <- 'exact'
no_joint <- FALSE
no_indiv <- FALSE
save_file <- paste0("demo2_", rank_spec, "_", no_joint, "_", no_indiv)
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

sim_iter <- 50
models <- c("jive", "ajive", "dcca", "slide", "unifac", "proposed", "proposed_subsampling")
iters <- foreach(i = 1:sim_iter) %dopar% {
  # load models
  source("src/models_2_views.R")
  # compute error
  source("src/metrics.R")
  # data generator
  source("src/generate_data_2_views.R")
  data <- generate_data(m, n1, n2, 
                        rj, ri1, ri2, rank_spec, 
                        signal_strength1, signal_strength2, 
                        sigma1, sigma2,
                        no_joint, no_indiv, 
                        phi_max)
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
# process and save
results <- lapply(results, function(x) colMeans(x))
results <- do.call(rbind, results)
results <- as.data.frame(results)
colnames(results) <- c("fdr.P1", "tpr.P1",
                      "fdr.P2", "tpr.P2",
                      "fdr.Pjoint", "tpr.Pjoint",
                      "fdr.Pindiv1", "tpr.Pindiv1",
                      "fdr.Pindiv2", "tpr.Pindiv2")
# compute precision = 1 - fdr
results <- results %>% mutate(precision.Pjoint = 1 - fdr.Pjoint,
                              precision.Pindiv1 = 1 - fdr.Pindiv1,
                              precision.Pindiv2 = 1 - fdr.Pindiv2)
# avg precision and avg tpr
results$avg_precision <- (results$precision.Pjoint +
                            results$precision.Pindiv1 +
                            results$precision.Pindiv2) / 3
results$avg_tpr <- (results$tpr.Pjoint +
                      results$tpr.Pindiv1 +
                      results$tpr.Pindiv2) / 3
# compute F1 score
results$f1 <- 2 * (results$avg_precision * results$avg_tpr) / (results$avg_precision + results$avg_tpr)
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
save(results.save, file=paste0("results/demo2_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RData"))
