my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(mixOmics)
library(foreach)
library(doParallel)
library(ggplot2)
source('src/utils.R')
source('src/models_2_views.R')
source('src/regression_2_views.R')


##########################################################
#----------- Data preprocessing.-------------------------#
##########################################################
data(nutrimouse)
Y1 <- nutrimouse$gene
Y2 <- nutrimouse$lipid
diet <- nutrimouse$diet
subtype <- nutrimouse$genotype
rank1 <- 3
rank2 <- 4

# process lipids
Y2 <- Y2 * 100
Y2[Y2 == 0] <- 0.375 / (100 + 0.75)

# Standardize columns
Y1 <- scale(Y1)
Y2 <- scale(Y2)
data <- list(Y1 = Y1, Y2 = Y2)

##########################################################
#----------- Model fitting.------------------------------#
##########################################################
nsim <- 1
packages <- c(
  'reticulate', 'ajive','r.jive', 
  'SLIDE','Ckmeans.1d.dp', 'pracma', 
  'PRIMME', 'denoiseR', 'RMTstat', 
  'MASS', 'nnet'
)
numCores <- detectCores()-1  # Leave one core for system processes
cl <- makeCluster(numCores)
clusterEvalQ(cl, .libPaths("./multiview_rlibs"))
registerDoParallel(cl)
iters <- foreach(i = 1:nsim, .packages=packages) %dopar% {
  naive <- function(Y1, Y2, rank1, rank2, return_scores = FALSE){
    joint <- matrix(1, nrow = nrow(Y1), ncol = nrow(Y1))
    indiv1 <- svd(Y1)$u[, 1:rank1]
    indiv2 <- svd(Y2)$u[, 1:rank2]
    return(list(joint = joint, indiv1 = indiv1, indiv2 = indiv2))
  }
  methods <- list(
    naive = naive,
    slide = slide_func,
    jive = jive_func,
    ajive = ajive_func,
    dcca = dcca_func,
    unifac = unifac_func,
    proposed = proposed_func
  )
  # split into training and test data
  set.seed(i)
  trainID = 1:40
  testID = 1:40
  
  # run models
  results <- list()
  for (method_name in names(methods)){
    method <- methods[[method_name]]
    out <- reduce_dimensions(method, 
                             list(data$Y1[trainID, ], data$Y2[trainID, ]), 
                             list(data$Y1[testID, ], data$Y2[testID, ]), 
                             c(rank1, rank2))
    model_results <- c(dim(out$train_joint)[2], dim(out$train_indiv1)[2], dim(out$train_indiv2)[2])
    model_results <- c(model_results,
                       (acos(svd(t(out$train_joint) %*% out$train_indiv1)$d[1]) / pi * 180),
                       (acos(svd(t(out$train_joint) %*% out$train_indiv2)$d[1]) / pi * 180),
                       (acos(svd(t(out$train_indiv1) %*% out$train_indiv2)$d[1]) / pi * 180))
    for (component in c("joint", "indiv1", "indiv2")){
      # genotype
      train_data <- data.frame(subtype = subtype[trainID], out[[paste0("train_", component)]])
      test_data <- data.frame(subtype = subtype[testID], out[[paste0("test_", component)]])
      colnames(train_data) <- colnames(test_data)
      predictive_model <- multinom(subtype ~ ., data = train_data)
      predicted_labels <- predict(predictive_model, test_data, subtype = "class")
      error <- sum(predicted_labels != subtype[testID]) / length(subtype[testID])
      model_results <- c(model_results, error)
      
      # diet
      train_data <- data.frame(diet = diet[trainID], out[[paste0("train_", component)]])
      test_data <- data.frame(diet = diet[testID], out[[paste0("test_", component)]])
      colnames(train_data) <- colnames(test_data)
      predictive_model <- multinom(diet ~ ., data = train_data)
      predicted_labels <- predict(predictive_model, test_data, diet = "class")
      error <- sum(predicted_labels != diet[testID]) / length(diet[testID])
      model_results <- c(model_results, error)
    }
    results[[method_name]] <- model_results
  }
  results
}
# shut down the parallel backend
stopCluster(cl)
# go through list and row bind vectors for each method
collapsed <- list()
for (method_name in names(iters[[1]])) {
  method_results <- lapply(iters, function(sim) sim[[method_name]])
  
  # Combine the list of vectors into a matrix where each row is a simulation iteration
  collapsed[[method_name]] <- do.call(rbind, method_results)
}
save(collapsed, file = "data/MICEdata_results.rda")


str <- ''
for (method_name in names(collapsed)){
  results <- collapsed[[method_name]]
  ranks <- results[1:3]
  angles <- results[4:6]
  errors <- results[7:12]
  str <- paste0(str, method_name, 
                ' & ', ranks[1], ' & ', round(errors[1]*100, 1), ' & ', round(errors[2]*100, 1),
                ' & ', ranks[2], ' & ', round(errors[3]*100, 1), ' & ', round(errors[4]*100, 1),
                ' & ', ranks[3], ' & ', round(errors[5]*100, 1), ' & ', round(errors[6]*100, 1),
                ' & ', round(angles[3], 0 ),
                '\n')
}
cat(str)
