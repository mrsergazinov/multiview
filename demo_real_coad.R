my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(foreach)
library(doParallel)
library(ggplot2)
source('src/utils.R')
source('src/models_2_views.R')
source('src/regression_2_views.R')
load("data/COADdata.rda")


##########################################################
#----------- Data preprocessing.-------------------------#
##########################################################

# clean response data: delete 0 valued data and create class idx matrix
ID1 = (!is.na(COAD$subtype)) & (!apply(COAD$X1, 1, function(o) sum(is.na(o)))) & 
  (!apply(COAD$X2, 1, function(o) sum(is.na(o))))
# check the numbers of obeservations for each case
sum(ID1 == T) # 167 - no missing
# extract subjects that have no missing information
rnaData = COAD$X1[ID1, ]
mRNAData = COAD$X2[ID1, ]
# center and scale
# rnaData = scale(rnaData)
# mRNAData = scale(mRNAData)
data = list(Y1 = rnaData, 
            Y2 = mRNAData)
type = COAD$subtype[ID1]
# estimate marginal ranks
rank1 <- 16
# optishrink(rnaData)$nb.eigen
rank2 <- 16
  # optishrink(mRNAData)$nb.eigen


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
  # id = sample(1:167, 167)
  trainID = 1:167 
  testID = 1:167

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
      train_data <- data.frame(subtype = type[trainID], out[[paste0("train_", component)]])
      test_data <- data.frame(subtype = type[testID], out[[paste0("test_", component)]])
      colnames(train_data) <- colnames(test_data)
      predictive_model <- multinom(subtype ~ ., data = train_data)
      predicted_labels <- predict(predictive_model, test_data, type = "class")
      error <- sum(predicted_labels != type[testID]) / length(type[testID])
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
save(collapsed, file = "data/COADdata_results.rda")


str <- ''
for (method_name in names(collapsed)){
  results <- collapsed[[method_name]]
  ranks <- results[1:3]
  angles <- results[4:6]
  errors <- results[7:9]
  str <- paste0(str, method_name, 
                 ' & ', ranks[1], ' & ', round(errors[1]*100, 1),
                 ' & ', ranks[2], ' & ', round(errors[2]*100, 1),
                 ' & ', ranks[3], ' & ', round(errors[3]*100, 1),
                ' & ', round(angles[1], 0), ' & ', round(angles[2], 0), ' & ', round(angles[3], 0 ),
                '\n')
}
cat(str)

