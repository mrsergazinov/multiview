my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(foreach)
library(doParallel)
library(denoiseR)
library(MASS)
library(Ckmeans.1d.dp)
library(ajive)
library(nnet)

source('src/generate_data_2_views.R')
source('src/utils.R')
source('src/models_2_views.R')
source('src/metrics.R')
load("data/COADdata.rda")

# clean response data: delete 0 valued data and create class idx matrix
# all available 
ID1 = (!is.na(COAD$subtype)) & (!apply(COAD$X1, 1, function(o) sum(is.na(o)))) & 
  (!apply(COAD$X2, 1, function(o) sum(is.na(o))))
# check the numbers of obeservations for each case
sum(ID1 == T) # 167 - no missing

# extract subjects that have no missing information
rnaData = COAD$X1[ID1, ]
mRNAData = COAD$X2[ID1, ]
predictor = list(rnaData, mRNAData)
typeData = COAD$subtype[ID1]

result = matrix(nrow=100, ncol=3)
for (i in 1:100){
  # split into training and test data
  set.seed(i)
  id = sample(1:167, 167)
  trainID = id[1:132] # 80% of the data 
  testID = id[133:167] # 20% of the data
  # training data
  trainData = list(rnaData[trainID, ], mRNAData[trainID, ])
  trainType = typeData[trainID]
  # test data
  testData = list(rnaData[testID, ], mRNAData[testID, ])
  testType = typeData[testID]
  
  # scale the data
  mu <- lapply(trainData, function(x) apply(x, 2, mean))
  s <- lapply(trainData, function(x) apply(x, 2, sd))
  trainData <- lapply(1:2, function(d) scale(trainData[[d]], center = mu[[d]], scale = s[[d]]))
  testData <- lapply(1:2, function(d) scale(testData[[d]], center = mu[[d]], scale = s[[d]]))
  
  # estimate rank
  rank1 <- optishrink(trainData[[1]])$nb.eigen
  rank2 <- optishrink(trainData[[2]])$nb.eigen
  
  # supervised PCA
  trainPCA <- list(
    trainData[[1]], trainData[[2]]
    # plsda(trainData[[1]], trainType, ncomp = rank1),
    # plsda(trainData[[2]], trainType, ncomp = rank2)
  )
  testPCA <- list(
    testData[[1]], testData[[2]]
    # predict(trainPCA[[1]], testData[[1]])$variates,
    # predict(trainPCA[[2]], testData[[2]])$variates
  )
  # trainPCA <- lapply(trainPCA, function(model) model$variates$X)
  
  
  # Proposed
  out <- proposed_func(trainPCA[[1]], trainPCA[[2]], rank1, rank2, return_scores = TRUE)
  trainPCA.joint <- lapply(trainPCA, function(x) out$joint)
  trainPCA.indiv1 <- lapply(trainPCA, function(x) out$indiv1)
  trainPCA.indiv2 <- lapply(trainPCA, function(x) out$indiv2)
  testPCA.joint <- list()
  testPCA.indiv1 <- list()
  testPCA.indiv2 <- list()
  for (j in 1:2) {
    inv.X <- ginv(trainPCA[[j]])
    testPCA.joint[[j]] <- testPCA[[j]] %*% inv.X %*% out$joint 
    testPCA.indiv1[[j]] <- testPCA[[j]] %*% inv.X %*% out$indiv1
    testPCA.indiv2[[j]] <- testPCA[[j]] %*% inv.X %*% out$indiv2
  }
  
  # train models
  data <- data.frame(subtype = trainType, cbind(trainPCA.joint[[1]], trainPCA.joint[[2]]))
  test_data <- data.frame(subtype = testType, cbind(testPCA.joint[[1]], testPCA.joint[[2]]))
  model <- multinom(subtype ~ ., data = data)
  predicted_labels <- predict(model, test_data, type = "class")
  error_joint = sum(predicted_labels != testType) / length(testType)
  
  data <- data.frame(subtype = trainType, cbind(trainPCA.indiv1[[1]], trainPCA.indiv1[[2]]))
  test_data <- data.frame(subtype = testType, cbind(testPCA.indiv1[[1]], testPCA.indiv1[[2]]))
  model <- multinom(subtype ~ ., data = data)
  predicted_labels <- predict(model, test_data, type = "class")
  error_indiv1 = sum(predicted_labels != testType) / length(testType)
  
  data <- data.frame(subtype = trainType, cbind(trainPCA.indiv2[[1]], trainPCA.indiv2[[2]]))
  test_data <- data.frame(subtype = testType, cbind(testPCA.indiv2[[1]], testPCA.indiv2[[2]]))
  model <- multinom(subtype ~ ., data = data)
  predicted_labels <- predict(model, test_data, type = "class")
  error_indiv2 = sum(predicted_labels != testType) / length(testType)
  
  result[i, ] <- c(error_joint, error_indiv1, error_indiv2)
}

prod.sing.vals <- svd(out$test$prod)$d
p <- ggplot(data.frame(sing.vals = prod.sing.vals), aes(x = sing.vals)) +
  geom_histogram(binwidth = 0.1, aes(y = after_stat(count / sum(count)))) +
  ggtitle('Spectrum of Product of Projections') +
  xlab('Singular values') +
  ylab('Frequency') + 
  geom_vline(xintercept = out$test$lam, linetype = 'dashed', color='red') +
  theme_minimal()