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
  
  # Proposed
  out <- proposed_func(trainData[[1]], trainData[[2]], rank1, rank2, return_scores = TRUE)
  W1.joint <- ginv(trainData[[1]]) %*% out$joint
  W2.joint <- ginv(trainData[[2]]) %*% out$joint
  W1.indiv <- ginv(trainData[[1]]) %*% out$indiv1
  W2.indiv <- ginv(trainData[[2]]) %*% out$indiv2
  
  # Joint
  data <- data.frame(subtype = trainType, cbind(out$joint, out$joint))
  test_data = data.frame(subtype = testType, cbind(testData[[1]] %*% W1.joint, testData[[2]] %*% W2.joint))
  model <- multinom(subtype ~ ., data = data)
  predicted_labels <- predict(model, test_data, type = "class")
  error_joint = sum(predicted_labels != testType) / length(testType)
  
  # Individual
  data <- data.frame(subtype = trainType, out$indiv1)
  test_data = data.frame(subtype = testType, testData[[1]] %*% W1.indiv)
  model <- multinom(subtype ~ ., data = data)
  predicted_labels <- predict(model, test_data, type = "class")
  error_indiv1 = sum(predicted_labels != testType) / length(testType)
  
  data <- data.frame(subtype = trainType, out$indiv2)
  test_data = data.frame(subtype = testType, testData[[2]] %*% W2.indiv)
  model <- multinom(subtype ~ ., data = data)
  predicted_labels <- predict(model, test_data, type = "class")
  error_indiv2 = sum(predicted_labels != testType) / length(testType)
  
  # Save errors
  result[i, ] <- c(error_joint, error_indiv1, error_indiv2)
}