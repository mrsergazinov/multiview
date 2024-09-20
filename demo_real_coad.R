my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(reticulate)
library(denoiseR)
library(MASS)
library(ggplot2)

library(nnet)

library(ajive)
library(SLIDE)
library(r.jive)
library(PRIMME)
library(pracma)
library(Ckmeans.1d.dp)

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
# center and scale
rnaData = scale(rnaData)
mRNAData = scale(mRNAData)
data = list(Y1 = rnaData, 
            Y2 = mRNAData)
type = COAD$subtype[ID1]

# estimate marginal ranks
rank1 <- optishrink(rnaData)$nb.eigen
rank2 <- optishrink(mRNAData)$nb.eigen

# extract joint and individual 
models <- c("naive")
            # "proposed", "slide", "jive", "ajive", "unifac", "dcca")
naive <- function(Y1, Y2, rank1, rank2, return_scores = FALSE){
  joint <- matrix(1, nrow = nrow(Y1), ncol = nrow(Y1))
  indiv1 <- svd(Y1)$u[, 1:rank1]
  indiv2 <- svd(Y2)$u[, 1:rank2]
  return(list(joint = joint, indiv1 = indiv1, indiv2 = indiv2))
}
compute <- list(naive = naive,
                slide = slide_func,
                jive = jive_func,
                ajive = ajive_func,
                dcca = dcca_func,
                unifac = unifac_func,
                proposed = proposed_func)
for (model in models) {
  out <- compute[[model]](data$Y1, data$Y2, rank1, rank2, return_scores = TRUE)
  data[[paste0(model, "_joint")]] <- out$joint
  data[[paste0(model, "_indiv1")]] <- out$indiv1
  data[[paste0(model, "_indiv2")]] <- out$indiv2
}
# save(data, file = "data/COADdata_processed.rda")
# load data
# load("data/COADdata_processed.rda")

n_sim <- 100
results <- list()
for (model in models){
  results[[model]] <- matrix(NA, nrow = n_sim, ncol = 3)
}
for (i in 1:n_sim){
  # split into training and test data
  set.seed(i)
  id = sample(1:167, 167)
  trainID = id[1:132] # 80% of the data 
  testID = id[133:167] # 20% of the data
  
  # run models
  for (model in models){
    model_results <- c()
    for (component in c("_joint", "_indiv1", "_indiv2")){
      train_data <- data.frame(subtype = type[trainID], data[[paste0(model, component)]][trainID, ])
      test_data <- data.frame(subtype = type[testID], data[[paste0(model, component)]][testID, ])
      colnames(train_data) <- colnames(test_data)
      predictive_model <- multinom(subtype ~ ., data = train_data)
      predicted_labels <- predict(predictive_model, test_data, type = "class")
      error <- sum(predicted_labels != type[testID]) / length(type[testID])
      model_results <- c(model_results, error)
    }
    results[[model]][i, ] <- model_results
  }
}

for (model in models){
  means <- colMeans(results[[model]], na.rm = TRUE)
  sds <- apply(results[[model]], 2, sd, na.rm = TRUE)
  str <- paste0(model)
  for (i in 1:3){
    str <- paste0(str, " & ", round(means[i], 2), " (", round(sds[i], 2), ")")
  }
  print(str)
}