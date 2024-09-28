my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(reticulate)
library(denoiseR)
library(RMTstat)
library(MASS)
library(ggplot2)
library(gridExtra)

library(mixOmics)

library(nnet)

library(ajive)
library(SLIDE)
library(r.jive)
library(PRIMME)
library(pracma)
library(Ckmeans.1d.dp)

source('src/utils.R')
source('src/models_2_views.R')

data(nutrimouse)
Y1 <- nutrimouse$gene
Y2 <- nutrimouse$lipid
diet <- nutrimouse$diet
rank1 <- 3
rank2 <- 4

# process lipids
Y2 <- Y2 * 100
Y2[Y2 == 0] <- 0.375 / (100 + 0.75)

# Standardize columns
Y1 <- scale(Y1)
Y2 <- scale(Y2)
data <- list(Y1 = Y1, Y2 = Y2)

# extract joint and individual 
models <- c("naive", "jive", "ajive", "slide",  "dcca", "unifac", "proposed")
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
save.proposed.out <- NA
for (model in models) {
  out <- compute[[model]](data$Y1, data$Y2, rank1, rank2, return_scores = TRUE)
  if (model == "proposed"){
    save.proposed.out <- out
  }
  data[[paste0(model, "_joint")]] <- out$joint
  data[[paste0(model, "_indiv1")]] <- out$indiv1
  data[[paste0(model, "_indiv2")]] <- out$indiv2
}
save(data, file = "data/MICEdata_processed.rda")
# load data
load("data/MICEdata_processed.rda")

n_sim <- 100
results <- list()
for (model in models){
  results[[model]] <- matrix(NA, nrow = n_sim, ncol = 3)
}
for (i in 1:n_sim){
  # split into training and test data
  set.seed(i)
  id = sample(1:40, 40)
  trainID = id[1:35] # 80% of the data
  testID = id[36:40] # 20% of the data
  
  # run models
  for (model in models){
    model_results <- c()
    for (component in c("_joint", "_indiv1", "_indiv2")){
      train_data <- data.frame(diet = diet[trainID], data[[paste0(model, component)]][trainID, ])
      test_data <- data.frame(diet = diet[testID], data[[paste0(model, component)]][testID, ])
      colnames(train_data) <- colnames(test_data)
      predictive_model <- multinom(diet ~ ., data = train_data)
      predicted_labels <- predict(predictive_model, test_data, type = "class")
      error <- sum(predicted_labels != diet[testID]) / length(diet[testID])
      model_results <- c(model_results, error)
    }
    results[[model]][i, ] <- model_results
  }
}

for (model in models){
  ranks <- c(dim(data[[paste0(model, "_joint")]])[2], 
             dim(data[[paste0(model, "_indiv1")]])[2], 
             dim(data[[paste0(model, "_indiv2")]])[2])
  means <- colMeans(results[[model]], na.rm = TRUE)
  sds <- apply(results[[model]], 2, sd, na.rm = TRUE)
  str <- paste0(model)
  for (i in 1:3){
    str <- paste0(str, " & $", ranks[i] , "$ & $", round(means[i], 2), " pm ", round(sds[i], 2), "$")
  }
  print(str)
}

for (model in models){
  angles_joint_indiv1 <- round(acos(svd(t(data[[paste0(model, "_joint")]]) %*% data[[paste0(model, "_indiv1")]])$d) / pi * 180)
  angles_joint_indiv2 <- round(acos(svd(t(data[[paste0(model, "_joint")]]) %*% data[[paste0(model, "_indiv2")]])$d) / pi * 180)
  angles_indiv1_indiv2 <- round(acos(svd(t(data[[paste0(model, "_indiv1")]]) %*% data[[paste0(model, "_indiv2")]])$d) / pi * 180)
  str <- paste0(model)
  str <- paste0(str, " & $[", min(angles_joint_indiv1), "^circ, ", max(angles_joint_indiv1), "^circ]$")
  str <- paste0(str, " & $[", min(angles_joint_indiv2), "^circ, ", max(angles_joint_indiv2), "^circ]$")
  str <- paste0(str, " & $[", min(angles_indiv1_indiv2), "^circ, ", max(angles_indiv1_indiv2), "^circ]$")
  print(str)
}

# Plot proposed
plot.data <- data.frame(diet = nutrimouse$diet, genotype = nutrimouse$genotype)
plot.data[['Joint 1']] <- data$proposed_joint[, 1]
plot.data[['Joint 2']] <- data$proposed_joint[, 2]
plot.data[['Individual 1']] <- data$proposed_indiv2[, 1]
plot.data[['Individual 2']] <- data$proposed_indiv2[, 2]
p1 <- ggplot(plot.data, aes(x = `Joint 1`, y = `Joint 2`, color = diet, shape = genotype)) +
  geom_point(size=3) +
  xlab('Component 1') + 
  ylab('Component 2') +
  ggtitle('Joint View: Gene and Lipid') +
  theme_minimal()

p2 <- ggplot(plot.data, aes(x = `Individual 1`, y = `Individual 2`, color = diet, shape = genotype)) +
  geom_point(size=3) +
  xlab('Component 1') + 
  ylab('Component 2') +
  ggtitle('Individual View 2: Lipid') +
  theme_minimal()

# Plot histogram of singular values
prod.sing.vals <- svd(save.proposed.out$test$prod)$d
p3 <- ggplot(data.frame(sing.vals = prod.sing.vals), aes(x = sing.vals)) +
  geom_histogram(binwidth = 0.1, aes(y = after_stat(count / sum(count)))) +
  ggtitle('Spectrum of Product of Projections') +
  xlab('Singular values') +
  ylab('Frequency') + 
  geom_vline(xintercept = out$test$lam, linetype = 'dashed', color='red') +
  theme_minimal()

# facet grid of 6 plots
grid.arrange(grobs = list(p3, p1, p2), ncol = 3)

