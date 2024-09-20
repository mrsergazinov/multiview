my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(reticulate)
library(denoiseR)
library(MASS)
library(ggplot2)
library(gridExtra)

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

out <- proposed_func(data$Y1, data$Y2, rank1, rank2, rank_joint = 4, return_scores = TRUE)  
data[[paste0("proposed", "_joint")]] <- out$joint
data[[paste0("proposed", "_indiv1")]] <- out$indiv1
data[[paste0("proposed", "_indiv2")]] <- out$indiv2

n_sim <- 100
results <- matrix(NA, nrow = n_sim, ncol = 3)
for (i in 1:n_sim){
  # split into training and test data
  set.seed(i)
  id = sample(1:167, 167)
  trainID = id[1:132] # 80% of the data 
  testID = id[133:167] # 20% of the data
  
  model_results <- c()
  for (component in c("_joint", "_indiv1", "_indiv2")){
    train_data <- data.frame(subtype = type[trainID], data[[paste0("proposed", component)]][trainID, ])
    test_data <- data.frame(subtype = type[testID], data[[paste0("proposed", component)]][testID, ])
    colnames(train_data) <- colnames(test_data)
    predictive_model <- multinom(subtype ~ ., data = train_data)
    predicted_labels <- predict(predictive_model, test_data, type = "class")
    error <- sum(predicted_labels != type[testID]) / length(type[testID])
    model_results <- c(model_results, error)
    }
  results[i, ] <- model_results
}


means <- colMeans(results, na.rm = TRUE)
sds <- apply(results, 2, sd, na.rm = TRUE)
str = "proposed"
for (i in 1:3){
    str <- paste0(str, " & ", round(means[i], 2), " (", round(sds[i], 2), ")")
}
print(str)

prod.eigenvals <- out$test$svd.prod$d^2
p <- ggplot(data.frame(prod.eigenvals = prod.eigenvals), aes(x = prod.eigenvals)) +
  geom_histogram(binwidth = 0.1, aes(y = ..density..)) +
  ggtitle('Spectrum of Product of Projections') +
  xlab('Singular values') +
  ylab('Frequency') +
  geom_vline(xintercept = out$test$lam^2, linetype = 'dashed', color='red') +
  theme_minimal()

# Define the function f(lambda)
f_lambda <- function(lambda, q1, q2) {
  # Compute lambda_+ and lambda_-
  lambda_plus <- q1 + q2 - 2 * q1 * q2 + 2 * sqrt(q1 * q2 * (1 - q1) * (1 - q2))
  lambda_minus <- q1 + q2 - 2 * q1 * q2 - 2 * sqrt(q1 * q2 * (1 - q1) * (1 - q2))
  
  # Calculate f(lambda)
  numerator <- sqrt((lambda_plus - lambda) * (lambda - lambda_minus))
  denominator <- 2 * pi * lambda * (1 - lambda)
  
  # Ensure we handle cases where denominator might be zero
  result <- ifelse(denominator != 0, numerator / denominator, NA)
  
  return(result)
}

# Set values for q1 and q2
m <- dim(data$Y1)[1]
q1 <- rank1 / m  # Example value
q2 <- rank2 / m  # Example value

# Generate a sequence of lambda values between 0 and 1 (excluding 0 and 1 to avoid division by zero)
lambda_values <- seq(0.01, 0.99, by = 0.01)

# Compute f(lambda) for each lambda value
f_values <- sapply(lambda_values, f_lambda, q1 = q1, q2 = q2)

# Create a data frame for plotting
data <- data.frame(lambda = lambda_values, f_lambda = f_values)

# Plot the function using ggplot2
p <- p + 
  geom_line(data=data, aes(x = lambda, y = f_lambda), color = "blue", size = 1)


sing.vals.indiv1 <- out$svd.indiv1$d
p2 <- ggplot(data.frame(sing.vals.indiv1 = sing.vals.indiv1), aes(x = sing.vals.indiv1)) +
  geom_histogram(binwidth = 0.1, aes(y = ..density..)) +
  ggtitle('Spectrum of Individual 1 Projections') +
  xlab('Singular values') +
  ylab('Frequency') +
  geom_vline(xintercept = out$lam1^2, linetype = 'dashed', color='red') +
  theme_minimal()

p3 <- ggplot(data.frame(sing.vals.indiv2 = out$svd.indiv2$d), aes(x = out$svd.indiv2$d)) +
  geom_histogram(binwidth = 0.1, aes(y = ..density..)) +
  ggtitle('Spectrum of Individual 2 Projections') +
  xlab('Singular values') +
  ylab('Frequency') +
  geom_vline(xintercept = out$lam2^2, linetype = 'dashed', color='red') +
  theme_minimal()

grid.arrange(grobs = list(p, p2, p3), ncol = 3)

