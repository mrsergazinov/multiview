my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(Ckmeans.1d.dp)
library(mixOmics)
library(ajive)
library(ggplot2)

source('src/models_2_views.R')


data(nutrimouse)
Y1 <- nutrimouse$gene
Y2 <- nutrimouse$lipid
rank1 <- 3
rank2 <- 4

# process lipids
Y2 <- Y2 * 100
Y2[Y2 == 0] <- 0.375 / (100 + 0.75)

# Standardize columns
Y1 <- scale(Y1)
Y2 <- scale(Y2)

out <- proposed_func(Y1, Y2, rank1, rank2, return_scores = TRUE)

data <- data.frame(diet = nutrimouse$diet, genotype = nutrimouse$genotype)
data[['Joint 1']] <- out$joint[, 1]
data[['Joint 2']] <- out$joint[, 2]
data[['Individual 1']] <- out$indiv2[, 1]
data[['Individual 2']] <- out$indiv2[, 2]

# Plot
p1 <- ggplot(data, aes(x = `Joint 1`, y = `Joint 2`, color = diet, shape = genotype)) +
  geom_point(size=3) +
  xlab('Component 1') + 
  ylab('Component 2') +
  ggtitle('Joint View: Gene and Lipid') +
  theme_minimal()

p2 <- ggplot(data, aes(x = `Individual 1`, y = `Individual 2`, color = diet, shape = genotype)) +
  geom_point(size=3) +
  xlab('Component 1') + 
  ylab('Component 2') +
  ggtitle('Individual View 2: Lipid') +
  theme_minimal()

# Plot histogram of singular values
prod.sing.vals <- svd(out$test$prod)$d
p3 <- ggplot(data.frame(sing.vals = prod.sing.vals), aes(x = sing.vals)) +
  geom_histogram(binwidth = 0.1, aes(y = after_stat(count / sum(count)))) +
  ggtitle('Spectrum of Product of Projections') +
  xlab('Singular values') +
  ylab('Frequency') + 
  geom_vline(xintercept = out$test$lam, linetype = 'dashed', color='red') +
  theme_minimal()

# facet grid of 6 plots
grid.arrange(grobs = list(p3, p1, p2), ncol = 3)

