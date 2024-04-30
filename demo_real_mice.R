my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)
source('src/models_2_views.R')
library(Ckmeans.1d.dp)
library(mixOmics)

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

out <- proposed_subsampling_func(Y1, Y2, rank1, rank2, return_scores = TRUE)
out.ajive <- ajive_func(Y1, Y2, rank1, rank2, return_scores = TRUE)

data <- data.frame(diet = nutrimouse$diet, genotype = nutrimouse$genotype)
data[['Joint 1']] <- out$joint[, 1]
data[['Joint 2']] <- out$joint[, 2]
data[['Indiv 1']] <- out$indiv2[, 1]
data[['Indiv 2']] <- out$indiv2[, 2]
data[['Joint 1 AJIVE']] <- out.ajive$joint[, 1]
data[['Joint 2 AJIVE']] <- out.ajive$joint[, 2]
data[['Indiv 1 AJIVE']] <- out.ajive$indiv2[, 1]
data[['Indiv 2 AJIVE']] <- out.ajive$indiv2[, 2]

# Plot
library(ggplot2)
ggplot(data, aes(x = `Joint 1`, y = `Joint 2`, color = diet, shape = genotype)) +
  geom_point(size=2) +
  theme_minimal()

ggplot(data, aes(x = `Joint 1 AJIVE`, y = `Joint 2 AJIVE`, color = diet, shape = genotype)) +
  geom_point(size=2) +
  theme_minimal()

ggplot(data, aes(x = `Indiv 1`, y = `Indiv 2`, color = diet, shape = genotype)) +
  geom_point(size=2) +
  theme_minimal()

ggplot(data, aes(x = `Indiv 1 AJIVE`, y = `Indiv 2 AJIVE`, color = diet, shape = genotype)) +
  geom_point(size=2) +
  theme_minimal()

