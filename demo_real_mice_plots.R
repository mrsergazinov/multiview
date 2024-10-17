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

# process lipids
Y2 <- Y2 * 100
Y2[Y2 == 0] <- 0.375 / (100 + 0.75)

# Standardize columns
Y1 <- scale(Y1)
Y2 <- scale(Y2)
data <- list(Y1 = Y1, Y2 = Y2)

########################################################## 
#----------- Plotting -----------------------------------#
##########################################################
sing.vals1 <- svd(scale(Y1))$d
sing.vals2 <- svd(scale(Y2))$d
rank1 <- optishrink(Y1)$nb.eigen
rank2 <- optishrink(Y2)$nb.eigen
df <- data.frame(x = 1:length(sing.vals1), 
                 sing.vals1 = sing.vals1)
plt1 <- ggplot(df, aes(x = x, y = sing.vals1, color="Gene")) +
  geom_point() +
  geom_point(data = df[df$x == rank1, ], aes(x = x, y = sing.vals1), color = "blue", size = 3) +
  geom_point(data = df[df$x == 3, ], aes(x = x, y = sing.vals1), color = "red", size = 3) +
  geom_line() +
  xlab("Singular Value Index") +
  ylab("Singular Value") +
  ggtitle("Scree plot") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
df <- data.frame(x = 1:length(sing.vals2), 
                 sing.vals2 = sing.vals2)
plt1 <- plt1 + geom_point(data = df, aes(x = x, y = sing.vals2, color = "Lipid")) + 
  geom_line(data=df, aes(x = x, y = sing.vals2, color = "Lipid")) +
  geom_point(data = df[df$x == rank2, ], aes(x = x, y = sing.vals2), color = "blue", size = 3) +
  geom_point(data = df[df$x == 4, ], aes(x = x, y = sing.vals2), color = "red", size = 3)


save.out <- proposed_func(data$Y1, data$Y2, rank1, rank2, return_scores = TRUE)
plt3 <- ggplot(data.frame(sing.vals = save.out$test$svd.prod$d), aes(x = sing.vals)) +
  geom_histogram(closed = "right") +
  xlab("Singular Value") +
  ylab("Frequency") +
  ggtitle("Diagnostic plot, estimated ranks") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
max_height <- max(ggplot_build(plt3)$data[[1]]$y)
plt3 <- plt3 + annotate("rect",
                        xmin = 0, xmax = save.out$test$lam,
                        ymin = 0, ymax = max_height,
                        alpha=0.3, fill='blue')
plt3 <- plt3 + annotate("rect",
                        xmin = 1-mean(save.out$bootstrap$epsilon1), xmax = 1,
                        ymin = 0, ymax = max_height,
                        alpha=0.3, fill='green')

save.out <- proposed_func(data$Y1, data$Y2, 3, 4, return_scores = TRUE)
plt2 <- ggplot(data.frame(sing.vals = save.out$test$svd.prod$d), aes(x = sing.vals)) +
  geom_histogram(closed = "right") +
  xlab("Singular Value") +
  ylab("Frequency") +
  ggtitle("Diagnostic plot, manual ranks") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
max_height <- max(ggplot_build(plt2)$data[[1]]$y)
plt2 <- plt2 + annotate("rect",
                        xmin = 0, xmax = save.out$test$lam,
                        ymin = 0, ymax = max_height,
                        alpha=0.3, fill='blue')
plt2 <- plt2 + annotate("rect",
                        xmin = 1-mean(save.out$bootstrap$epsilon1), xmax = 1,
                        ymin = 0, ymax = max_height,
                        alpha=0.3, fill='green')

# plot all together
grid.arrange(plt1, plt3, plt2, ncol = 3)


