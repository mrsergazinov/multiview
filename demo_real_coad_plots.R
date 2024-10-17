my_lib_path <- "./multiview_rlibs"
.libPaths(my_lib_path)

library(Ckmeans.1d.dp)
library(MASS)
library(RMTstat)
library(denoiseR)

library(ggplot2)
library(gridExtra)

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
data = list(Y1 = rnaData, 
            Y2 = mRNAData)
type = COAD$subtype[ID1]

########################################################## 
#----------- Plotting -----------------------------------#
##########################################################
sing.vals.rna <- svd(scale(rnaData))$d
sing.vals.mRNA <- svd(scale(mRNAData))$d
df.sing.vals <- data.frame(x = c(1:length(sing.vals.rna), 1:length(sing.vals.mRNA)),
                           sing.vals = c(sing.vals.rna, sing.vals.mRNA), 
                           View = c(rep("RNA", length(sing.vals.rna)), 
                                    rep("mRNA", length(sing.vals.mRNA)))
)
plt1 <- ggplot(df.sing.vals, aes(x = x, y = sing.vals, color = View)) +
  geom_point() +
  geom_line() +
  xlab("Singular Value Index") +
  ylab("Singular Value") +
  ggtitle("Scree plot") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
plt1 <- plt1 + geom_point(data = df.sing.vals[df.sing.vals$x == 16, ], 
                        aes(x = x, y = sing.vals), color = "red", size = 3)
plt1 <- plt1 + geom_point(data = df.sing.vals[df.sing.vals$x == 42, ],
                        aes(x = x, y = sing.vals), color = "blue", size = 3)


save.out <- proposed_func(data$Y1, data$Y2, 42, 42, return_scores = TRUE)
plt3 <- ggplot(data.frame(sing.vals = save.out$test$svd.prod$d), aes(x = sing.vals)) +
  geom_histogram(closed = "right") +
  xlab("Singular Value") +
  ylab("Frequency") +
  ggtitle("Diagnostic plot, estimated ranks") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
max_height <- max(ggplot_build(plt)$data[[1]]$y)
plt3 <- plt3 + annotate("rect",
                        xmin = 0, xmax = save.out$test$lam,
                        ymin = 0, ymax = max_height,
                        alpha=0.3, fill='blue')
plt3 <- plt3 + annotate("rect",
                        xmin = 1-mean(save.out$bootstrap$epsilon1), xmax = 1,
                        ymin = 0, ymax = max_height,
                        alpha=0.3, fill='green')

save.out <- proposed_func(data$Y1, data$Y2, 16, 16, return_scores = TRUE)
plt2 <- ggplot(data.frame(sing.vals = save.out$test$svd.prod$d), aes(x = sing.vals)) +
  geom_histogram(closed = "right") +
  xlab("Singular Value") +
  ylab("Frequency") +
  ggtitle("Diagnostic plot, manual ranks") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
max_height <- max(ggplot_build(plt)$data[[1]]$y)
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
