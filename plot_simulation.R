# load ggplot 
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)

# load results
load("results.RData")

snrs <- unique(results$SNR_X1)
# get FD and TP to be between [0, 1]
results$MeanJointFD <- results$MeanJointFD / (min(c(results$ncol1, results$ncol2)) - results$RankJoint)
results$MeanJointTP <- results$MeanJointTP / results$RankJoint
results$MeanIndiv1FD <- results$MeanIndiv1FD / (min(c(results$ncol1, results$ncol2)) - results$RankIndiv1)
results$MeanIndiv1TP <- results$MeanIndiv1TP / results$RankIndiv1
results$MeanIndiv2FD <- results$MeanIndiv2FD / (min(c(results$ncol1, results$ncol2)) - results$RankIndiv2)
results$MeanIndiv2TP <- results$MeanIndiv2TP / results$RankIndiv2
# concatenate method and alpha if alpha is not NA
# results$Method <- ifelse(is.na(results$alpha), results$Method, paste(results$Method, results$alpha, sep="_"))
for (snr in snrs) {
    n <- results[results$SNR_X1 == snr, "nrow"][1]
    ncol1 <- results[results$SNR_X1 == snr, "ncol1"][1]
    ncol2 <- results[results$SNR_X1 == snr, "ncol2"][1]
    rj <- results[results$SNR_X1 == snr, "RankJoint"][1]
    ri1 <- results[results$SNR_X1 == snr, "RankIndiv1"][1]
    ri2 <- results[results$SNR_X1 == snr, "RankIndiv2"][1]
    plot1 <- ggplot(data = results %>% filter(SNR_X1 == snr), 
            aes(x = MeanJointFD, y = MeanJointTP, color=Method)) +
        geom_point() +
        geom_line() +
        xlab("FD") +
        ylab("TP")
    plot1 <- plot1 + geom_text(data = results %>% filter(SNR_X1 == snr, 
                                                         Method == "fd_control_joint" | 
                                                           Method == "fd_control_joint_mult"),
                               aes(x = MeanJointFD, y = MeanJointTP, label=alpha), 
                               hjust=0, vjust=0)
    plot2 <- ggplot(data = results %>% filter(SNR_X1 == snr), 
            aes(x = MeanIndiv1FD, y = MeanIndiv1TP, color=Method)) +
        geom_point() +
        geom_line() +
        xlab("FD") +
        ylab("TP")
    plot2 <- plot2 + geom_text(data = results %>% filter(SNR_X1 == snr, 
                                                           Method == "fd_control_joint" | 
                                                             Method == "fd_control_joint_mult"),
                               aes(x = MeanIndiv1FD, y = MeanIndiv1TP, label=alpha), 
                               hjust=0, vjust=0)
    plot3 <- ggplot(data = results %>% filter(SNR_X1 == snr), 
            aes(x = MeanIndiv2FD, y = MeanIndiv2TP, color=Method)) +
        geom_point() +
        geom_line() +
        xlab("FD") +
        ylab("TP")
    plot3 <- plot3 + geom_text(data = results %>% filter(SNR_X1 == snr, 
                                                         Method == "fd_control_joint" | 
                                                           Method == "fd_control_joint_mult"),
                               aes(x = MeanIndiv2FD, y = MeanIndiv2TP, label=alpha), 
                               hjust=0, vjust=0)
    # save plot on a grid
    g <- grid.arrange(plot1, plot2, plot3, ncol=3, 
                      top = textGrob(paste(c("SNR =", round(snr, 3), 
                                             "nrow =", n, "ncol1 =", ncol1, 
                                             "ncol2 =", ncol2, "rj =", rj,
                                             "ri1 =", ri1, "ri2 =", ri2), collapse=" "),
                                     gp=gpar(fontsize=14)))
    ggsave(paste("plots/snr", round(snr, 1), ".png", sep=""), g,
           width = 18, height = 3, dpi = 300,)
}