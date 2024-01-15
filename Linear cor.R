library(rstudioapi)
library(ggplot2)
library(tidyr)
library(dplyr)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
Sim <- read.csv("Spiked CCA.csv", header = FALSE)
Sim <- Sim[,-5]
Sim <- cbind(Sim[1:12,], Sim[13:24,], Sim[25:36,])

df = data.frame(matrix(nrow = (nrow(Sim)), ncol = 7)) 
df[,1] <- Sim[,1]
df[,2] <- rowMeans(subset(Sim, select=2:4))
df[,3] <- rowMeans(subset(Sim, select=6:8))
df[,4] <- rowMeans(subset(Sim, select=10:12))
df[,5] <- apply(subset(Sim, select=2:4), 1, sd, na.rm=TRUE)
df[,6] <- apply(subset(Sim, select=6:8), 1, sd, na.rm=TRUE)
df[,7] <- apply(subset(Sim, select=10:12), 1, sd, na.rm=TRUE)
df <- df %>% 
  rename(
    ug_ml = X1,
    G_score_1 = X2,
    G_score_2 = X3,
    G_score_3 = X4,
    Standard_Deviation_1 = X5,
    Standard_Deviation_2 = X6,
    Standard_Deviation_3 = X7
  )

plot1 <- ggplot(data = df, aes(x = ug_ml)) + 
  geom_line(aes(y = (G_score_1), color = "Batch 1"), size = 1) + 
  geom_line(aes(y = (G_score_2), color = "Batch 2"), size = 1) + 
  geom_line(aes(y = (G_score_3), color = "Batch 3"), size = 1) + 
  geom_ribbon(aes(y = (G_score_1), ymin = (G_score_1 - Standard_Deviation_1), ymax = (G_score_1 + Standard_Deviation_1), fill = "Batch 1", colour = "Batch 1"), alpha = .2) +
  geom_ribbon(aes(y = (G_score_2), ymin = (G_score_2 - Standard_Deviation_2), ymax = (G_score_2 + Standard_Deviation_2), fill = "Batch 2", colour = "Batch 2"), alpha = .2) +
  geom_ribbon(aes(y = (G_score_3), ymin = (G_score_3 - Standard_Deviation_3), ymax = (G_score_3 + Standard_Deviation_3), fill = "Batch 3", colour = "Batch 3"), alpha = .2) +
  scale_x_log10()+
  theme_bw() +
  xlab("Concentration (ug_ml)") + 
  ylab("G Score") + 
  ggtitle("G Score/conetration for\nCCA diagnostic test")

Sim <- read.csv("Spiked REC-CCA.csv", header = FALSE)
Sim <- Sim[,-5]
Sim <- cbind(Sim[1:12,], Sim[13:24,], Sim[25:36,])

df = data.frame(matrix(nrow = (nrow(Sim)), ncol = 7)) 
df[,1] <- Sim[,1]
df[,2] <- rowMeans(subset(Sim, select=2:4))
df[,3] <- rowMeans(subset(Sim, select=6:8))
df[,4] <- rowMeans(subset(Sim, select=10:12))
df[,5] <- apply(subset(Sim, select=2:4), 1, sd, na.rm=TRUE)
df[,6] <- apply(subset(Sim, select=6:8), 1, sd, na.rm=TRUE)
df[,7] <- apply(subset(Sim, select=10:12), 1, sd, na.rm=TRUE)
df <- df %>% 
  rename(
    ug_ml = X1,
    G_score_1 = X2,
    G_score_2 = X3,
    G_score_3 = X4,
    Standard_Deviation_1 = X5,
    Standard_Deviation_2 = X6,
    Standard_Deviation_3 = X7
  )

plot2 <- ggplot(data = df, aes(x = ug_ml)) + 
  geom_line(aes(y = (G_score_1), color = "Batch 1"), size = 1) + 
  geom_line(aes(y = (G_score_2), color = "Batch 2"), size = 1) + 
  geom_line(aes(y = (G_score_3), color = "Batch 3"), size = 1) + 
  geom_ribbon(aes(y = (G_score_1), ymin = (G_score_1 - Standard_Deviation_1), ymax = (G_score_1 + Standard_Deviation_1), fill = "Batch 1", colour = "Batch 1"), alpha = .2) +
  geom_ribbon(aes(y = (G_score_2), ymin = (G_score_2 - Standard_Deviation_2), ymax = (G_score_2 + Standard_Deviation_2), fill = "Batch 2", colour = "Batch 2"), alpha = .2) +
  geom_ribbon(aes(y = (G_score_3), ymin = (G_score_3 - Standard_Deviation_3), ymax = (G_score_3 + Standard_Deviation_3), fill = "Batch 3", colour = "Batch 3"), alpha = .2) +
  scale_x_log10()+
  theme_bw() +
  xlab("Concentration (ug_ml)") + 
  ylab("G Score") + 
  ggtitle("G Score/conetration for\nrecCCA diagnostic test")


library(cowplot)
library(gridExtra)

comb_plot <- plot_grid(
  plot1 + theme(legend.position = "none",
  plot.title = element_text(size = rel(1.3)),     # Increase plot title size
  axis.title.x = element_text(size = rel(1.2)),   # Increase x axis title size
  axis.title.y = element_text(size = rel(1.2)),   # Increase y axis title size
  axis.text.x = element_text(size = rel(1.2)),    # Increase x axis text size
  axis.text.y = element_text(size = rel(1.2))),
  plot2+ theme(plot.title = element_text(size = rel(1.3)),     # Increase plot title size
               axis.title.x = element_text(size = rel(1.2)),   # Increase x axis title size
               axis.title.y = element_text(size = rel(1.2)),   # Increase y axis title size
               axis.text.x = element_text(size = rel(1.2)),    # Increase x axis text size
               axis.text.y = element_text(size = rel(1.2)),
               legend.title = element_text(size = 12),  # Adjust size for legend title
               legend.text = element_text(size = 10)),
  ncol = 2, 
  nrow = 1,
  align = "h",
  rel_widths = c(1, 1.3)  
)
png("ul.png", 700, 340)
comb_plot
dev.off()
