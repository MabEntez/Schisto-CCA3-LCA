rm(list = ls())
library(rstudioapi)
library(runjags)
library(ggplot2)
library(rjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

dt_t <- read.csv("CORNTD_recPOCCCA_2023_Tororo.csv")
dt_m <- read.csv("CORNTD_recPOCCCA_2023_Mayuge.csv")

dt_t <- dt_t[,c(9,10,19,20,29,30,39,41,43,45,47,49,52,55,59,63)]
dt_m <- dt_m[,c(8,9,18,19,28,29, 38,40,42,44,46 ,48,51,54,58,62)]

t_kk1 <- length(which(rowSums(dt_t[,1:2])>0))/length(na.omit(rowSums(dt_t[,1:2])>0))
t_kk3 <- length(which(rowSums(dt_t[,1:6])>0))/length(na.omit(rowSums(dt_t[,1:6])>0))

#### Prev Calc ####

t_cca1_g2 <- length(which((dt_t[, 7] >= 2))) / length(na.omit((dt_t[, 7] >= 2)))
t_cca3_g2 <- length(which((dt_t[, 12] >= 2))) / length(na.omit((dt_t[, 12] >= 2)))
t_cca1_g3 <- length(which((dt_t[, 7] >= 3))) / length(na.omit((dt_t[, 7] >= 3)))
t_cca3_g3 <- length(which((dt_t[, 12] >= 3))) / length(na.omit((dt_t[, 12] >= 3)))
t_cca1_g4 <- length(which((dt_t[, 7] >= 4))) / length(na.omit((dt_t[, 7] >= 4)))
t_cca3_g4 <- length(which((dt_t[, 12] >= 4))) / length(na.omit((dt_t[, 12] >= 4)))

t_all_cca1_g2 <- length(which(rowMeans(dt_t[, 7:11]) >= 2)) / length(na.omit(rowMeans(dt_t[, 7:11]) >= 2))
t_all_cca3_g2 <- length(which(rowMeans(dt_t[, 12:16]) >= 2)) / length(na.omit(rowMeans(dt_t[, 12:16]) >= 2))
t_all_cca1_g2.5 <- length(which(rowMeans(dt_t[, 7:11]) >= 2.5)) / length(na.omit(rowMeans(dt_t[, 7:11]) >= 2.5))
t_all_cca3_g2.5 <- length(which(rowMeans(dt_t[, 12:16]) >= 2.5)) / length(na.omit(rowMeans(dt_t[, 12:16]) >= 2.5))
t_all_cca1_g3 <- length(which(rowMeans(dt_t[, 7:11]) >= 3)) / length(na.omit(rowMeans(dt_t[, 7:11]) >= 3))
t_all_cca3_g3 <- length(which(rowMeans(dt_t[, 12:16]) >= 3)) / length(na.omit(rowMeans(dt_t[, 12:16]) >= 3))
t_all_cca1_g4 <- length(which(rowMeans(dt_t[, 7:11]) >= 4)) / length(na.omit(rowMeans(dt_t[, 7:11]) >= 4))
t_all_cca3_g4 <- length(which(rowMeans(dt_t[, 12:16]) >= 4)) / length(na.omit(rowMeans(dt_t[, 12:16]) >= 4))

m_kk1 <- length(which(rowSums(dt_m[,1:2])>0))/length(na.omit(rowSums(dt_m[,1:2])>0))
m_kk3 <- length(which(rowSums(dt_m[,1:6])>0))/length(na.omit(rowSums(dt_m[,1:6])>0))

m_cca1_g2 <- length(which((dt_m[, 7] >= 2))) / length(na.omit((dt_m[, 7] >= 2)))
m_cca3_g2 <- length(which((dt_m[, 12] >= 2))) / length(na.omit((dt_m[, 12] >= 2)))
m_cca1_g3 <- length(which((dt_m[, 7] >= 3))) / length(na.omit((dt_m[, 7] >= 3)))
m_cca3_g3 <- length(which((dt_m[, 12] >= 3))) / length(na.omit((dt_m[, 12] >= 3)))
m_cca1_g4 <- length(which((dt_m[, 7] >= 4))) / length(na.omit((dt_m[, 7] >= 4)))
m_cca3_g4 <- length(which((dt_m[, 12] >= 4))) / length(na.omit((dt_m[, 12] >= 4)))

m_all_cca3_g2 <- length(which(rowMeans(dt_m[, 12:16]) >= 2)) / length(na.omit(rowMeans(dt_m[, 12:16]) >= 2))
m_all_cca1_g2 <- length(which(rowMeans(dt_m[, 7:11]) >= 2)) / length(na.omit(rowMeans(dt_m[, 7:11]) >= 2))
m_all_cca3_g2.5 <- length(which(rowMeans(dt_m[, 12:16]) >= 2.5)) / length(na.omit(rowMeans(dt_m[, 12:16]) >= 2.5))
m_all_cca1_g2.5 <- length(which(rowMeans(dt_m[, 7:11]) >= 2.5)) / length(na.omit(rowMeans(dt_m[, 7:11]) >= 2.5))
m_all_cca1_g3 <- length(which(rowMeans(dt_m[, 7:11]) >= 3)) / length(na.omit(rowMeans(dt_m[, 7:11]) >= 3))
m_all_cca3_g3 <- length(which(rowMeans(dt_m[, 12:16]) >= 3)) / length(na.omit(rowMeans(dt_m[, 12:16]) >= 3))
m_all_cca1_g4 <- length(which(rowMeans(dt_m[, 7:11]) >= 4)) / length(na.omit(rowMeans(dt_m[, 7:11]) >= 4))
m_all_cca3_g4 <- length(which(rowMeans(dt_m[, 12:16]) >= 4)) / length(na.omit(rowMeans(dt_m[, 12:16]) >= 4))

df <- data.frame(
  Category = c("kk1", "kk3", "cca1_g2", "cca3_g2", "cca1_g3", "cca3_g3", "cca1_g4", "cca3_g4", "cca1_all_g2", "cca3_all_g2","cca1_all_g2.5", "cca3_all_g2.5", "cca1_all_g3", "cca3_all_g3", "cca1_all_g4", "cca3_all_g4"),
  T_Values = c(t_kk1, t_kk3, t_cca1_g2, t_cca3_g2, t_cca1_g3, t_cca3_g3, t_cca1_g4, t_cca3_g4, t_all_cca1_g2, t_all_cca3_g2, t_all_cca1_g2.5, t_all_cca3_g2.5, t_all_cca1_g3, t_all_cca3_g3, t_all_cca1_g4, t_all_cca3_g4),
  M_Values = c(m_kk1, m_kk3, m_cca1_g2, m_cca3_g2, m_cca1_g3, m_cca3_g3, m_cca1_g4, m_cca3_g4, m_all_cca1_g2, m_all_cca3_g2, m_all_cca1_g2.5, m_all_cca3_g2.5, m_all_cca1_g3, m_all_cca3_g3, m_all_cca1_g4, m_all_cca3_g4)
)

# Save the data frame to a CSV file
write.csv(df, "sampled prev.csv", row.names = FALSE)

mean(unlist(dt_t[7:11]), na.rm = T)



#### Misc Calc ####
