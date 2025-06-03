rm(list = ls())
library(rstudioapi)
library(runjags)
library(ggplot2)
library(rjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

dt <- read.csv("CORNTD_recPOCCCA_2023_Tororo.csv")

KK_day1 <- rowMeans(dt[, c(9, 10)], na.rm = TRUE)
KK_day2 <- rowMeans(dt[, c(19, 20)], na.rm = TRUE)
KK_day3 <- rowMeans(dt[, c(29, 30)], na.rm = TRUE)
KK_daily <- c(KK_day1, KK_day2, KK_day3)

CCA1 <- unlist(dt[,c(39,41,43)]) - 1
CCA2 <- unlist(dt[,c(45)]) - 1
CCA3 <- unlist(dt[,c(47)]) - 1
CCA <- c(CCA1, CCA2, CCA3)

REC_CCA1 <- unlist(dt[,c(49,52,55)]) - 1
REC_CCA2 <- unlist(dt[,c(59)]) - 1
REC_CCA3 <- unlist(dt[,c(63)]) - 1

REC_CCA <- c(REC_CCA1, REC_CCA2, REC_CCA3)

# Perform Spearman's rank correlation test
spearman_result1 <- cor.test(REC_CCA, CCA, method = "spearman")
print(spearman_result1)


# Compute daily averages for KK (each day has 2 replicates)
KK_day1 <- rowMeans(dt[, c(9, 10)], na.rm = TRUE)
KK_day2 <- rowMeans(dt[, c(19, 20)], na.rm = TRUE)
KK_day3 <- rowMeans(dt[, c(29, 30)], na.rm = TRUE)
KK_daily <- c(KK_day1, KK_day2, KK_day3)

# Compute daily averages for REC_CCA
# Day 1: average of triplicate measurements (subtracting 1 from each)
REC_CCA_day1 <- rowMeans(dt[, c(49, 52, 55)] - 1, na.rm = TRUE)
# Day 2 and Day 3: single measurement (subtract 1)
REC_CCA_day2 <- dt[, 59] - 1
REC_CCA_day3 <- dt[, 63] - 1
REC_CCA_daily <- c(REC_CCA_day1, REC_CCA_day2, REC_CCA_day3)

# Run Spearman's rank correlation test between the daily averages
spearman_result2 <- cor.test(REC_CCA_daily, KK_daily, method = "spearman")
print(spearman_result2)


dt <- read.csv("CORNTD_recPOCCCA_2023_Mayuge.csv")

KK <- unlist(dt[,c(8,9,18,19,28,29)])

CCA1 <- unlist(dt[,c(38,40,42)]) - 1
CCA2 <- unlist(dt[,c(44)]) - 1
CCA3 <- unlist(dt[,c(46)]) - 1
CCA <- c(CCA1, CCA2, CCA3)

REC_CCA1 <- unlist(dt[,c(48,51,54)]) - 1
REC_CCA2 <- unlist(dt[,58]) - 1
REC_CCA3 <- unlist(dt[,62]) - 1
REC_CCA <- c(REC_CCA1, REC_CCA2, REC_CCA3)

# Perform Spearman's rank correlation test
spearman_result3 <- cor.test(REC_CCA, CCA, method = "spearman")
print(spearman_result3)
# Compute daily averages for KK (each day has 2 replicates)
KK_day1 <- rowMeans(dt[, c(8,9)], na.rm = TRUE)
KK_day2 <- rowMeans(dt[, c(18,19)], na.rm = TRUE)
KK_day3 <- rowMeans(dt[, c(28,29)], na.rm = TRUE)
KK_daily <- c(KK_day1, KK_day2, KK_day3)

# Compute daily averages for REC_CCA
# Day 1: average of triplicate measurements (subtracting 1 from each)
REC_CCA_day1 <- rowMeans(dt[, c(38,40,42)] - 1, na.rm = TRUE)
# Day 2 and Day 3: single measurement (subtract 1)
REC_CCA_day2 <- dt[, 58] - 1
REC_CCA_day3 <- dt[, 62] - 1
REC_CCA_daily <- c(REC_CCA_day1, REC_CCA_day2, REC_CCA_day3)

# Run Spearman's rank correlation test between the daily averages
spearman_result4 <- cor.test(REC_CCA_daily, KK_daily, method = "spearman")
print(spearman_result4)
