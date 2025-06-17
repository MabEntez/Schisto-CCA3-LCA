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

dt_t[,c(1:6)] <- dt_t[,c(1:6)] * 24
dt_m[,c(1:6)] <- dt_m[,c(1:6)] * 24

mean(unlist(dt_t), na.rm = T)
mean(unlist(dt_m), na.rm = T)

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



t_kk1 <- length(which(rowSums(dt_t[,1:2])>0))/length(na.omit(rowSums(dt_t[,1:2])>0))
t_kk3 <- length(which(rowSums(dt_t[,1:6])>0))/length(na.omit(rowSums(dt_t[,1:6])>0))

#### Misc Calc ####
library(tidyverse)
#### G-score KK mean ####
# 1. Compute KK_mean and cut into a factor with all four levels
dt_t_a <- dt_t %>%
  mutate(
    KK_mean = rowMeans(select(., 1:6), na.rm = TRUE),
    intensity = cut(
      KK_mean,
      breaks = c(-Inf, 0, 100, 400, Inf),
      labels = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram"),
      right = TRUE,
      include.lowest = TRUE
    ),
    # force the factor to keep those four levels in order
    intensity = factor(intensity,
                       levels = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram"))
  ) %>%
  filter(!is.na(intensity))

# 2. Compute mean G-scores
dt_t_a <- dt_t_a %>%
  mutate(
    mean_CCA  = rowMeans(select(., 7:11),  na.rm = TRUE),
    mean_CCA3 = rowMeans(select(., 12:16), na.rm = TRUE)
  )

# 3. Pivot to long form
dt_t_a_long <- dt_t_a %>%
  select(intensity, mean_CCA, mean_CCA3) %>%
  pivot_longer(
    cols      = starts_with("mean_"),
    names_to  = "assay",
    values_to = "mean_Gscore"
  ) %>%
  mutate(
    assay = recode(assay,
                   mean_CCA  = "POC-CCA",
                   mean_CCA3 = "POC-CCA3")
  )

# 4. Plot, telling ggplot not to drop empty categories
t_avgcomp_plot <- ggplot(dt_t_a_long, aes(x = intensity, y = mean_Gscore, fill = intensity)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~ assay, ncol = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    values = grey.colors(4, start = 1, end = 0, rev = FALSE),
    # ensure the names line up with your factor levels:
    breaks = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram")
  ) +
  scale_y_continuous(
    breaks = seq(2, 10, 2),
    limits = c(1, 10)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x        = "Kato-Katz eggs per Gram",
    y        = "Mean G-score",
    title    = "Tororo - Mean POC-CCA and POC-CCA3 G-scores\nby Kato-Katz eggs per gram"
  )

# 1. Compute KK_mean and cut into a factor with all four levels
dt_m_a <- dt_m %>%
  mutate(
    KK_mean = rowMeans(select(., 1:6), na.rm = TRUE),
    intensity = cut(
      KK_mean,
      breaks = c(-Inf, 0, 100, 400, Inf),
      labels = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram"),
      right = TRUE,
      include.lowest = TRUE
    ),
    # force the factor to keep those four levels in order
    intensity = factor(intensity,
                       levels = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram"))
  ) %>%
  filter(!is.na(intensity))

# 2. Compute mean G-scores
dt_m_a <- dt_m_a %>%
  mutate(
    mean_CCA  = rowMeans(select(., 7:11),  na.rm = TRUE),
    mean_CCA3 = rowMeans(select(., 12:16), na.rm = TRUE)
  )

# 3. Pivot to long form
dt_m_a_long <- dt_m_a %>%
  select(intensity, mean_CCA, mean_CCA3) %>%
  pivot_longer(
    cols      = starts_with("mean_"),
    names_to  = "assay",
    values_to = "mean_Gscore"
  ) %>%
  mutate(
    assay = recode(assay,
                   mean_CCA  = "POC-CCA",
                   mean_CCA3 = "POC-CCA3")
  )

# 4. Plot, telling ggplot not to drop empty categories
m_avgcomp_plot <- ggplot(dt_m_a_long, aes(x = intensity, y = mean_Gscore, fill = intensity)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~ assay, ncol = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    values = grey.colors(4, start = 1, end = 0, rev = FALSE),
    # ensure the names line up with your factor levels:
    breaks = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram")
  ) +
  scale_y_continuous(
    breaks = seq(0, 10, 2),
    limits = c(1, 10)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x        = "Kato-Katz eggs per Gram",
    y        = "Mean G-score",
    title    = "Mayuge - Mean POC-CCA and POC-CCA3 G-scores\nby Kato-Katz eggs per gram"
  )

#### G-score KK single ####

# 1. Compute KK_mean and cut into a factor with all four levels
dt_t_s <- dt_t %>%
  mutate(
    KK_mean = rowMeans(select(., 1), na.rm = TRUE),
    intensity = cut(
      KK_mean,
      breaks = c(-Inf, 0, 100, 400, Inf),
      labels = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram"),
      right = TRUE,
      include.lowest = TRUE
    ),
    # force the factor to keep those four levels in order
    intensity = factor(intensity,
                       levels = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram"))
  ) %>%
  filter(!is.na(intensity))

# 2. Compute mean G-scores
dt_t_s <- dt_t_s %>%
  mutate(
    mean_CCA  = rowMeans(select(., 7),  na.rm = TRUE),
    mean_CCA3 = rowMeans(select(., 12), na.rm = TRUE)
  )

# 3. Pivot to long form
dt_t_s_long <- dt_t_s %>%
  select(intensity, mean_CCA, mean_CCA3) %>%
  pivot_longer(
    cols      = starts_with("mean_"),
    names_to  = "assay",
    values_to = "mean_Gscore"
  ) %>%
  mutate(
    assay = recode(assay,
                   mean_CCA  = "POC-CCA",
                   mean_CCA3 = "POC-CCA3")
  )

# 4. Plot, telling ggplot not to drop empty categories
t_singcomp_plot <- ggplot(dt_t_s_long, aes(x = intensity, y = mean_Gscore, fill = intensity)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~ assay, ncol = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    values = grey.colors(4, start = 1, end = 0, rev = FALSE),
    # ensure the names line up with your factor levels:
    breaks = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram")
  ) +
  scale_y_continuous(
    breaks = seq(0, 10, 2),
    limits = c(1, 10)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x        = "Kato-Katz eggs per Gram",
    y        = "Mean G-score",
    title    = "Tororo - Single POC-CCA and POC-CCA3 G-scores\nby Kato-Katz eggs per gram"
  )

# 1. Compute KK_mean and cut into a factor with all four levels
dt_m_s <- dt_m %>%
  mutate(
    KK_mean = rowMeans(select(., 1), na.rm = TRUE),
    intensity = cut(
      KK_mean,
      breaks = c(-Inf, 0, 100, 400, Inf),
      labels = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram"),
      right = TRUE,
      include.lowest = TRUE
    ),
    # force the factor to keep those four levels in order
    intensity = factor(intensity,
                       levels = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram"))
  ) %>%
  filter(!is.na(intensity))

# 2. Compute mean G-scores
dt_m_s <- dt_m_s %>%
  mutate(
    mean_CCA  = rowMeans(select(., 7),  na.rm = TRUE),
    mean_CCA3 = rowMeans(select(., 12), na.rm = TRUE)
  )

# 3. Pivot to long form
dt_m_s_long <- dt_m_s %>%
  select(intensity, mean_CCA, mean_CCA3) %>%
  pivot_longer(
    cols      = starts_with("mean_"),
    names_to  = "assay",
    values_to = "mean_Gscore"
  ) %>%
  mutate(
    assay = recode(assay,
                   mean_CCA  = "POC-CCA",
                   mean_CCA3 = "POC-CCA3")
  )

# 4. Plot, telling ggplot not to drop empty categories
m_singcomp_plot <- ggplot(dt_m_s_long, aes(x = intensity, y = mean_Gscore, fill = intensity)) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~ assay, ncol = 1) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    values = grey.colors(4, start = 1, end = 0, rev = FALSE),
    # ensure the names line up with your factor levels:
    breaks = c("Negative", "1–100 Eggs per gram", "101–400 Eggs per gram", ">400 Eggs per gram")
  ) +
  scale_y_continuous(
    breaks = seq(0, 10, 2),
    limits = c(1, 10)
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x        = "Kato-Katz eggs per Gram",
    y        = "Mean G-score",
    title    = "Mayuge - Single POC-CCA and POC-CCA3 G-scores\nby Kato-Katz eggs per gram"
  )

t_avgcomp_plot
t_singcomp_plot
m_avgcomp_plot
m_singcomp_plot

# install.packages("gridExtra")  # if you don’t already have it
library(gridExtra)

# Draw to screen:
grid.arrange(
  t_avgcomp_plot, t_singcomp_plot,
  m_avgcomp_plot, m_singcomp_plot,
  ncol = 2
)

# If you also want to save as a single 800×500 PNG:
png("all_plots_grid.png", width = 1200, height = 900)
grid.arrange(
  t_avgcomp_plot, m_avgcomp_plot,
  ncol = 1
)
dev.off()


