rm(list = ls())
library(rstudioapi)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(cowplot)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

labCCA <- read.csv("lab tau - CCA.csv", header = TRUE)
labCCA3 <- read.csv("lab tau - CCA3.csv", header = TRUE)
mCCA <- read.csv("m tau - CCA.csv", header = TRUE)
mCCA3 <- read.csv("m tau - CCA3.csv", header = TRUE)
tCCA <- read.csv("t tau - CCA.csv", header = TRUE)
tCCA3 <- read.csv("t tau - CCA3.csv", header = TRUE)

colnames(labCCA) <- c("Batch_1", "Batch_2", "Batch_3")
colnames(labCCA3) <- c("Batch_1", "Batch_2", "Batch_3")
colnames(mCCA) <- c("Batch_1", "Batch_2", "Batch_3")
colnames(mCCA3) <- c("Batch_1", "Batch_2", "Batch_3")
colnames(tCCA) <- c("Batch_2", "Batch_1", "Batch_3")
colnames(tCCA3) <- c("Batch_2", "Batch_1", "Batch_3")

# Convert the data frames to long format
convert_to_long <- function(df, dataset_name) {
  df %>%
    pivot_longer(cols = everything(), names_to = "Batch", values_to = "Value") %>%
    mutate(Dataset = dataset_name)
}

# Reformat data frames
labCCA_long <- convert_to_long(labCCA, "Lab CCA")
labCCA3_long <- convert_to_long(labCCA3, "Lab CCA3")
mCCA_long <- convert_to_long(mCCA, "Mayuge CCA")
mCCA3_long <- convert_to_long(mCCA3, "Mayuge CCA3")
tCCA_long <- convert_to_long(tCCA, "Tororo CCA")
tCCA3_long <- convert_to_long(tCCA3, "Tororo CCA3")

# Plot function
create_density_plot1 <- function(data, title) {
  ggplot(data, aes(x = sqrt(1 / Value), fill = Batch, y = ..scaled..)) +
    geom_density(alpha = 0.5) +
    labs(x = "Standard Deviation", y = "Density", title = title) +
    scale_fill_manual(
      values = c("#2b6f04", "#f1f1f1", "#0ade4c"),
      name = "CCA Batches"
    ) +
    theme_bw() +
    xlim(0, 3.5) +
    theme(legend.position = "right")
}

create_density_plot2 <- function(data, title) {
  ggplot(data, aes(x = sqrt(1 / Value), fill = Batch, y = ..scaled..)) +
    geom_density(alpha = 0.5) +
    labs(x = "Standard Deviation", y = "Density", title = title) +
    scale_fill_manual(
      values = c("#003f5c", "#bc5090", "#ffa600"),
      name = "CCA3 Batches"
    ) +
    theme_bw() +
    xlim(0, 3.5) +
    theme(legend.position = "right")
}
# Create the plots
plot_labCCA <- create_density_plot1(labCCA_long, "POC-CCA - Lab")
plot_labCCA3 <- create_density_plot2(labCCA3_long, "POC-CCA3 - Lab")
plot_mCCA <- create_density_plot1(mCCA_long, "POC-CCA - Mayuge")
plot_mCCA3 <- create_density_plot2(mCCA3_long, "POC-CCA3 - Mayuge")
plot_tCCA <- create_density_plot1(tCCA_long, "POC-CCA - Tororo")
plot_tCCA3 <- create_density_plot2(tCCA3_long, "POC-CCA3 - Tororo")

# Extract legends
legend_CCA <- get_legend(plot_labCCA)
legend_CCA3 <- get_legend(plot_labCCA3)

combined_grid <- plot_grid(
  plot_grid(
    plot_labCCA + theme(legend.position = "none"),
    plot_labCCA3 + theme(legend.position = "none"),
    legend_CCA,
    ncol = 3,
    rel_widths = c(1, 1, 0.3)
  ),
  plot_grid(
    plot_tCCA + theme(legend.position = "none"),
    plot_tCCA3 + theme(legend.position = "none"),
    legend_CCA3,
    ncol = 3,
    rel_widths = c(1, 1, 0.3)
  ),
  plot_grid(
    plot_mCCA + theme(legend.position = "none"),
    plot_mCCA3 + theme(legend.position = "none"),
    NULL,
    ncol = 3,
    rel_widths = c(1, 1, 0.3)
  ),
  nrow = 3
)

# Display the grid
png("Tau All.png", width = 1000, height = 700)
print(combined_grid)
dev.off()

minterCCA <- sqrt(1/read.csv("m inter tau - CCA.csv", header = TRUE))
minterCCA3 <- sqrt(1/read.csv("m inter tau - CCA3.csv", header = TRUE))
tinterCCA <- sqrt(1/read.csv("t inter tau - CCA.csv", header = TRUE))
tinterCCA3 <- sqrt(1/read.csv("t inter tau - CCA3.csv", header = TRUE))

# Create density plots for each dataset using the first column (adjust if needed)
plot_minterCCA <- ggplot(minterCCA, aes(x = minterCCA[[1]])) +
  geom_density(fill = "grey25", colour = "grey25", alpha = 0.5) +
  labs(title = "CCA Inter-Sample Variation - Mayuge", x = "Standard Deviation", y = "Density") +
  xlim(1, 2) +
  theme_minimal() +
  theme(legend.position = "none")

plot_minterCCA3 <- ggplot(minterCCA3, aes(x = minterCCA3[[1]])) +
  geom_density(fill = "grey25", colour = "grey25", alpha = 0.5) +
  labs(title = "CCA3 Inter-Sample Variation - Mayuge", x = "Standard Deviation", y = "Density") +
  xlim(1, 2) +
  theme_minimal() +
  theme(legend.position = "none")

plot_tinterCCA <- ggplot(tinterCCA, aes(x = tinterCCA[[1]])) +
  geom_density(fill = "grey25", colour = "grey25", alpha = 0.5) +
  labs(title = "CCA Inter-Sample Variation - Tororo", x = "Standard Deviation", y = "Density") +
  xlim(1, 2) +
  theme_minimal() +
  theme(legend.position = "none")

plot_tinterCCA3 <- ggplot(tinterCCA3, aes(x = tinterCCA3[[1]])) +
  geom_density(fill = "grey25", colour = "grey25", alpha = 0.5) +
  labs(title = "CCA3 Inter-Sample Variation - Tororo", x = "Standard Deviation", y = "Density") +
  xlim(1, 2) +
  theme_minimal() +
  theme(legend.position = "none")

# Combine the four plots into a 2x2 grid
combined_grid <- plot_grid(
  plot_tinterCCA, plot_tinterCCA3,
  plot_minterCCA, plot_minterCCA3,
  ncol = 2
)

# Display the combined grid
combined_grid

png("Tau inter All.png", width = 1000, height = 500)
print(combined_grid)
dev.off()
