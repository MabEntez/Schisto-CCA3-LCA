rm(list = ls())
library(rstudioapi)
library(dplyr)
library(tidyverse)
library(runjags)
library(ggplot2)
library(rjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

#### Jags Model ####
dt <- read.csv("Spiked CCA.csv")
dt <- dt[, 1:4]
colnames(dt) <- c("Dosage", "Batch_1", "Batch_2", "Batch_3")

CCA1 <- dt[, 2] -1
CCA2 <- dt[, 3] -1
CCA3 <- dt[, 4] -1

Spiked <- as.vector(dt[,1])

model_str = "model {
  for (i in 1:n){
    ## CCA component
    x1[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(log(y[i]+ 1e-6)-multiParam[2]))),tau1)
    x2[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(log(y[i]+ 1e-6)-multiParam[2]))),tau2)
    x3[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(log(y[i]+ 1e-6)-multiParam[2]))),tau3)
  }
  tau1 ~ dgamma(0.001,0.001)
  tau2 ~ dgamma(0.001,0.001)
  tau3 ~ dgamma(0.001,0.001)
  multiParam[1] ~ dnorm(0, 0.0001)  # non-informative prior for multiParam[1]
  multiParam[2] ~ dnorm(0, 0.0001)  # non-informative prior for multiParam[2]
}"

data_jags = list(
  n = length(CCA1),
  x1 = CCA1,
  x2 = CCA2,
  x3 = CCA3,
  y = Spiked
)

method <- list(
  "burnin" = 3000,
  "sample" = 10000
)

# Run the model
results2 <- run.jags(model_str,
                    data = data_jags,
                    monitor = c("tau1", "tau2", "tau3", "multiParam[1]", "multiParam[2]"),
                    n.chains = 4,
                    burnin = method$burnin,
                    sample = method$sample
)

# Summary statistics
summary(results2)

# Plotting
plot(results2)

mcmc_samples <- as.mcmc(results2)


#### tau plotting ####
tau_samples <- as.mcmc(results2, variable.names = "tau")[,1:3]

tau_df2 <- data.frame(tau_samples)

t2 <- ggplot() +
  geom_density(data=tau_df2, aes(x= 1/ sqrt(tau1)), fill="#003f5c", alpha=0.5) +
  geom_density(data=tau_df2, aes(x= 1/ sqrt(tau2)), fill="#bc5090", alpha=0.5) +
  geom_density(data=tau_df2, aes(x= 1/ sqrt(tau3)), fill="#ffa600", alpha=0.5) +
  labs(x="Standard deviation", y="Density", title="POC-CCA")+
  theme_bw()

library(cowplot)
library(gridExtra)
library(scales)

combined_plot <- plot_grid(
  t2+ theme(plot.title = element_text(size = rel(1.3)),
            axis.title.x = element_text(size = rel(1.2)),
            axis.title.y = element_text(size = rel(1.2)),
            axis.text.x = element_text(size = rel(1.3)),
            axis.text.y = element_text(size = rel(1.3)),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)) + 
    ylim(0, 6.5) +
    xlim(0, 2),
  t1+ theme(plot.title = element_text(size = rel(1.3)),
            axis.title.x = element_text(size = rel(1.2)),
            axis.title.y = element_text(size = rel(1.2)),
            axis.text.x = element_text(size = rel(1.3)),
            axis.text.y = element_text(size = rel(1.3)),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)) + 
    ylim(0, 6.5) +
    xlim(0, 2),
  ncol = 2, 
  nrow = 1,
  align = "h"
)
png("tau.png", 800, 350)
combined_plot
dev.off()

svg("tau.svg", width = 8, height = 3.5)
print(combined_plot)
dev.off()
####Fitting gamma distribution####
library(fitdistrplus)
fit.gamma <- fitdist(tau_df2[,2], distr = "gamma", method = "mle")
plot(fit.gamma)
summary(fit.gamma)

#### Getting the Curve ####

sampled_particles <- sample(1:nrow(mcmc_samples), 1000)
MP <- mcmc_samples[ sampled_particles, c("multiParam[1]", "multiParam[2]")]
Spiked <- dt[1:12,1]

xaxis <- seq(0, 10000, by = 10)
G_score <- matrix(NA,nrow = length(xaxis), ncol = 1001)
for (i in 1:length(xaxis)){
  G_score[i,1] <- xaxis[i] 
  for (j in 1:nrow(MP)){
    G_score[i, j+1] <- 9 / (1 + exp(-MP[j,1]*(log(xaxis[i])-MP[j,2])))
  }
}
G_score_long <- pivot_longer(as.data.frame(G_score), cols = -V1, values_to = "Value")

ggplot(G_score_long, aes(x = V1, y = Value, group = name)) +
  geom_line() +
  scale_x_log10()+
  theme_bw() +
  xlab("ng/ml") + 
  ylab("G Score") + 
  ggtitle("POC-CCA")

Sim <- read.csv("Spiked CCA.csv", header = FALSE)
Sim <- Sim[,-5]
colnames(Sim) <- c("ngml", "Batch 1", "Batch 2", "Batch 3")
Sim_long <- pivot_longer(as.data.frame(Sim), cols = -"ngml", values_to = "G_score")
colnames(Sim_long) <- c("ngml", "Lab_Data", "G_score")

base_plot <- ggplot() +
  scale_x_log10() +
  theme_bw() +
  xlab("Spiked CCA Concentration (ng/ml)") +
  ylab("G Score") +
  ggtitle("POC-CCA")

mean_trend <- G_score_long %>%
  group_by(V1) %>%
  summarise(mean_value = mean(Value + 1, na.rm = TRUE))
# Adding elements from the first plot
CCAplot <- base_plot +
  geom_line(data = G_score_long, aes(x = V1 * 0.03, y = Value + 1, group = name), linewidth = 1, color = "black", alpha = 0.1) +
  geom_line(data = mean_trend, aes(x = V1 * 0.03, y = mean_value), linewidth = 1.3, color = "grey", linetype = "dashed")
# Adding elements from the second plot with correct color and shape mappings
CCAplot <- CCAplot +
  annotation_logticks(sides = "b") +
  geom_point(data = Sim_long, aes(x = ngml * 0.03, y = G_score, color = Lab_Data, shape = Lab_Data), alpha = 0.5, size = 5) +
  scale_color_manual(values = c("Batch 1" = "#46c19a", "Batch 2" = "#7f62b8", "Batch 3" = "#b84c7d")) +
  scale_shape_manual(values = c("Batch 1" = 16, "Batch 2" = 17, "Batch 3" = 18))+
  labs(color = "", shape = "")

print(CCAplot)

library(cowplot)
library(gridExtra)

comb_plot <- plot_grid(
  CCAplot + theme(plot.title = element_text(size = rel(1.3)),     # Increase plot title size
                  axis.title.x = element_text(size = rel(1.2)),   # Increase x axis title size
                  axis.title.y = element_text(size = rel(1.2)),   # Increase y axis title size
                  axis.text.x = element_text(size = rel(1.2)),    # Increase x axis text size
                  axis.text.y = element_text(size = rel(1.2))),
  CCA3plot+ theme(plot.title = element_text(size = rel(1.3)),     # Increase plot title size
                  axis.title.x = element_text(size = rel(1.2)),   # Increase x axis title size
                  axis.title.y = element_text(size = rel(1.2)),   # Increase y axis title size
                  axis.text.x = element_text(size = rel(1.2)),    # Increase x axis text size
                  axis.text.y = element_text(size = rel(1.2)),
                  legend.title = element_text(size = 12),  # Adjust size for legend title
                  legend.text = element_text(size = 10)),
  ncol = 2, 
  nrow = 1,
  align = "h",
  rel_widths = c(1, 1)  
)
png("ng.png", 800, 350)
comb_plot
dev.off()


df$G_mean <- rowMeans(df[,2:4])
df$sd_mean <- rowMeans(df[,5:7])

plot3 <- ggplot(data = df, aes(x = ug_ml)) + 
  geom_line(aes(y = (G_mean), color = "Mean"), linewidth = 1) +
  geom_ribbon(aes(y = (G_mean), ymin = (G_mean - sd_mean), ymax = (G_mean + sd_mean), fill = "Mean", colour = "Mean"), alpha = .2) +
  scale_x_log10()+
  theme_bw() +
  xlab("ng/ml") + 
  ylab("G Score") + 
  ggtitle("recPOC-CCA")
