rm(list = ls())
library(rstudioapi)
library(dplyr)
library(tidyverse)
library(runjags)
library(ggplot2)
library(rjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

dt <- read.csv("Spiked CCA.csv")
dt <- dt[, 1:4]
colnames(dt) <- c("Dosage", "Batch 1", "Batch 2", "Batch 3")

CCA1 <- dt[, 2]
CCA2 <- dt[, 3]
CCA3 <- dt[, 4]

Spiked <- as.vector(dt[,1])

model_str = "model {
  for (i in 1:n){
    ## CCA component
    x1[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(y[i]-multiParam[2]))),tau1)
    x2[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(y[i]-multiParam[2]))),tau2)
    x3[i] ~ dnorm(9 / (1 + exp(-multiParam[1]*(y[i]-multiParam[2]))),tau3)
  }
  tau1 ~ dgamma(0.001,0.001)
  tau2 ~ dgamma(0.001,0.001)
  tau3 ~ dgamma(0.001,0.001)
  multiParam[1] ~ dnorm(0, 0.0001)  # non-informative prior for multiParam[1]
  multiParam[2] ~ dnorm(0, 0.0001)  # non-informative prior for multiParam[2]
}
"

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

tau_samples <- as.mcmc(results2, variable.names = "tau")[,1:3]

tau_df2 <- data.frame(tau_samples)

t2 <- ggplot() +
  geom_density(data=tau_df2, aes(x=tau1), fill="#003f5c", alpha=0.5) +
  geom_density(data=tau_df2, aes(x=tau2), fill="#bc5090", alpha=0.5) +
  geom_density(data=tau_df2, aes(x=tau3), fill="#ffa600", alpha=0.5) +
  labs(x="Tau", y="Density", title="CCA")+
  theme_bw()

library(cowplot)
library(gridExtra)

combined_plot <- plot_grid(
  t2+ theme(plot.title = element_text(size = rel(1.3)),
            axis.title.x = element_text(size = rel(1.2)),
            axis.title.y = element_text(size = rel(1.2)),
            axis.text.x = element_text(size = rel(1.3)),
            axis.text.y = element_text(size = rel(1.3)),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)) + 
    ylim(0, 2.75) +
    xlim(0, 6),
  t1+ theme(plot.title = element_text(size = rel(1.3)),
            axis.title.x = element_text(size = rel(1.2)),
            axis.title.y = element_text(size = rel(1.2)),
            axis.text.x = element_text(size = rel(1.3)),
            axis.text.y = element_text(size = rel(1.3)),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)) + 
    ylim(0, 2.75) +
    xlim(0, 6),
  ncol = 2, 
  nrow = 1,
  align = "h"
)
png("tau.png", 800, 340)
combined_plot
dev.off()
