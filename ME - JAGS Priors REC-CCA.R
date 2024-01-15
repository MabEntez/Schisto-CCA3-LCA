rm(list = ls())
library(rstudioapi)
library(dplyr)
library(tidyverse)
library(runjags)
library(ggplot2)
library(rjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

dt <- read.csv("Spiked REC-CCA.csv")
colnames(dt) <- c("Dosage", "Batch 1", "Batch 2", "Batch 3")

REC_CCA1 <- dt[,2]
REC_CCA2 <- dt[,3]
REC_CCA3 <- dt[,4]

Spiked <- as.vector(dt[,1])
Spiked <- c(Spiked, Spiked, Spiked)

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
  n = length(REC_CCA1),
  x1 = REC_CCA1,
  x2 = REC_CCA2,
  x3 = REC_CCA3,
  y = Spiked
)

method <- list(
  "burnin" = 3000,
  "sample" = 10000
)

# Run the model
results1 <- run.jags(model_str,
                    data = data_jags,
                    monitor = c("tau1", "tau2", "tau3", "multiParam[1]", "multiParam[2]"),
                    n.chains = 4,
                    burnin = method$burnin,
                    sample = method$sample
)

# Summary statistics
summary(results1)

# Plotting
plot(results1)

mcmc_samples <- as.mcmc(results1)

param_samples <- mcmc_samples[, c("multiParam[1]", "multiParam[2]")]
covariance_matrix <- cov(param_samples)
write.csv(covariance_matrix, "covariance_matrix.csv")
mean_matrix <- colMeans(mcmc_samples)[2:3]
write.csv(mean_matrix, "mean_matrix.csv")

tau_samples <- as.mcmc(results1, variable.names = "tau")[,1:3]


# Convert to a dataframe for ggplot
tau_df <- data.frame(tau_samples)

t1 <- ggplot() +
  geom_density(data=tau_df, aes(x=tau1), fill="#003f5c", alpha=0.5) +
  geom_density(data=tau_df, aes(x=tau2), fill="#bc5090", alpha=0.5) +
  geom_density(data=tau_df, aes(x=tau3), fill="#ffa600", alpha=0.5) +
  labs(x="Tau", y="Density", title="recCCA")+
  theme_bw()
