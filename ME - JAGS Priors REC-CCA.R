rm(list = ls())
library(rstudioapi)
library(dplyr)
library(tidyverse)
library(runjags)
library(ggplot2)
library(rjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

#### Jags Model ####

dt <- read.csv("Spiked REC-CCA.csv", header = FALSE)
colnames(dt) <- c("Dosage", "Batch_1", "Batch_2", "Batch_3")

REC_CCA1 <- dt[,2] - 1
REC_CCA2 <- dt[,3] - 1
REC_CCA3 <- dt[,4] - 1

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

#### Saving outputs ####

param_samples <- mcmc_samples[, c("multiParam[1]", "multiParam[2]")]
covariance_matrix <- cov(param_samples)
write.csv(covariance_matrix, "covariance_matrixCCA3.csv")
mean_matrix <- colMeans(mcmc_samples)[2:3]
write.csv(mean_matrix, "mean_matrixCCA3.csv")

#### tau plotting ####
tau_samples <- as.mcmc(results1, variable.names = "tau")[,1:3]


# Convert to a dataframe for ggplot
tau_df <- data.frame(tau_samples)

t1 <- ggplot() +
  geom_density(data=tau_df, aes(x= 1/ sqrt(tau1)), fill="#003f5c", alpha=0.5) +
  geom_density(data=tau_df, aes(x= 1/ sqrt(tau2)), fill="#bc5090", alpha=0.5) +
  geom_density(data=tau_df, aes(x= 1/ sqrt(tau3)), fill="#ffa600", alpha=0.5) +
  labs(x="Standard deviation", y="Density", title="POC-CCA3")+
  theme_bw()
t1

####Fitting gamma distribution####
library(fitdistrplus)
fit.gamma <- fitdist(tau_df[,2], distr = "gamma", method = "mle")
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

Sim <- read.csv("Spiked REC-CCA.csv", header = FALSE)
Sim <- Sim[,-5]
colnames(Sim) <- c("ngml", "Batch 1", "Batch 2", "Batch 3")
Sim_long <- pivot_longer(as.data.frame(Sim), cols = -"ngml", values_to = "G_score")
colnames(Sim_long) <- c("ngml", "Lab_Data", "G_score")

base_plot <- ggplot() +
  scale_x_log10() +
  theme_bw() +
  xlab("Spiked CCA Concentration (ng/ml)") +
  ylab("G Score") +
  ggtitle("POC-CCA3")

mean_trend <- G_score_long %>%
  group_by(V1) %>%
  summarise(mean_value = mean(Value + 1, na.rm = TRUE))

# Adding elements from the first plot
CCA3plot <- base_plot +
  geom_line(data = G_score_long, aes(x = V1 * 0.03, y = Value + 1, group = name), linewidth = 1, color = "black", alpha = 0.1) +
  geom_line(data = mean_trend, aes(x = V1 * 0.03, y = mean_value), linewidth = 1.3, color = "grey", linetype = "dashed")

# Adding elements from the second plot with correct color and shape mappings
CCA3plot <- CCA3plot +
  annotation_logticks(sides = "b") +
  geom_point(data = Sim_long, aes(x = ngml * 0.03, y = G_score, color = Lab_Data, shape = Lab_Data), alpha = 0.5, size = 5) +
  scale_color_manual(values = c("Batch 1" = "#b8533c", "Batch 2" = "#bc9b3c", "Batch 3" = "#6ca24d")) +
  scale_shape_manual(values = c("Batch 1" = 16, "Batch 2" = 17, "Batch 3" = 18))+
  labs(color = "", shape = "")

print(CCA3plot)
