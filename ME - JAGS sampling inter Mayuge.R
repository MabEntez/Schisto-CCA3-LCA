rm(list = ls())
library(rstudioapi)
library(runjags)
library(ggplot2)
library(rjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

dt <- read.csv("CORNTD_recPOCCCA_2023_Mayuge.csv")

KK <- as.matrix(dt[,c(8,9,18,19,28,29)])

CCA1 <- as.matrix(dt[,c(38,40,42)]) - 1
CCA2 <- as.matrix(dt[,c(44)]) - 1
CCA3 <- as.matrix(dt[,c(46)]) - 1
CCA <- cbind(CCA1, CCA2, CCA3)

REC_CCA1 <- as.matrix(dt[,c(48,51,54)]) - 1
REC_CCA2 <- as.matrix(dt[,58]) - 1
REC_CCA3 <- as.matrix(dt[,62]) - 1
REC_CCA <- cbind(REC_CCA1, REC_CCA2, REC_CCA3)

CCA_batch <- as.matrix(dt[,c(39,41,43,45,47)])
REC_CCA_batch <- as.matrix(dt[,c(49,52,55,59,63)])

CCA_batch[CCA_batch == 220701075] <- 1
CCA_batch[CCA_batch == 220902098] <- 2
CCA_batch[CCA_batch == 221117133] <- 3

REC_CCA_cov <- as.matrix(read.csv("covariance_matrix.csv")[,2:3])
REC_CCA_m <- read.csv("mean_matrix.csv")[,2]

CCA_batch[is.na(CCA_batch)] <- 4
REC_CCA_batch[is.na(REC_CCA_batch)] <- 4

dtprior <- read.csv("kints.csv")

## Fit a multivariate normal distribution for the two parameters
library(MGMM)
fitParams <- FitGMM(data = as.matrix(dtprior[,c(4,5)]))

multim <- fitParams@Mean
multicov <- fitParams@Covariance

## Set seed ##
.RNG.seed <- function(chain)
  return( switch(chain, "1"= 1, "2"= 2) )
.RNG.name <- function(chain)
  return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )

## Inits ##
N = nrow(dt)
status = rep(1,N)
nCCA = ncol(CCA1)
nKK = ncol(KK)

m <- "model {
  # Prior random walk #
  
  for (n in 1:N){
  
    ## status
    status[n] ~ dbern(P)
    
    infect_intensity[n] ~ dgamma(sh,rt)
    
    ## KK component
     for (r in 1:nKK){
      KK[n,r] ~ dnegbin(k/(infect_intensity[n]*status[n] +k),k)
    }
    
    for(d in 1:3){
      CCAintensity[n,d,1] <- 0
      CCAintensity[n,d,2] <- (9 / (1 + exp(-multiParamCCA[1]*(infect_intensity[n]-multiParamCCA[2]))))
      
      RECCCAintensity[n,d,1] <- 0
      RECCCAintensity[n,d,2] <- (9 / (1 + exp(-multiParamRECCCA[1]*(infect_intensity[n]-multiParamRECCCA[2]))))
      
    X_CCA[n,d,1] ~ dnorm(CCAintensity[n,d,1], tauCCA_inter)
    X_CCA[n,d,2] ~ dnorm(CCAintensity[n,d,2], tauCCA_inter)
    X_RECCCA[n,d,1] ~ dnorm(RECCCAintensity[n,d,1], tauRECCCA_inter)
    X_RECCCA[n,d,2] ~ dnorm(RECCCAintensity[n,d,2], tauRECCCA_inter)
}


    ## CCA component
    for (r in 1:nCCA){
      CCA1[n,r] ~ dnorm(X_CCA[n,1,status[n]+1],tauCCA_batch[CCA_batch[n,r]]) 
    }
    CCA2[n,1] ~ dnorm(X_CCA[n,2,status[n]+1],tauCCA_batch[CCA_batch[n,4]])
    CCA3[n,1] ~ dnorm(X_CCA[n,3,status[n]+1],tauCCA_batch[CCA_batch[n,5]])
    
    ## RECCCA component
    for (r in 1:nCCA){
      REC_CCA1[n,r] ~ dnorm(X_RECCCA[n,1,status[n]+1],tauRECCCA_batch[REC_CCA_batch[n,r]]) 
    }
    REC_CCA2[n,1] ~ dnorm(X_RECCCA[n,2,status[n]+1],tauRECCCA_batch[REC_CCA_batch[n,4]])
    REC_CCA3[n,1] ~ dnorm(X_RECCCA[n,3,status[n]+1],tauRECCCA_batch[REC_CCA_batch[n,5]])

  }
  
  ## Priors ##
  P ~ dbeta(1, 1)
  k  ~ dgamma(0.001,0.001)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  
  tauCCA_inter ~ dgamma(0.001,0.001)
  tauRECCCA_inter ~ dgamma(0.001,0.001)
  tauCCA_batch[4] ~ dgamma(0.001,0.001) # DUMMY for NAs
  tauCCA_batch[1] ~ dgamma(0.001,0.001)
  tauCCA_batch[2] ~ dgamma(0.001,0.001)
  tauCCA_batch[3] ~ dgamma(0.001,0.001)
  tauRECCCA_batch[4] ~ dgamma(0.001,0.001) # DUMMY for NAs
  tauRECCCA_batch[1] ~ dgamma(0.001,0.001)
  tauRECCCA_batch[2] ~ dgamma(0.001,0.001)
  tauRECCCA_batch[3] ~ dgamma(0.001,0.001)
  
  multiParamCCA[1:2] ~ dmnorm.vcov(multim, multicov)
  multiParamRECCCA[1] ~ dnorm(multim[1],1/multicov[1,1])
  multiParamRECCCA[2] ~ dnorm(multim[2],1/multicov[2,2])
  
  #inits# .RNG.seed, .RNG.name, status
  #data# N, KK, CCA1, CCA2, CCA3, REC_CCA1, REC_CCA2, REC_CCA3, CCA_batch, REC_CCA_batch, nKK, nCCA, multim, multicov
  #monitor# P, sh, rt, tauCCA_inter, tauCCA_batch, tauRECCCA_inter, tauRECCCA_batch, k, multiParamCCA, multiParamRECCCA
}"


Results <- run.jags(m, burnin=1000, sample=5000, thin=1, n.chains=2, jags.refresh = 1, method = 'parallel',
                    plots = F, silent.jags = F)

plot(Results)

mcmc_samples <- as.mcmc(Results)

save.image("Results_inter_m.RData")

load("Results_inter_m.RData")

# Load necessary libraries
library(coda)
library(ggplot2)
library(reshape2)

# Extract relevant columns
m_CCAtau <- as.mcmc(Results)[, c(5, 6, 7)]
m_RECCCAtau <- as.mcmc(Results)[, c(10, 11, 12)]

# Convert data to a data frame for ggplot
m_CCAtau_df <- data.frame(m_CCAtau)
m_RECCCAtau_df <- data.frame(m_RECCCAtau)

write.csv(m_CCAtau_df, "m tau - CCA.csv", row.names = FALSE)
write.csv(m_RECCCAtau_df, "m tau - CCA3.csv", row.names = FALSE)

m_interCCAtau <- as.mcmc(Results)[, 4]
m_interRECCCAtau <- as.mcmc(Results)[, 9]

write.csv(m_interCCAtau, "m inter tau - CCA.csv", row.names = FALSE)
write.csv(m_interRECCCAtau, "m inter tau - CCA3.csv", row.names = FALSE)
