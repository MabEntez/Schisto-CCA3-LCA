rm(list = ls())
library(rstudioapi)
library(runjags)
library(ggplot2)
library(rjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

dt <- read.csv("CORNTD_recPOCCCA_2023_Mayuge.csv")

KK <- as.matrix(dt[,c(8,9,18,19,28,29)])

CCA <- as.matrix(dt[,c(38,39,40,42,44)]) - 1

REC_CCA <- as.matrix(dt[,c(46,49,52,56,60)]) - 1

REC_CCA_cov <- as.matrix(read.csv("covariance_matrix.csv")[,2:3])
REC_CCA_m <- read.csv("mean_matrix.csv")[,2]

plot(rowMeans(REC_CCA, na.rm = TRUE), rowMeans(CCA, na.rm = TRUE)) 
plot(rowMeans(REC_CCA, na.rm = TRUE), rowMeans(KK, na.rm = TRUE)) 
plot(rowMeans(CCA, na.rm = TRUE), rowMeans(KK, na.rm = TRUE))
abline(a = 0, b = 1)

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
nCCA = ncol(CCA)
nKK = ncol(KK)
CCAintensity <- array(cbind(CCA,CCA),dim = c(N,5,2))
RECCCAintensity <- array(cbind(REC_CCA,REC_CCA),dim = c(N,5,2))

m <- "model {
  # Prior random walk #
  
  for (n in 1:N){
  
    ## status
    status[n] ~ dbern(P)
    
    infect_intensity[n,2] ~ dgamma(sh,rt)
    infect_intensity[n,1] <- 0
    
    ## KK component
     for (r in 1:nKK){
      KK[n,r] ~ dnegbin(k/(infect_intensity[n,status[n] + 1] +k),k)
    }
    
    for(d in 1:5){
      CCAintensity[n,d,1] <- 0
      CCAintensity[n,d,2] <- (9 / (1 + exp(-multiParamCCA[1]*(infect_intensity[n,status[n] + 1]-multiParamCCA[2]))))
      
      RECCCAintensity[n,d,1] <- 0
      RECCCAintensity[n,d,2] <- (9 / (1 + exp(-multiParamRECCCA[1]*(infect_intensity[n,status[n] + 1]-multiParamRECCCA[2]))))
      
      CCA[n,d] ~ dnorm(CCAintensity[n,d,status[n]+1],tauCCA_inter)T(0,9)
      REC_CCA[n,d] ~ dnorm(RECCCAintensity[n,d,status[n]+1],tauRECCCA_inter)T(0,9)
    }
  }
  
  ## Priors ##
  ssh ~ dgamma(0.001,0.001)
  rrt <- ssh/0.001
  mu ~ dunif(0,5)
  P ~ dbeta(1, 1)
  k  ~ dgamma(0.001,0.001)
  tauCCA ~ dgamma(0.001,0.001)
  tauRECCCA ~ dgamma(0.001,0.001)
  tauCCA_inter ~ dgamma(204.6194,956.8694)T(0,)
  tauCCA_intra ~ dnorm(1.24,10.70141)T(0,)
  tauRECCCA_inter ~ dgamma(204.6194,956.8694)T(0,) #dnorm(1.009404,1.128064)T(0,)
  tauRECCCA_intra ~ dnorm(1.55,6.786854)T(0,)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  multiParamCCA[1:2] ~ dmnorm.vcov(multim, multicov)
  multiParamRECCCA[1] ~ dnorm(multim[1],1/multicov[1,1])
  multiParamRECCCA[2] ~ dnorm(multim[2],1/multicov[2,2])
  
  #inits# .RNG.seed, .RNG.name, status 
  #, RECCCAintensity, CCAintensity
  #data# N, KK, CCA, REC_CCA, nKK, multim, multicov
  #monitor# P, sh, rt, tauCCA_inter, tauRECCCA_inter, k, multiParamCCA, multiParamRECCCA
}"

Results <- run.jags(m, burnin=1000, sample=5000, thin=1, n.chains=2, jags.refresh = 1, method = 'parallel',
                    plots = F, silent.jags = F)

plot(Results)

mcmc_samples <- as.mcmc(Results)

save.image("Results1.RData")

load("Results1.RData")

day1_data <- as.matrix(dt[, c(46, 49, 52)]) - 1
day2_data <- as.matrix(dt[, c(56)]) - 1
day3_data <- as.matrix(dt[, c(60)]) - 1

# Assuming that day1_data might contain multiple measurements per day,
# we'll take the mean for each subject for day 1
day1_mean <- rowMeans(day1_data, na.rm = T)

# Combine the day means for each subject into a new matrix
# where each row represents a subject and each column represents a day
all_days_data <- cbind(day1_mean, day2_data, day3_data)

# Calculate the standard deviation across the days for each subject
subject_sds <- apply(all_days_data, 1, sd)

# Now, if you want the mean of these standard deviations, you can calculate that as well
mean_subject_sd <- mean(subject_sds, na.rm = T)
sd_subject_sd <- sd(subject_sds, na.rm = T)
# Print the mean of the subject-wise standard deviations
print(mean_subject_sd)
print(sd_subject_sd)
# If you want to explore the individual standard deviations
print(subject_sds)

#### Sensitivity Analysis ####

## Draw samples from the posteriors
set.seed(1)
ndraws <- 1000
IDs <- sample.int(5000,ndraws)

## calculate vars for CCA
var_inter_RECCCA <- 1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs]
var_intra_RECCCA <- 1/c(Results$mcmc[[1]][,"tauRECCCA_intra"],Results$mcmc[[2]][,"tauRECCCA_intra"])[IDs]

mean(var_intra_RECCCA/var_inter_RECCCA, na.rm = TRUE) #0.26 times the variance of intra is higher than inter
sd(var_intra_RECCCA/var_inter_RECCCA, na.rm = TRUE) #0.022 sd in the mean

## calculate vars for CCA
var_inter_CCA <- 1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs]
var_intra_CCA <- 1/c(Results$mcmc[[1]][,"tauCCA_intra"],Results$mcmc[[2]][,"tauCCA_intra"])[IDs]

mean(var_intra_CCA/var_inter_CCA, na.rm = TRUE) #0.26 times the variance of intra is higher than inter
sd(var_intra_CCA/var_inter_CCA, na.rm = TRUE) #0.022 sd in the mean

## Calculate 95% quantiles
## KK * 24 to transform from eggs to Eggs per gram (EPGs)
quantile(var_inter_RECCCA,c(0.025,0.5,0.975), na.rm = TRUE)
quantile(var_intra_RECCCA,c(0.025,0.5,0.975), na.rm = TRUE)
quantile(var_inter_CCA,c(0.025,0.5,0.975), na.rm = TRUE)
quantile(var_intra_CCA,c(0.025,0.5,0.975), na.rm = TRUE)

####################################################################
## Do plot
library(scales)

pdf("SDKKandCCAMayuge.pdf",width=7*2,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfcol=c(1,2))

plot(NA,xlim = c(0,30),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)
polygon(density(var_intra_RECCCA)$x,density(var_intra_RECCCA)$y/max(density(var_intra_RECCCA)$y),
        col = alpha("grey",alpha=0.5), border = NA)
polygon(density(var_inter_RECCCA)$x,density(var_inter_RECCCA)$y/max(density(var_inter_RECCCA)$y),
        col = alpha("red",alpha=0.5), border = NA)

axis(1, at = seq(0,100,by=20), labels = seq(0,100,by=20))
axis(2)

mtext("Scaled Density",side=2,cex=1,line=1.2)
mtext("Variance of KK (epg)",side=1,cex=1,line=1.2)

legend("topright",c("Inter-day","Intra-day"),col=alpha(c("red","grey"),alpha=0.5),
       bty='n',cex=0.75,lty=1)

### G-Score
plot(NA,xlim = c(0,8),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)
polygon(density(var_intra_CCA)$x,density(var_intra_CCA)$y/max(density(var_intra_CCA)$y),
        col = alpha("grey",alpha=0.5), border = NA)
polygon(density(var_inter_CCA)$x,density(var_inter_CCA)$y/max(density(var_inter_CCA)$y),
        col = alpha("red",alpha=0.5), border = NA)

axis(1)
axis(2)

mtext("Scaled Density",side=2,cex=1,line=1.2)
mtext("Variance of POC-CCA",side=1,cex=1,line=1.2)

legend("topright",c("Inter-day","Intra-day"),col=alpha(c("red","grey"),alpha=0.5),
       bty='n',cex=0.75,lty=1)

dev.off()

#### Prevalence Figure ####

KKprev <- mean(sapply(1:nrow(KK1),function(x) 
{ifelse(sum(c(KK1[x,],
              KK2[x,],
              KK3[x,]))>0,1,0)}
)
, na.rm = T)

CCAprevH <- mean(sapply(1:nrow(KK1),function(x)
{ifelse(mean(c(CCA1[x,],
               CCA2[x,],
               CCA3[x,]), na.rm = T)>3,1,0)})
, na.rm = T)

CCAprevL <- mean(sapply(1:nrow(KK1),function(x)
{ifelse(mean(c(CCA1[x,],
               CCA2[x,],
               CCA3[x,]), na.rm = T)>2,1,0)})
, na.rm = T)

ModelP <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
ModelPrev <- quantile(ModelP,c(0.025,0.5,0.975))

pdf("PrevTororo.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim = c(0,2),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)

points(0,KKprev, pch=19,col="black")
points(1,CCAprevH, pch = 19, col = "red")
points(1,CCAprevL, pch = 19, col = "green")
points(2,ModelPrev[2], pch = 19)
arrows(2,ModelPrev[1],2,ModelPrev[3],angle=90,code = 3, length = 0.1)

axis(1, at= 0:2, labels = c("Kato-Katz", "POC-CCA", "Model"))
axis(2, at=seq(0,1,by=.2), labels = seq(0,1,by=.2)*100)

mtext("Prevalence (%)",side=2,cex=1,line=1.2)
#mtext("Variance of G-score",side=1,cex=1,line=1.2)

dev.off()