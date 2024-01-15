###### Simulations to evaluate the error #####
## By M Entezami + JM Prada
rm(list = ls())

## Load library and set workspace
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

## Load data
load("Results2.RData")
Results2 <- Results
load("Results1.RData")
source("SimulationErrorFunctions_ME.R")

## Set draws
set.seed(100)
ndraws <- 1000
IDs <- sample.int(10000,ndraws)

#### Mayuge
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]

## Simulate prevalence estimates for Mayuge
PrevEstimatesMayuge <- sapply(1:length(IDs),function(x){modelData(1000,x)})

## Calculate quantiles for squared error
sderH <- apply((P-PrevEstimatesMayuge)^2,1,quantile,c(.025,.5,.975))

#### Tororo
interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauCCA_inter"],Results2$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]

## Simulate prevalence estimates for Tororo
PrevEstimatesTororo <- sapply(1:length(IDs),function(x){modelData(1000,x)})

## Calculate squared error
sderL <- apply((P-PrevEstimatesTororo)^2,1,quantile,c(.025,.5,.975))


######################################################
### Plot
pdf("ErrorPlot.pdf",width=7*2,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfcol=c(1,2))

plot(NA,xlim = c(1,12),ylim = c(0,0.1),axes = F, xlab = "", ylab = "")
points(sderH[2,], pch=19)
arrows(x0=1:12,y0=sderH[1,],y1=sderH[3,],angle = 90, code = 3, length = .1)
abline(v=3.5, lty = "dashed")
abline(v=6.5, lty = "dashed")
abline(v=9.5, lty = "dashed")

axis(1, at = 1:12, labels = rep(1:3,4))
axis(2)

mtext("Squared error",side=2,cex=1,line=1.2)
mtext("Number of days of sampling",side=1,cex=1,line=1.2)

text(2,0.095,labels = "Threshold = G2")
text(5,0.095,labels = "Threshold = G2.5")
text(8,0.095,labels = "Threshold = G3")
text(11,0.095,labels = "Threshold = G4")

#### Tororo
plot(NA,xlim = c(1,12),ylim = c(0,0.1),axes = F, xlab = "", ylab = "")
points(sderL[2,], pch=19)
arrows(x0=1:12,y0=sderL[1,],y1=sderL[3,],angle = 90, code = 3, length = .1)
abline(v=3.5, lty = "dashed")
abline(v=6.5, lty = "dashed")
abline(v=9.5, lty = "dashed")

axis(1, at = 1:12, labels = rep(1:3,4))
axis(2)

mtext("Squared error",side=2,cex=1,line=1.2)
mtext("Number of days of sampling",side=1,cex=1,line=1.2)

text(2,0.095,labels = "Threshold = G2")
text(5,0.095,labels = "Threshold = G2.5")
text(8,0.095,labels = "Threshold = G3")
text(11,0.095,labels = "Threshold = G4")

dev.off()


#############################################
#### Mayuge Sens/Spec
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]
multim1_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[1]"],Results$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[2]"],Results$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[1]"],Results$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[2]"],Results$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

## Simulate Sensitivity and Specificity for Mayuge 
SensSpecEstimatesMayugeCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
SensSpecEstimatesMayugeREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(1000,x)})
## Calculate quantiles
QsenspcMayugeCCA <- apply(SensSpecEstimatesMayugeCCA,1,quantile,c(.025,.5,.975))
QsenspcMayugeREC <- apply(SensSpecEstimatesMayugeREC,1,quantile,c(.025,.5,.975))

#### Tororo
interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]
multim1_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[1]"],Results2$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[2]"],Results2$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[1]"],Results2$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[2]"],Results2$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

## Simulate Sensitivity and Specificity for Tororo
SensSpecEstimatesTororoCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
SensSpecEstimatesTororoREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(1000,x)})
## Calculate quantiles
QsenspcTororoCCA <- apply(SensSpecEstimatesTororoCCA,1,quantile,c(.025,.5,.975))
QsenspcTororoREC <- apply(SensSpecEstimatesTororoREC,1,quantile,c(.025,.5,.975))
### Plot ROC curve ###

png("Figure 2 - ROCcurve.png",width=800,height=500)
par(font=2, cex.axis=1, lwd=2, mar=c(4, 4, 2, 1), mgp=c(3, 0.4, 0)) # Adjust margins
par(mfcol=c(2,2))

## Mayuge CCA
plot(NA,xlim = c(0,.6),ylim = c(.4,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcMayugeCCA[2,c(2,4,6,8)],QsenspcMayugeCCA[2,c(1,3,5,7)], pch=19)
lines(1-QsenspcMayugeCCA[2,c(2,4,6,8)],QsenspcMayugeCCA[2,c(1,3,5,7)])
arrows(x0=1-QsenspcMayugeCCA[2,c(2,4,6,8)],y0=QsenspcMayugeCCA[1,c(1,3,5,7)],
       y1=QsenspcMayugeCCA[3,c(1,3,5,7)], angle=90, code = 3,length = .05)
arrows(x0=1-QsenspcMayugeCCA[1,c(2,4,6,8)],y0=QsenspcMayugeCCA[2,c(1,3,5,7)],
       x1=1-QsenspcMayugeCCA[3,c(2,4,6,8)], angle=90, code = 3,length = .05)

## 2 Days
points(1-QsenspcMayugeCCA[2,c(10,12,14,16)],QsenspcMayugeCCA[2,c(9,11,13,15)], pch=19, col="red")
lines(1-QsenspcMayugeCCA[2,c(10,12,14,16)],QsenspcMayugeCCA[2,c(9,11,13,15)], col="red")
arrows(x0=1-QsenspcMayugeCCA[2,c(10,12,14,16)],y0=QsenspcMayugeCCA[1,c(9,11,13,15)],
       y1=QsenspcMayugeCCA[3,c(9,11,13,15)], angle=90, code = 3,length = .05, col="red")
arrows(x0=1-QsenspcMayugeCCA[1,c(10,12,14,16)],y0=QsenspcMayugeCCA[2,c(9,11,13,15)],
       x1=1-QsenspcMayugeCCA[3,c(10,12,14,16)], angle=90, code = 3,length = .05, col="red")

## 3 Days
points(1-QsenspcMayugeCCA[2,c(18,20,22,24)],QsenspcMayugeCCA[2,c(17,19,21,23)], pch=19, col="darkgreen")
lines(1-QsenspcMayugeCCA[2,c(18,20,22,24)],QsenspcMayugeCCA[2,c(17,19,21,23)], col="darkgreen")
arrows(x0=1-QsenspcMayugeCCA[2,c(18,20,22,24)],y0=QsenspcMayugeCCA[1,c(17,19,21,23)],
       y1=QsenspcMayugeCCA[3,c(17,19,21,23)], angle=90, code = 3,length = .05, col="darkgreen")
arrows(x0=1-QsenspcMayugeCCA[1,c(18,20,22,24)],y0=QsenspcMayugeCCA[2,c(17,19,21,23)],
       x1=1-QsenspcMayugeCCA[3,c(18,20,22,24)], angle=90, code = 3,length = .05, col="darkgreen")

axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)

mtext("Sensitivity",side=2,cex=1.25,line=2)
mtext("1 - Specificity",side=1,cex=1.25,line=2)
title("A", cex.main=1.25)

## Mayuge REC
plot(NA,xlim = c(0,.6),ylim = c(.4,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcMayugeREC[2,c(2,4,6,8)],QsenspcMayugeREC[2,c(1,3,5,7)], pch=19)
lines(1-QsenspcMayugeREC[2,c(2,4,6,8)],QsenspcMayugeREC[2,c(1,3,5,7)])
arrows(x0=1-QsenspcMayugeREC[2,c(2,4,6,8)],y0=QsenspcMayugeREC[1,c(1,3,5,7)],
       y1=QsenspcMayugeREC[3,c(1,3,5,7)], angle=90, code = 3,length = .05)
arrows(x0=1-QsenspcMayugeREC[1,c(2,4,6,8)],y0=QsenspcMayugeREC[2,c(1,3,5,7)],
       x1=1-QsenspcMayugeREC[3,c(2,4,6,8)], angle=90, code = 3,length = .05)

## 2 Days
points(1-QsenspcMayugeREC[2,c(10,12,14,16)],QsenspcMayugeREC[2,c(9,11,13,15)], pch=19, col="red")
lines(1-QsenspcMayugeREC[2,c(10,12,14,16)],QsenspcMayugeREC[2,c(9,11,13,15)], col="red")
arrows(x0=1-QsenspcMayugeREC[2,c(10,12,14,16)],y0=QsenspcMayugeREC[1,c(9,11,13,15)],
       y1=QsenspcMayugeREC[3,c(9,11,13,15)], angle=90, code = 3,length = .05, col="red")
arrows(x0=1-QsenspcMayugeREC[1,c(10,12,14,16)],y0=QsenspcMayugeREC[2,c(9,11,13,15)],
       x1=1-QsenspcMayugeREC[3,c(10,12,14,16)], angle=90, code = 3,length = .05, col="red")

## 3 Days
points(1-QsenspcMayugeREC[2,c(18,20,22,24)],QsenspcMayugeREC[2,c(17,19,21,23)], pch=19, col="darkgreen")
lines(1-QsenspcMayugeREC[2,c(18,20,22,24)],QsenspcMayugeREC[2,c(17,19,21,23)], col="darkgreen")
arrows(x0=1-QsenspcMayugeREC[2,c(18,20,22,24)],y0=QsenspcMayugeREC[1,c(17,19,21,23)],
       y1=QsenspcMayugeREC[3,c(17,19,21,23)], angle=90, code = 3,length = .05, col="darkgreen")
arrows(x0=1-QsenspcMayugeREC[1,c(18,20,22,24)],y0=QsenspcMayugeREC[2,c(17,19,21,23)],
       x1=1-QsenspcMayugeREC[3,c(18,20,22,24)], angle=90, code = 3,length = .05, col="darkgreen")

axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)

mtext("Sensitivity",side=2,cex=1.25,line=2)
mtext("1 - Specificity",side=1,cex=1.25,line=2)
title("C", cex.main=1.25)

## Tororo
plot(NA,xlim = c(0,.6),ylim = c(.4,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcTororoCCA[2,c(2,4,6,8)],QsenspcTororoCCA[2,c(1,3,5,7)], pch=19)
lines(1-QsenspcTororoCCA[2,c(2,4,6,8)],QsenspcTororoCCA[2,c(1,3,5,7)])
arrows(x0=1-QsenspcTororoCCA[2,c(2,4,6,8)],y0=QsenspcTororoCCA[1,c(1,3,5,7)],
       y1=QsenspcTororoCCA[3,c(1,3,5,7)], angle=90, code = 3,length = .05)
arrows(x0=1-QsenspcTororoCCA[1,c(2,4,6,8)],y0=QsenspcTororoCCA[2,c(1,3,5,7)],
       x1=1-QsenspcTororoCCA[3,c(2,4,6,8)], angle=90, code = 3,length = .05)

## 2 Days
points(1-QsenspcTororoCCA[2,c(10,12,14,16)],QsenspcTororoCCA[2,c(9,11,13,15)], pch=19, col="red")
lines(1-QsenspcTororoCCA[2,c(10,12,14,16)],QsenspcTororoCCA[2,c(9,11,13,15)], col="red")
arrows(x0=1-QsenspcTororoCCA[2,c(10,12,14,16)],y0=QsenspcTororoCCA[1,c(9,11,13,15)],
       y1=QsenspcTororoCCA[3,c(9,11,13,15)], angle=90, code = 3,length = .05, col="red")
arrows(x0=1-QsenspcTororoCCA[1,c(10,12,14,16)],y0=QsenspcTororoCCA[2,c(9,11,13,15)],
       x1=1-QsenspcTororoCCA[3,c(10,12,14,16)], angle=90, code = 3,length = .05, col="red")

## 3 Days
points(1-QsenspcTororoCCA[2,c(18,20,22,24)],QsenspcTororoCCA[2,c(17,19,21,23)], pch=19, col="darkgreen")
lines(1-QsenspcTororoCCA[2,c(18,20,22,24)],QsenspcTororoCCA[2,c(17,19,21,23)], col="darkgreen")
arrows(x0=1-QsenspcTororoCCA[2,c(18,20,22,24)],y0=QsenspcTororoCCA[1,c(17,19,21,23)],
       y1=QsenspcTororoCCA[3,c(17,19,21,23)], angle=90, code = 3,length = .05, col="darkgreen")
arrows(x0=1-QsenspcTororoCCA[1,c(18,20,22,24)],y0=QsenspcTororoCCA[2,c(17,19,21,23)],
       x1=1-QsenspcTororoCCA[3,c(18,20,22,24)], angle=90, code = 3,length = .05, col="darkgreen")

axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)

mtext("Sensitivity",side=2,cex=1.25,line=2)
mtext("1 - Specificity",side=1,cex=1.25,line=2)
title("B", cex.main=1.25)

## Tororo
plot(NA,xlim = c(0,.6),ylim = c(.4,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcTororoREC[2,c(2,4,6,8)],QsenspcTororoREC[2,c(1,3,5,7)], pch=19)
lines(1-QsenspcTororoREC[2,c(2,4,6,8)],QsenspcTororoREC[2,c(1,3,5,7)])
arrows(x0=1-QsenspcTororoREC[2,c(2,4,6,8)],y0=QsenspcTororoREC[1,c(1,3,5,7)],
       y1=QsenspcTororoREC[3,c(1,3,5,7)], angle=90, code = 3,length = .05)
arrows(x0=1-QsenspcTororoREC[1,c(2,4,6,8)],y0=QsenspcTororoREC[2,c(1,3,5,7)],
       x1=1-QsenspcTororoREC[3,c(2,4,6,8)], angle=90, code = 3,length = .05)

## 2 Days
points(1-QsenspcTororoREC[2,c(10,12,14,16)],QsenspcTororoREC[2,c(9,11,13,15)], pch=19, col="red")
lines(1-QsenspcTororoREC[2,c(10,12,14,16)],QsenspcTororoREC[2,c(9,11,13,15)], col="red")
arrows(x0=1-QsenspcTororoREC[2,c(10,12,14,16)],y0=QsenspcTororoREC[1,c(9,11,13,15)],
       y1=QsenspcTororoREC[3,c(9,11,13,15)], angle=90, code = 3,length = .05, col="red")
arrows(x0=1-QsenspcTororoREC[1,c(10,12,14,16)],y0=QsenspcTororoREC[2,c(9,11,13,15)],
       x1=1-QsenspcTororoREC[3,c(10,12,14,16)], angle=90, code = 3,length = .05, col="red")

## 3 Days
points(1-QsenspcTororoREC[2,c(18,20,22,24)],QsenspcTororoREC[2,c(17,19,21,23)], pch=19, col="darkgreen")
lines(1-QsenspcTororoREC[2,c(18,20,22,24)],QsenspcTororoREC[2,c(17,19,21,23)], col="darkgreen")
arrows(x0=1-QsenspcTororoREC[2,c(18,20,22,24)],y0=QsenspcTororoREC[1,c(17,19,21,23)],
       y1=QsenspcTororoREC[3,c(17,19,21,23)], angle=90, code = 3,length = .05, col="darkgreen")
arrows(x0=1-QsenspcTororoREC[1,c(18,20,22,24)],y0=QsenspcTororoREC[2,c(17,19,21,23)],
       x1=1-QsenspcTororoREC[3,c(18,20,22,24)], angle=90, code = 3,length = .05, col="darkgreen")

axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)

mtext("Sensitivity",side=2,cex=1.25,line=2)
mtext("1 - Specificity",side=1,cex=1.25,line=2)
title("D", cex.main=1.25)
legend("bottomright",c("1 Day","2 Days","3 Days"),
       col=c("black","red","darkgreen"),
       bty='n',cex=1.25,lty=1)

dev.off()

#### Mayuge Sens/Spec
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]
multim1_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[1]"],Results$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[2]"],Results$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[1]"],Results$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[2]"],Results$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

AUCMayugeCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
AUCMayugeREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(1000,x)})

interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]
multim1_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[1]"],Results2$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[2]"],Results2$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[1]"],Results2$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[2]"],Results2$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

AUCTororoCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
AUCTororoREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(1000,x)})


#### AUC calculation ####

#### Mayuge AUC ####
### CCA 
## 1 Day 
# Initialize an empty vector to store AUC values
AUCMayugeCCA <- t(AUCMayugeCCA)
auc_values <- numeric(nrow(AUCMayugeCCA))
M_CCA_1 <- c()
# Loop through each repeat (row)
for (j in 1:nrow(AUCMayugeCCA)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCMayugeCCA[j, c(1, 3, 5, 7)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeCCA[j, c(2, 4, 6, 8)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  M_CCA_1[j] <- auc
}


## 2 Day 
# Loop through each repeat (row)
M_CCA_2 <- c()
for (j in 1:nrow(AUCMayugeCCA)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCMayugeCCA[j, c(9,11,13,15)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeCCA[j, c(10,12,14,16)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  M_CCA_2[j] <- auc
}

## 3 Day 
# Loop through each repeat (row)
M_CCA_3 <- c()
for (j in 1:nrow(AUCMayugeCCA)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCMayugeCCA[j, c(17,19,21,23)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeCCA[j, c(18,20,22,24)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  M_CCA_3[j] <- auc
}

### REC 
## 1 Day 
# Initialize an empty vector to store AUC values
M_REC_1 <- c()
AUCMayugeREC <- t(AUCMayugeREC)
auc_values <- numeric(nrow(AUCMayugeREC))
# Loop through each repeat (row)
for (j in 1:nrow(AUCMayugeREC)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCMayugeREC[j, c(1, 3, 5, 7)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeREC[j, c(2, 4, 6, 8)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  M_REC_1[j] <- auc
}


## 2 Day 
# Loop through each repeat (row)
M_REC_2 <- c()
for (j in 1:nrow(AUCMayugeREC)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCMayugeREC[j, c(9,11,13,15)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeREC[j, c(10,12,14,16)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  M_REC_2[j] <- auc
}

## 3 Day 
# Loop through each repeat (row)
M_REC_3 <- c()
for (j in 1:nrow(AUCMayugeREC)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCMayugeREC[j, c(17,19,21,23)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeREC[j, c(18,20,22,24)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  M_REC_3[j] <- auc
}

#### Tororo AUC ####
### CCA 
## 1 Day 
# Initialize an empty vector to store AUC values
AUCTororoCCA <- t(AUCTororoCCA)
auc_values <- numeric(nrow(AUCTororoCCA))
T_CCA_1 <- c()
# Loop through each repeat (row)
for (j in 1:nrow(AUCTororoCCA)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCTororoCCA[j, c(1, 3, 5, 7)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoCCA[j, c(2, 4, 6, 8)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  T_CCA_1[j] <- auc
}


## 2 Day 
# Loop through each repeat (row)
T_CCA_2 <- c()
for (j in 1:nrow(AUCTororoCCA)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCTororoCCA[j, c(9,11,13,15)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoCCA[j, c(10,12,14,16)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  T_CCA_2[j] <- auc
}

## 3 Day 
# Loop through each repeat (row)
T_CCA_3 <- c()
for (j in 1:nrow(AUCTororoCCA)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCTororoCCA[j, c(17,19,21,23)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoCCA[j, c(18,20,22,24)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  T_CCA_3[j] <- auc
}

### REC 
## 1 Day 
# Initialize an empty vector to store AUC values
AUCTororoREC <- t(AUCTororoREC)
auc_values <- numeric(nrow(AUCTororoREC))
T_REC_1 <- c()
# Loop through each repeat (row)
for (j in 1:nrow(AUCTororoREC)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCTororoREC[j, c(1, 3, 5, 7)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoREC[j, c(2, 4, 6, 8)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  T_REC_1[j] <- auc
}


## 2 Day 
# Loop through each repeat (row)
T_REC_2 <- c()
for (j in 1:nrow(AUCTororoREC)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCTororoREC[j, c(9,11,13,15)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoREC[j, c(10,12,14,16)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  T_REC_2[j] <- auc
}

## 3 Day 
# Loop through each repeat (row)
T_REC_3 <- c()
for (j in 1:nrow(AUCTororoREC)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCTororoREC[j, c(17,19,21,23)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoREC[j, c(18,20,22,24)], 1)
  
  # Sort the data by 1-specificity if not already sorted
  sorted_indices = order(one_minus_specificity)
  sensitivity = sensitivity[sorted_indices]
  one_minus_specificity = one_minus_specificity[sorted_indices]
  
  # Calculate AUC using the trapezoidal rule for this repeat
  auc = 0
  for (i in 2:length(sensitivity)) {
    height = one_minus_specificity[i] - one_minus_specificity[i - 1]
    avg_base = (sensitivity[i] + sensitivity[i - 1]) / 2
    auc = auc + (height * avg_base)
  }
  
  # Store the calculated AUC
  T_REC_3[j] <- auc
}

#### Format AUC ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(scales)

QAUCMayugeCCA <- apply(t(data.frame(M_CCA_1, M_CCA_2, M_CCA_3)),1,quantile,c(.025,.5,.975))
QAUCMayugeREC <- apply(t(data.frame(M_REC_1, M_REC_2, M_REC_3)),1,quantile,c(.025,.5,.975))
QAUCTororoCCA <- apply(t(data.frame(T_CCA_1, T_CCA_2, T_CCA_3)),1,quantile,c(.025,.5,.975))
QAUCTororoREC <- apply(t(data.frame(T_REC_1, T_REC_2, T_REC_3)),1,quantile,c(.025,.5,.975))

### Target profile ####
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]
multim1_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[1]"],Results$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[2]"],Results$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[1]"],Results$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[2]"],Results$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

TPMayugeCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(100,x)})
TPMayugeREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(100,x)})
TPMayugeCCA <- t(TPMayugeCCA)
TPMayugeREC <- t(TPMayugeREC)

interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]
multim1_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[1]"],Results2$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[2]"],Results2$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[1]"],Results2$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[2]"],Results2$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

TPTororoCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(100,x)})
TPTororoREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(100,x)})
TPTororoCCA <- t(TPTororoCCA)
TPTororoREC <- t(TPTororoREC)


# Condition values
x <- 0.6 # Your threshold for column1
y <- 0.95 # Your threshold for column2

# Calculate the percentage
MC_p_2 <- sum(TPMayugeCCA[,1] >= x & TPMayugeCCA[,2] >= y) / nrow(TPMayugeCCA) * 100
MC_p_25 <- sum(TPMayugeCCA[,3] >= x & TPMayugeCCA[,4] >= y) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,5] >= x & TPMayugeCCA[,6] >= y) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,7] >= x & TPMayugeCCA[,8] >= y) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,1] >= x & TPMayugeREC[,2] >= y) / nrow(TPMayugeREC) * 100
MR_p_25 <- sum(TPMayugeREC[,3] >= x & TPMayugeREC[,4] >= y) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,5] >= x & TPMayugeREC[,6] >= y) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,7] >= x & TPMayugeREC[,8] >= y) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,1] >= x & TPTororoCCA[,2] >= y) / nrow(TPTororoCCA) * 100
TC_p_25 <- sum(TPTororoCCA[,3] >= x & TPTororoCCA[,4] >= y) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,5] >= x & TPTororoCCA[,6] >= y) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,7] >= x & TPTororoCCA[,8] >= y) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,1] >= x & TPTororoREC[,2] >= y) / nrow(TPTororoREC) * 100
TR_p_25 <- sum(TPTororoREC[,3] >= x & TPTororoREC[,4] >= y) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,5] >= x & TPTororoREC[,6] >= y) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,7] >= x & TPTororoREC[,8] >= y) / nrow(TPTororoREC) * 100

d1_TPP <- data.frame(c(MC_p_2, MC_p_25, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_25, MR_p_3, MR_p_4),
                     c(TC_p_2, TC_p_25, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_25, TR_p_3, TR_p_4),
                     c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d1_TPP) <- c("A", "B", "C", "D", "Treshold")


d1_TPP <- melt(d1_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d1_m <- ggplot(d1_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "2.5 G score threshold" = "grey33",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "1 days sampling", 
       y = "Percentage achieving TPP", 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

MC_p_2 <- sum(TPMayugeCCA[,9] >= x & TPMayugeCCA[,10] >= y) / nrow(TPMayugeCCA) * 100
MC_p_25 <- sum(TPMayugeCCA[,11] >= x & TPMayugeCCA[,12] >= y) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,13] >= x & TPMayugeCCA[,14] >= y) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,15] >= x & TPMayugeCCA[,16] >= y) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,9] >= x & TPMayugeREC[,10] >= y) / nrow(TPMayugeREC) * 100
MR_p_25 <- sum(TPMayugeREC[,11] >= x & TPMayugeREC[,12] >= y) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,13] >= x & TPMayugeREC[,14] >= y) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,15] >= x & TPMayugeREC[,16] >= y) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,9] >= x & TPTororoCCA[,10] >= y) / nrow(TPTororoCCA) * 100
TC_p_25 <- sum(TPTororoCCA[,11] >= x & TPTororoCCA[,12] >= y) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,13] >= x & TPTororoCCA[,14] >= y) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,15] >= x & TPTororoCCA[,16] >= y) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,9] >= x & TPTororoREC[,10] >= y) / nrow(TPTororoREC) * 100
TR_p_25 <- sum(TPTororoREC[,11] >= x & TPTororoREC[,12] >= y) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,13] >= x & TPTororoREC[,14] >= y) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,15] >= x & TPTororoREC[,16] >= y) / nrow(TPTororoREC) * 100

d2_TPP <- data.frame(c(MC_p_2, MC_p_25, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_25, MR_p_3, MR_p_4),
                     c(TC_p_2, TC_p_25, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_25, TR_p_3, TR_p_4),
                     c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d2_TPP) <- c("A", "B", "C", "D", "Treshold")

d2_TPP <- melt(d2_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d2_m <- ggplot(d2_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "2.5 G score threshold" = "grey33",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "2 days sampling", 
       y = "Percentage achieving TPP", 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

MC_p_2 <- sum(TPMayugeCCA[,17] >= x & TPMayugeCCA[,18] >= y) / nrow(TPMayugeCCA) * 100
MC_p_25 <- sum(TPMayugeCCA[,19] >= x & TPMayugeCCA[,20] >= y) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,21] >= x & TPMayugeCCA[,22] >= y) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,23] >= x & TPMayugeCCA[,24] >= y) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,17] >= x & TPMayugeREC[,18] >= y) / nrow(TPMayugeREC) * 100
MR_p_25 <- sum(TPMayugeREC[,19] >= x & TPMayugeREC[,20] >= y) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,21] >= x & TPMayugeREC[,22] >= y) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,23] >= x & TPMayugeREC[,24] >= y) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,17] >= x & TPTororoCCA[,18] >= y) / nrow(TPTororoCCA) * 100
TC_p_25 <- sum(TPTororoCCA[,19] >= x & TPTororoCCA[,20] >= y) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,21] >= x & TPTororoCCA[,22] >= y) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,23] >= x & TPTororoCCA[,24] >= y) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,17] >= x & TPTororoREC[,18] >= y) / nrow(TPTororoREC) * 100
TR_p_25 <- sum(TPTororoREC[,19] >= x & TPTororoREC[,20] >= y) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,21] >= x & TPTororoREC[,22] >= y) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,23] >= x & TPTororoREC[,24] >= y) / nrow(TPTororoREC) * 100

d3_TPP <- data.frame(c(MC_p_2, MC_p_25, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_25, MR_p_3, MR_p_4),
                     c(TC_p_2, TC_p_25, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_25, TR_p_3, TR_p_4),
                     c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d3_TPP) <- c("A", "B", "C", "D", "Treshold")

d3_TPP <- melt(d3_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d3_m <- ggplot(d3_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "2.5 G score threshold" = "grey33",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "3 days sampling", 
       y = "Percentage achieving TPP", 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

# Condition values
x <- 0.75 # Your threshold for column1
y <- 0.965 # Your threshold for column2

# Calculate the percentage
MC_p_2 <- sum(TPMayugeCCA[,1] >= x & TPMayugeCCA[,2] >= y) / nrow(TPMayugeCCA) * 100
MC_p_25 <- sum(TPMayugeCCA[,3] >= x & TPMayugeCCA[,4] >= y) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,5] >= x & TPMayugeCCA[,6] >= y) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,7] >= x & TPMayugeCCA[,8] >= y) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,1] >= x & TPMayugeREC[,2] >= y) / nrow(TPMayugeREC) * 100
MR_p_25 <- sum(TPMayugeREC[,3] >= x & TPMayugeREC[,4] >= y) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,5] >= x & TPMayugeREC[,6] >= y) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,7] >= x & TPMayugeREC[,8] >= y) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,1] >= x & TPTororoCCA[,2] >= y) / nrow(TPTororoCCA) * 100
TC_p_25 <- sum(TPTororoCCA[,3] >= x & TPTororoCCA[,4] >= y) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,5] >= x & TPTororoCCA[,6] >= y) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,7] >= x & TPTororoCCA[,8] >= y) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,1] >= x & TPTororoREC[,2] >= y) / nrow(TPTororoREC) * 100
TR_p_25 <- sum(TPTororoREC[,3] >= x & TPTororoREC[,4] >= y) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,5] >= x & TPTororoREC[,6] >= y) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,7] >= x & TPTororoREC[,8] >= y) / nrow(TPTororoREC) * 100

d1_TPP <- data.frame(c(MC_p_2, MC_p_25, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_25, MR_p_3, MR_p_4),
                     c(TC_p_2, TC_p_25, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_25, TR_p_3, TR_p_4),
                     c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d1_TPP) <- c("A", "B", "C", "D", "Treshold")


d1_TPP <- melt(d1_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d1_o <- ggplot(d1_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "2.5 G score threshold" = "grey33",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "1 days sampling", 
       y = "Percentage achieving TPP", 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

MC_p_2 <- sum(TPMayugeCCA[,9] >= x & TPMayugeCCA[,10] >= y) / nrow(TPMayugeCCA) * 100
MC_p_25 <- sum(TPMayugeCCA[,11] >= x & TPMayugeCCA[,12] >= y) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,13] >= x & TPMayugeCCA[,14] >= y) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,15] >= x & TPMayugeCCA[,16] >= y) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,9] >= x & TPMayugeREC[,10] >= y) / nrow(TPMayugeREC) * 100
MR_p_25 <- sum(TPMayugeREC[,11] >= x & TPMayugeREC[,12] >= y) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,13] >= x & TPMayugeREC[,14] >= y) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,15] >= x & TPMayugeREC[,16] >= y) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,9] >= x & TPTororoCCA[,10] >= y) / nrow(TPTororoCCA) * 100
TC_p_25 <- sum(TPTororoCCA[,11] >= x & TPTororoCCA[,12] >= y) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,13] >= x & TPTororoCCA[,14] >= y) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,15] >= x & TPTororoCCA[,16] >= y) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,9] >= x & TPTororoREC[,10] >= y) / nrow(TPTororoREC) * 100
TR_p_25 <- sum(TPTororoREC[,11] >= x & TPTororoREC[,12] >= y) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,13] >= x & TPTororoREC[,14] >= y) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,15] >= x & TPTororoREC[,16] >= y) / nrow(TPTororoREC) * 100

d2_TPP <- data.frame(c(MC_p_2, MC_p_25, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_25, MR_p_3, MR_p_4),
                     c(TC_p_2, TC_p_25, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_25, TR_p_3, TR_p_4),
                     c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d2_TPP) <- c("A", "B", "C", "D", "Treshold")

d2_TPP <- melt(d2_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d2_o <- ggplot(d2_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "2.5 G score threshold" = "grey33",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "2 days sampling", 
       y = "Percentage achieving TPP", 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

MC_p_2 <- sum(TPMayugeCCA[,17] >= x & TPMayugeCCA[,18] >= y) / nrow(TPMayugeCCA) * 100
MC_p_25 <- sum(TPMayugeCCA[,19] >= x & TPMayugeCCA[,20] >= y) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,21] >= x & TPMayugeCCA[,22] >= y) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,23] >= x & TPMayugeCCA[,24] >= y) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,17] >= x & TPMayugeREC[,18] >= y) / nrow(TPMayugeREC) * 100
MR_p_25 <- sum(TPMayugeREC[,19] >= x & TPMayugeREC[,20] >= y) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,21] >= x & TPMayugeREC[,22] >= y) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,23] >= x & TPMayugeREC[,24] >= y) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,17] >= x & TPTororoCCA[,18] >= y) / nrow(TPTororoCCA) * 100
TC_p_25 <- sum(TPTororoCCA[,19] >= x & TPTororoCCA[,20] >= y) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,21] >= x & TPTororoCCA[,22] >= y) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,23] >= x & TPTororoCCA[,24] >= y) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,17] >= x & TPTororoREC[,18] >= y) / nrow(TPTororoREC) * 100
TR_p_25 <- sum(TPTororoREC[,19] >= x & TPTororoREC[,20] >= y) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,21] >= x & TPTororoREC[,22] >= y) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,23] >= x & TPTororoREC[,24] >= y) / nrow(TPTororoREC) * 100

d3_TPP <- data.frame(c(MC_p_2, MC_p_25, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_25, MR_p_3, MR_p_4),
                     c(TC_p_2, TC_p_25, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_25, TR_p_3, TR_p_4),
                     c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d3_TPP) <- c("A", "B", "C", "D", "Treshold")

d3_TPP <- melt(d3_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d3_o <- ggplot(d3_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "2.5 G score threshold" = "grey33",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "2.5 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "3 days sampling", 
       y = "Percentage achieving TPP", 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))+
  theme(legend.title = element_text(size = 17),  # Adjust size for legend title
        legend.text = element_text(size = 15))

# Combined plot
library(cowplot)
library(gridExtra)
library(scales)

legend <- cowplot::get_legend(d3_o)

empty_plot <- ggplot() + 
  theme_void() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )

combined_plot <- plot_grid(
  d1_o + theme(legend.position = "none",
               plot.title = element_text(size = rel(1.5)),     # Increase plot title size
               axis.title.x = element_text(size = rel(1.5)),   # Increase x axis title size
               axis.title.y = element_text(size = rel(1.5)),   # Increase y axis title size
               axis.text.x = element_text(size = rel(1.5)),    # Increase x axis text size
               axis.text.y = element_text(size = rel(1.5))),
  d2_o + theme(legend.position = "none",
               plot.title = element_text(size = rel(1.5)),     # Increase plot title size
               axis.title.x = element_text(size = rel(1.5)),   # Increase x axis title size
               axis.title.y = element_text(size = rel(1.5)),   # Increase y axis title size
               axis.text.x = element_text(size = rel(1.5)),    # Increase x axis text size
               axis.text.y = element_text(size = rel(1.5))),
  d3_o + theme(legend.position = "none",
               plot.title = element_text(size = rel(1.5)),     # Increase plot title size
               axis.title.x = element_text(size = rel(1.5)),   # Increase x axis title size
               axis.title.y = element_text(size = rel(1.5)),   # Increase y axis title size
               axis.text.x = element_text(size = rel(1.5)),    # Increase x axis text size
               axis.text.y = element_text(size = rel(1.5))),
  legend,
  d1_m + theme(legend.position = "none",
               plot.title = element_text(size = rel(1.5)),     # Increase plot title size
               axis.title.x = element_text(size = rel(1.5)),   # Increase x axis title size
               axis.title.y = element_text(size = rel(1.5)),   # Increase y axis title size
               axis.text.x = element_text(size = rel(1.5)),    # Increase x axis text size
               axis.text.y = element_text(size = rel(1.5))),
  d2_m + theme(legend.position = "none",
               plot.title = element_text(size = rel(1.5)),     # Increase plot title size
               axis.title.x = element_text(size = rel(1.5)),   # Increase x axis title size
               axis.title.y = element_text(size = rel(1.5)),   # Increase y axis title size
               axis.text.x = element_text(size = rel(1.5)),    # Increase x axis text size
               axis.text.y = element_text(size = rel(1.5))),
  d3_m + theme(legend.position = "none",
               plot.title = element_text(size = rel(1.5)),     # Increase plot title size
               axis.title.x = element_text(size = rel(1.5)),   # Increase x axis title size
               axis.title.y = element_text(size = rel(1.5)),   # Increase y axis title size
               axis.text.x = element_text(size = rel(1.5)),    # Increase x axis text size
               axis.text.y = element_text(size = rel(1.5))),
  empty_plot,
  ncol = 4, 
  nrow = 2,
  align = "h",
  rel_widths = c(1, 1, 1, 0.5)
)
png("TTP.png", 1300, 500)
combined_plot
dev.off()

P_m <- as.data.frame(c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]) * 100
colnames(P_m) <- "variable"
P_t <- as.data.frame(c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]) * 100
colnames(P_t) <- "variable"

pm <- ggplot(P_m, aes(x = variable)) + 
  geom_histogram(binwidth = 1, fill = "grey45", color = "black") +
  ggtitle("A") +
  xlab("Prevalence") +
  ylab("Frequency") +
  scale_x_continuous(labels = label_percent(scale = 1)) +
  theme_minimal()

pt <- ggplot(P_t, aes(x = variable)) + 
  geom_histogram(binwidth = 1, fill = "grey45", color = "black") +
  ggtitle("B") +
  xlab("Prevalence") +
  ylab("Frequency") +
  scale_x_continuous(labels = label_percent(scale = 1)) +
  theme_minimal()

combined_plot <- plot_grid(
  pm+ theme(plot.title = element_text(size = rel(1.3)),     # Increase plot title size
            axis.title.x = element_text(size = rel(1.2)),   # Increase x axis title size
            axis.title.y = element_text(size = rel(1.2)),   # Increase y axis title size
            axis.text.x = element_text(size = rel(1.3)),    # Increase x axis text size
            axis.text.y = element_text(size = rel(1.3)),
            legend.title = element_text(size = 12),  # Adjust size for legend title
            legend.text = element_text(size = 10)),
  pt+ theme(plot.title = element_text(size = rel(1.3)),     # Increase plot title size
            axis.title.x = element_text(size = rel(1.2)),   # Increase x axis title size
            axis.title.y = element_text(size = rel(1.2)),   # Increase y axis title size
            axis.text.x = element_text(size = rel(1.3)),    # Increase x axis text size
            axis.text.y = element_text(size = rel(1.3)),
            legend.title = element_text(size = 12),  # Adjust size for legend title
            legend.text = element_text(size = 10)),
  ncol = 2, 
  nrow = 1,
  align = "h"
)
png("Predicted Prevalence.png", 800, 340)
combined_plot
dev.off()
