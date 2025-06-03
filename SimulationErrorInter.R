###### Simulations to evaluate the error #####
## By M Entezami + JM Prada
rm(list = ls())

## Load library and set workspace
library(rstudioapi)
library(runjags)
library(reshape2)
library(ggplot2)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

## Load data
load("Results_inter_t.RData")
Results2 <- Results
load("Results_inter_m.RData")
source("SimulationErrorFunctionsInter_ME.R")

## Set draws
set.seed(100)
ndraws <- 1000
IDs <- sample.int(10000,ndraws)

#### Error plot ####
#### Mayuge
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]
# Extract the three tauCCA_intra columns
tauCCA_intra_values <- rbind(Results$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                             Results$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
# Randomly select one of the three columns for each ID
random_selection_CCA <- apply(tauCCA_intra_values[IDs, ], 1, function(x) sample(x, 1))
intraCCAsd <- sqrt(1 / random_selection_CCA)
tauRECCCA_intra_values <- rbind(Results$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                                Results$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
# Randomly select one of the three columns for each ID
random_selection_RECCCA <- apply(tauRECCCA_intra_values[IDs, ], 1, function(x) sample(x, 1))
intraRECCCAsd <- sqrt(1 / random_selection_RECCCA)
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
# Extract the three tauCCA_intra columns
tauCCA_intra_values <- rbind(Results2$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                             Results2$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
# Randomly select one of the three columns for each ID
random_selection_CCA <- apply(tauCCA_intra_values[IDs, ], 1, function(x) sample(x, 1))
intraCCAsd <- sqrt(1 / random_selection_CCA)
tauRECCCA_intra_values <- rbind(Results2$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                                Results2$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
# Randomly select one of the three columns for each ID
random_selection_RECCCA <- apply(tauRECCCA_intra_values[IDs, ], 1, function(x) sample(x, 1))
intraRECCCAsd <- sqrt(1 / random_selection_RECCCA)
## Simulate prevalence estimates for Tororo
PrevEstimatesTororo <- sapply(1:length(IDs),function(x){modelData(1000,x)})
## Calculate squared error
sderL <- apply((P-PrevEstimatesTororo)^2,1,quantile,c(.025,.5,.975))


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


#### ROC plot ####
#### Mayuge Sens/Spec
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
K <- c(Results$mcmc[[1]][,"k"],Results$mcmc[[2]][,"k"])[IDs]
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]
# Extract the three tauCCA_intra columns
intraCCAsd <- rbind(Results$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                    Results$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
intraCCAsd <- sqrt(1 / intraCCAsd)
intraRECCCAsd <- rbind(Results$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                       Results$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
intraRECCCAsd <- sqrt(1 / intraRECCCAsd)
multim1_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[1]"],Results$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[2]"],Results$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[1]"],Results$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[2]"],Results$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

## Simulate Sensitivity and Specificity for Mayuge 
SensSpecEstimatesMayugeCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
SensSpecEstimatesMayugeREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(1000,x)})
SensSpecEstimatesMayugeKK <- sapply(1:length(IDs),function(x){modelSensSpecKK(1000,x)})
## Calculate quantiles
QsenspcMayugeCCA <- apply(SensSpecEstimatesMayugeCCA,1,quantile,c(.025,.5,.975))
QsenspcMayugeREC <- apply(SensSpecEstimatesMayugeREC,1,quantile,c(.025,.5,.975))
QsenspcMayugeKK <- apply(SensSpecEstimatesMayugeKK,1,quantile,c(.025,.5,.975))

#### Tororo
interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauCCA_inter"],Results2$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
K <- c(Results2$mcmc[[1]][,"k"],Results2$mcmc[[2]][,"k"])[IDs]
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]
# Extract the three tauCCA_intra columns
intraCCAsd <- rbind(Results2$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                    Results2$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
intraCCAsd <- sqrt(1 / intraCCAsd)
intraRECCCAsd <- rbind(Results2$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                       Results2$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
intraRECCCAsd <- sqrt(1 / intraRECCCAsd)
multim1_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[1]"],Results2$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[2]"],Results2$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[1]"],Results2$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[2]"],Results2$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

## Simulate Sensitivity and Specificity for Tororo
SensSpecEstimatesTororoCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
SensSpecEstimatesTororoREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(1000,x)})
SensSpecEstimatesTororoKK <- sapply(1:length(IDs),function(x){modelSensSpecKK(1000,x)})
## Calculate quantiles
QsenspcTororoCCA <- apply(SensSpecEstimatesTororoCCA,1,quantile,c(.025,.5,.975))
QsenspcTororoREC <- apply(SensSpecEstimatesTororoREC,1,quantile,c(.025,.5,.975))
QsenspcTororoKK <- apply(SensSpecEstimatesTororoKK,1,quantile,c(.025,.5,.975))
### Plot ROC curve ###

png("Figure 2 - ROCcurve.png",width=800,height=500)
par(font=2, cex.axis=1, lwd=2, mar=c(4, 4, 2, 1), mgp=c(3, 0.4, 0),  oma=c(0, 2, 1, 0)) # Adjust margins
par(mfcol=c(2,2))
## Tororo
plot(NA,xlim = c(0,.6),ylim = c(.7,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcTororoCCA[2,c(2,4,8)],QsenspcTororoCCA[2,c(1,3,7)], pch=19, col="#f78f0b")
lines(1-QsenspcTororoCCA[2,c(2,4,8)],QsenspcTororoCCA[2,c(1,3,7)], col="#f78f0b")
arrows(x0=1-QsenspcTororoCCA[2,c(2,4,8)],y0=QsenspcTororoCCA[1,c(1,3,7)],
       y1=QsenspcTororoCCA[3,c(1,3,7)], angle=90, code = 3,length = .05, col="#f78f0b")
arrows(x0=1-QsenspcTororoCCA[1,c(2,4,8)],y0=QsenspcTororoCCA[2,c(1,3,7)],
       x1=1-QsenspcTororoCCA[3,c(2,4,8)], angle=90, code = 3,length = .05, col="#f78f0b")

labels <- c("G2", "G3", "G4")
x_coords <- 1 - QsenspcTororoCCA[2, c(2, 4, 8)]
y_coords <- QsenspcTororoCCA[2, c(1, 3, 7)]
text(x_coords, y_coords - 0.03, labels, pos = 4, offset = 0.4)
## 2 Days
points(1-QsenspcTororoCCA[2,c(10,12,16)],QsenspcTororoCCA[2,c(9,11,15)], pch=19, col="#59ce48")
lines(1-QsenspcTororoCCA[2,c(10,12,16)],QsenspcTororoCCA[2,c(9,11,15)], col="#59ce48")
arrows(x0=1-QsenspcTororoCCA[2,c(10,12,16)],y0=QsenspcTororoCCA[1,c(9,11,15)],
       y1=QsenspcTororoCCA[3,c(9,11,15)], angle=90, code = 3,length = .05, col="#59ce48")
arrows(x0=1-QsenspcTororoCCA[1,c(10,12,16)],y0=QsenspcTororoCCA[2,c(9,11,15)],
       x1=1-QsenspcTororoCCA[3,c(10,12,16)], angle=90, code = 3,length = .05, col="#59ce48")

## 3 Days
points(1-QsenspcTororoCCA[2,c(18,20,24)],QsenspcTororoCCA[2,c(17,19,23)], pch=19, col="#8628af")
lines(1-QsenspcTororoCCA[2,c(18,20,24)],QsenspcTororoCCA[2,c(17,19,23)], col="#8628af")
arrows(x0=1-QsenspcTororoCCA[2,c(18,20,24)],y0=QsenspcTororoCCA[1,c(17,19,23)],
       y1=QsenspcTororoCCA[3,c(17,19,23)], angle=90, code = 3,length = .05, col="#8628af")
arrows(x0=1-QsenspcTororoCCA[1,c(18,20,24)],y0=QsenspcTororoCCA[2,c(17,19,23)],
       x1=1-QsenspcTororoCCA[3,c(18,20,24)], angle=90, code = 3,length = .05, col="#8628af")

axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)

mtext("Sensitivity",side=2,cex=1.25,line=2)
mtext("1 - Specificity",side=1,cex=1.25,line=2)

# Add extra labels
mtext("Tororo", side=3, outer=TRUE, at=0.25, las=1, cex=1.5, col="black", line = -0.35)
mtext("Mayuge", side=3, outer=TRUE, at=0.75, las=1, cex=1.5, col="black", line = -0.35)
mtext("POC-CCA", side=2, outer=TRUE, at=0.77, las=3, cex=1.5, col="black", line = -0.35)
mtext("POC-CCA3", side=2, outer=TRUE, at=0.27, las=3, cex=1.5, col="black", line = -0.35)


## Tororo REC
plot(NA,xlim = c(0,.6),ylim = c(.7,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcTororoREC[2,c(2,4,8)],QsenspcTororoREC[2,c(1,3,7)], pch=19, col="#f78f0b")
lines(1-QsenspcTororoREC[2,c(2,4,8)],QsenspcTororoREC[2,c(1,3,7)], col="#f78f0b")
arrows(x0=1-QsenspcTororoREC[2,c(2,4,8)],y0=QsenspcTororoREC[1,c(1,3,7)],
       y1=QsenspcTororoREC[3,c(1,3,7)], angle=90, code = 3,length = .05, col="#f78f0b")
arrows(x0=1-QsenspcTororoREC[1,c(2,4,8)],y0=QsenspcTororoREC[2,c(1,3,7)],
       x1=1-QsenspcTororoREC[3,c(2,4,8)], angle=90, code = 3,length = .05, col="#f78f0b")

## 2 Days
points(1-QsenspcTororoREC[2,c(10,12,16)],QsenspcTororoREC[2,c(9,11,15)], pch=19, col="#59ce48")
lines(1-QsenspcTororoREC[2,c(10,12,16)],QsenspcTororoREC[2,c(9,11,15)], col="#59ce48")
arrows(x0=1-QsenspcTororoREC[2,c(10,12,16)],y0=QsenspcTororoREC[1,c(9,11,15)],
       y1=QsenspcTororoREC[3,c(9,11,15)], angle=90, code = 3,length = .05, col="#59ce48")
arrows(x0=1-QsenspcTororoREC[1,c(10,12,16)],y0=QsenspcTororoREC[2,c(9,11,15)],
       x1=1-QsenspcTororoREC[3,c(10,12,16)], angle=90, code = 3,length = .05, col="#59ce48")

## 3 Days
points(1-QsenspcTororoREC[2,c(18,20,24)],QsenspcTororoREC[2,c(17,19,23)], pch=19, col="#8628af")
lines(1-QsenspcTororoREC[2,c(18,20,24)],QsenspcTororoREC[2,c(17,19,23)], col="#8628af")
arrows(x0=1-QsenspcTororoREC[2,c(18,20,24)],y0=QsenspcTororoREC[1,c(17,19,23)],
       y1=QsenspcTororoREC[3,c(17,19,23)], angle=90, code = 3,length = .05, col="#8628af")
arrows(x0=1-QsenspcTororoREC[1,c(18,20,24)],y0=QsenspcTororoREC[2,c(17,19,23)],
       x1=1-QsenspcTororoREC[3,c(18,20,24)], angle=90, code = 3,length = .05, col="#8628af")

axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)

mtext("Sensitivity",side=2,cex=1.25,line=2)
mtext("1 - Specificity",side=1,cex=1.25,line=2)

## Mayuge CCA
plot(NA,xlim = c(0,.6),ylim = c(.7,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcMayugeCCA[2,c(2,4,8)],QsenspcMayugeCCA[2,c(1,3,7)], pch=19, col="#f78f0b")
lines(1-QsenspcMayugeCCA[2,c(2,4,8)],QsenspcMayugeCCA[2,c(1,3,7)], col="#f78f0b")
arrows(x0=1-QsenspcMayugeCCA[2,c(2,4,8)],y0=QsenspcMayugeCCA[1,c(1,3,7)],
       y1=QsenspcMayugeCCA[3,c(1,3,7)], angle=90, code = 3,length = .05, col="#f78f0b")
arrows(x0=1-QsenspcMayugeCCA[1,c(2,4,8)],y0=QsenspcMayugeCCA[2,c(1,3,7)],
       x1=1-QsenspcMayugeCCA[3,c(2,4,8)], angle=90, code = 3,length = .05, col="#f78f0b")

## 2 Days
points(1-QsenspcMayugeCCA[2,c(10,12,16)],QsenspcMayugeCCA[2,c(9,11,15)], pch=19, col="#59ce48")
lines(1-QsenspcMayugeCCA[2,c(10,12,16)],QsenspcMayugeCCA[2,c(9,11,15)], col="#59ce48")
arrows(x0=1-QsenspcMayugeCCA[2,c(10,12,16)],y0=QsenspcMayugeCCA[1,c(9,11,15)],
       y1=QsenspcMayugeCCA[3,c(9,11,15)], angle=90, code = 3,length = .05, col="#59ce48")
arrows(x0=1-QsenspcMayugeCCA[1,c(10,12,16)],y0=QsenspcMayugeCCA[2,c(9,11,15)],
       x1=1-QsenspcMayugeCCA[3,c(10,12,16)], angle=90, code = 3,length = .05, col="#59ce48")

## 3 Days
points(1-QsenspcMayugeCCA[2,c(18,20,24)],QsenspcMayugeCCA[2,c(17,19,23)], pch=19, col="#8628af")
lines(1-QsenspcMayugeCCA[2,c(18,20,24)],QsenspcMayugeCCA[2,c(17,19,23)], col="#8628af")
arrows(x0=1-QsenspcMayugeCCA[2,c(18,20,24)],y0=QsenspcMayugeCCA[1,c(17,19,23)],
       y1=QsenspcMayugeCCA[3,c(17,19,23)], angle=90, code = 3,length = .05, col="#8628af")
arrows(x0=1-QsenspcMayugeCCA[1,c(18,20,24)],y0=QsenspcMayugeCCA[2,c(17,19,23)],
       x1=1-QsenspcMayugeCCA[3,c(18,20,24)], angle=90, code = 3,length = .05, col="#8628af")

axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)

mtext("Sensitivity",side=2,cex=1.25,line=2)
mtext("1 - Specificity",side=1,cex=1.25,line=2)

## Mayuge REC
plot(NA,xlim = c(0,.6),ylim = c(.7,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcMayugeREC[2,c(2,4,8)],QsenspcMayugeREC[2,c(1,3,7)], pch=19, col="#f78f0b")
lines(1-QsenspcMayugeREC[2,c(2,4,8)],QsenspcMayugeREC[2,c(1,3,7)], col="#f78f0b")
arrows(x0=1-QsenspcMayugeREC[2,c(2,4,8)],y0=QsenspcMayugeREC[1,c(1,3,7)],
       y1=QsenspcMayugeREC[3,c(1,3,7)], angle=90, code = 3,length = .05, col="#f78f0b")
arrows(x0=1-QsenspcMayugeREC[1,c(2,4,8)],y0=QsenspcMayugeREC[2,c(1,3,7)],
       x1=1-QsenspcMayugeREC[3,c(2,4,8)], angle=90, code = 3,length = .05, col="#f78f0b")

## 2 Days
points(1-QsenspcMayugeREC[2,c(10,12,16)],QsenspcMayugeREC[2,c(9,11,15)], pch=19, col="#59ce48")
lines(1-QsenspcMayugeREC[2,c(10,12,16)],QsenspcMayugeREC[2,c(9,11,15)], col="#59ce48")
arrows(x0=1-QsenspcMayugeREC[2,c(10,12,16)],y0=QsenspcMayugeREC[1,c(9,11,15)],
       y1=QsenspcMayugeREC[3,c(9,11,15)], angle=90, code = 3,length = .05, col="#59ce48")
arrows(x0=1-QsenspcMayugeREC[1,c(10,12,16)],y0=QsenspcMayugeREC[2,c(9,11,15)],
       x1=1-QsenspcMayugeREC[3,c(10,12,16)], angle=90, code = 3,length = .05, col="#59ce48")

## 3 Days
points(1-QsenspcMayugeREC[2,c(18,20,24)],QsenspcMayugeREC[2,c(17,19,23)], pch=19, col="#8628af")
lines(1-QsenspcMayugeREC[2,c(18,20,24)],QsenspcMayugeREC[2,c(17,19,23)], col="#8628af")
arrows(x0=1-QsenspcMayugeREC[2,c(18,20,24)],y0=QsenspcMayugeREC[1,c(17,19,23)],
       y1=QsenspcMayugeREC[3,c(17,19,23)], angle=90, code = 3,length = .05, col="#8628af")
arrows(x0=1-QsenspcMayugeREC[1,c(18,20,24)],y0=QsenspcMayugeREC[2,c(17,19,23)],
       x1=1-QsenspcMayugeREC[3,c(18,20,24)], angle=90, code = 3,length = .05, col="#8628af")

axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5)
mtext("Sensitivity",side=2,cex=1.25,line=2)
mtext("1 - Specificity",side=1,cex=1.25,line=2)

legend("bottomright",c("1 Day","2 Days","3 Days"),
       col=c("#f78f0b","#59ce48","#8628af"),
       bty='n',cex=1.25,lty=1)

dev.off()

# Define thresholds and time indices
thresholds <- c(2, 2.5, 3, 4)
time_cols <- list(
  "1 Day" = list(sens = c(1,3,5,7), spec = c(2,4,6,8)),
  "2 Days" = list(sens = c(9,11,13,15), spec = c(10,12,14,16)),
  "3 Days" = list(sens = c(17,19,21,23), spec = c(18,20,22,24))
)

# List of data arrays with corresponding settings and diagnostics
data_list <- list(
  list(array = QsenspcTororoREC, Setting = "Tororo", Diagnostic = "POC-CCA3"),
  list(array = QsenspcTororoCCA, Setting = "Tororo", Diagnostic = "POC-CCA"),
  list(array = QsenspcMayugeREC, Setting = "Mayuge", Diagnostic = "POC-CCA3"),
  list(array = QsenspcMayugeCCA, Setting = "Mayuge", Diagnostic = "POC-CCA")
)

# Initialise an empty data frame
result <- data.frame()

# Loop through each combination to extract metrics
for(d in data_list) {
  arr <- d$array
  for(time in names(time_cols)) {
    sens_idx <- time_cols[[time]]$sens
    spec_idx <- time_cols[[time]]$spec
    for(i in 1:4) {
      thresh <- thresholds[i]
      # Sensitivity values (mean, lower, upper)
      sens_mean  <- arr[2, sens_idx[i]]
      sens_lower <- arr[1, sens_idx[i]]
      sens_upper <- arr[3, sens_idx[i]]
      # Specificity values (note the conversion: specificity = 1 - false positive rate)
      spec_mean  <- arr[2, spec_idx[i]]
      spec_lower <- arr[1, spec_idx[i]]  # lower bound for specificity
      spec_upper <- arr[3, spec_idx[i]]  # upper bound for specificity
      
      result <- rbind(result, data.frame(
        Setting             = d$Setting,
        Diagnostic          = d$Diagnostic,
        Time                = time,
        Threshold           = thresh,
        Sensitivity_mean    = sens_mean,
        Sensitivity_lower   = sens_lower,
        Sensitivity_upper   = sens_upper,
        Specificity_mean    = spec_mean,
        Specificity_lower   = spec_lower,
        Specificity_upper   = spec_upper
      ))
    }
  }
}

# Format the sensitivity and specificity as mean% (lower% - upper%) with two decimal places
result$Sensitivity_fmt <- sprintf("%.1f%% (%.1f%% - %.1f%%)",
                                  result$Sensitivity_mean * 100,
                                  result$Sensitivity_lower * 100,
                                  result$Sensitivity_upper * 100)
result$Specificity_fmt <- sprintf("%.1f%% (%.1f%% - %.1f%%)",
                                  result$Specificity_mean * 100,
                                  result$Specificity_lower * 100,
                                  result$Specificity_upper * 100)

# Select and rename columns for the final table
final_table <- result[, c("Setting", "Diagnostic", "Time", "Threshold",
                          "Sensitivity_fmt", "Specificity_fmt")]
names(final_table) <- c("Setting", "Diagnostic", "Time", "Threshold",
                        "Sensitivity", "Specificity")

# Print the final table
print(final_table)


write.csv(final_table, file = "SeSf.csv", , row.names = F)
# Create a data frame for Tororo KK results
tororo_table <- data.frame(
  Setting = "Tororo",
  Day = c("1 Day", "2 Days", "3 Days"),
  Sensitivity = c(
    sprintf("%.1f%% (%.1f%% - %.1f%%)",
            QsenspcTororoKK[2, 1] * 100,
            QsenspcTororoKK[1, 1] * 100,
            QsenspcTororoKK[3, 1] * 100),
    sprintf("%.1f%% (%.1f%% - %.1f%%)",
            QsenspcTororoKK[2, 3] * 100,
            QsenspcTororoKK[1, 3] * 100,
            QsenspcTororoKK[3, 3] * 100),
    sprintf("%.1f%% (%.1f%% - %.1f%%)",
            QsenspcTororoKK[2, 5] * 100,
            QsenspcTororoKK[1, 5] * 100,
            QsenspcTororoKK[3, 5] * 100)
  ),
  stringsAsFactors = FALSE
)

# Create a data frame for Mayuge KK results
mayuge_table <- data.frame(
  Setting = "Mayuge",
  Day = c("1 Day", "2 Days", "3 Days"),
  Sensitivity = c(
    sprintf("%.1f%% (%.1f%% - %.1f%%)",
            QsenspcMayugeKK[2, 1] * 100,
            QsenspcMayugeKK[1, 1] * 100,
            QsenspcMayugeKK[3, 1] * 100),
    sprintf("%.1f%% (%.1f%% - %.1f%%)",
            QsenspcMayugeKK[2, 3] * 100,
            QsenspcMayugeKK[1, 3] * 100,
            QsenspcMayugeKK[3, 3] * 100),
    sprintf("%.1f%% (%.1f%% - %.1f%%)",
            QsenspcMayugeKK[2, 5] * 100,
            QsenspcMayugeKK[1, 5] * 100,
            QsenspcMayugeKK[3, 5] * 100)
  ),
  stringsAsFactors = FALSE
)

# Combine the two tables
combined_table <- rbind(tororo_table, mayuge_table)

write.csv(combined_table, file = "KKsesp.csv", , row.names = F)
#### Auc table ####
#### Mayuge Sens/Spec
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
K <- c(Results$mcmc[[1]][,"k"],Results$mcmc[[2]][,"k"])[IDs]
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]
# Extract the three tauCCA_intra columns
intraCCAsd <- rbind(Results$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                    Results$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
intraCCAsd <- sqrt(1 / intraCCAsd)
intraRECCCAsd <- rbind(Results$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                       Results$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
intraRECCCAsd <- sqrt(1 / intraRECCCAsd)
multim1_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[1]"],Results$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[2]"],Results$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[1]"],Results$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[2]"],Results$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

AUCMayugeCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
AUCMayugeREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(1000,x)})

interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauCCA_inter"],Results2$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
K <- c(Results2$mcmc[[1]][,"k"],Results2$mcmc[[2]][,"k"])[IDs]
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]
# Extract the three tauCCA_intra columns
intraCCAsd <- rbind(Results2$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                    Results2$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
intraCCAsd <- sqrt(1 / intraCCAsd)
intraRECCCAsd <- rbind(Results2$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                       Results2$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
intraRECCCAsd <- sqrt(1 / intraRECCCAsd)
multim1_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[1]"],Results2$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[2]"],Results2$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[1]"],Results2$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[2]"],Results2$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

AUCTororoCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
AUCTororoREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(1000,x)})

#### Mayuge AUC
### CCA 
## 1 Day 
# Initialize an empty vector to store AUC values
AUCMayugeCCA <- t(AUCMayugeCCA)
auc_values <- numeric(nrow(AUCMayugeCCA))
M_CCA_1 <- c()
# Loop through each repeat (row)
for (j in 1:nrow(AUCMayugeCCA)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCMayugeCCA[j, c(1, 3, 7)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeCCA[j, c(2, 4, 8)], 1)
  
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
  sensitivity = c(0, AUCMayugeCCA[j, c(9,11,15)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeCCA[j, c(10,12,16)], 1)
  
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
  sensitivity = c(0, AUCMayugeCCA[j, c(17,19,23)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeCCA[j, c(18,20,24)], 1)
  
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
  sensitivity = c(0, AUCMayugeREC[j, c(1, 3, 7)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeREC[j, c(2, 4, 8)], 1)
  
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
  sensitivity = c(0, AUCMayugeREC[j, c(9,11,15)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeREC[j, c(10,12,16)], 1)
  
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
  sensitivity = c(0, AUCMayugeREC[j, c(17,19,23)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCMayugeREC[j, c(18,20,24)], 1)
  
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

#### Tororo AUC
### CCA 
## 1 Day 
# Initialize an empty vector to store AUC values
AUCTororoCCA <- t(AUCTororoCCA)
auc_values <- numeric(nrow(AUCTororoCCA))
T_CCA_1 <- c()
# Loop through each repeat (row)
for (j in 1:nrow(AUCTororoCCA)) {
  # Extract sensitivity and 1-specificity values for this repeat
  sensitivity = c(0, AUCTororoCCA[j, c(1, 3, 7)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoCCA[j, c(2, 4, 8)], 1)
  
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
  sensitivity = c(0, AUCTororoCCA[j, c(9,11,15)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoCCA[j, c(10,12,16)], 1)
  
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
  sensitivity = c(0, AUCTororoCCA[j, c(17,19,23)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoCCA[j, c(18,20,24)], 1)
  
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
  sensitivity = c(0, AUCTororoREC[j, c(1, 3, 7)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoREC[j, c(2, 4, 8)], 1)
  
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
  sensitivity = c(0, AUCTororoREC[j, c(9,11,15)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoREC[j, c(10,12,16)], 1)
  
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
  sensitivity = c(0, AUCTororoREC[j, c(17,19,23)], 1)  # Adding points (0,0) and (1,1)
  one_minus_specificity = c(0, 1 - AUCTororoREC[j, c(18,20,24)], 1)
  
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

#Significance testing


T_all_CCA <- c(T_CCA_1, T_CCA_2, T_CCA_3)
T_all_REC <- c(T_REC_1, T_REC_2, T_REC_3)

wilcox.test(T_all_CCA, T_all_REC, paired = TRUE)

M_all_CCA <- c(M_CCA_1, M_CCA_2, M_CCA_3)
M_all_REC <- c(M_REC_1, M_REC_2, M_REC_3)

wilcox.test(M_all_CCA, M_all_REC, paired = TRUE)

# Day 1 comparison
wilcox.test(T_CCA_1, T_REC_1, paired = TRUE, alternative = "less")
wilcox.test(T_CCA_2, T_REC_2, paired = TRUE, alternative = "less")
wilcox.test(T_CCA_3, T_REC_3, paired = TRUE, alternative = "less")
wilcox.test(M_CCA_1, M_REC_1, paired = TRUE, alternative = "less")
wilcox.test(M_CCA_2, M_REC_2, paired = TRUE, alternative = "less")
wilcox.test(M_CCA_3, M_REC_3, paired = TRUE, alternative = "less")
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

AUC_table <- data.frame(
  `Number of sampling days` = c("1 day", "2 days", "3 days"),
  `Tororo POC-CCA` = c(
    paste0(round(QAUCTororoCCA["50%", "T_CCA_1"] * 100, 1), "% (",
           round(QAUCTororoCCA["2.5%", "T_CCA_1"] * 100, 1), "% - ",
           round(QAUCTororoCCA["97.5%", "T_CCA_1"] * 100, 1), "%)"),
    paste0(round(QAUCTororoCCA["50%", "T_CCA_2"] * 100, 1), "% (",
           round(QAUCTororoCCA["2.5%", "T_CCA_2"] * 100, 1), "% - ",
           round(QAUCTororoCCA["97.5%", "T_CCA_2"] * 100, 1), "%)"),
    paste0(round(QAUCTororoCCA["50%", "T_CCA_3"] * 100, 1), "% (",
           round(QAUCTororoCCA["2.5%", "T_CCA_3"] * 100, 1), "% - ",
           round(QAUCTororoCCA["97.5%", "T_CCA_3"] * 100, 1), "%)")
  ),
  `Tororo POC-CCA3` = c(
    paste0(round(QAUCTororoREC["50%", "T_REC_1"] * 100, 1), "% (",
           round(QAUCTororoREC["2.5%", "T_REC_1"] * 100, 1), "% - ",
           round(QAUCTororoREC["97.5%", "T_REC_1"] * 100, 1), "%)"),
    paste0(round(QAUCTororoREC["50%", "T_REC_2"] * 100, 1), "% (",
           round(QAUCTororoREC["2.5%", "T_REC_2"] * 100, 1), "% - ",
           round(QAUCTororoREC["97.5%", "T_REC_2"] * 100, 1), "%)"),
    paste0(round(QAUCTororoREC["50%", "T_REC_3"] * 100, 1), "% (",
           round(QAUCTororoREC["2.5%", "T_REC_3"] * 100, 1), "% - ",
           round(QAUCTororoREC["97.5%", "T_REC_3"] * 100, 1), "%)")
  ),
  `Mayuge POC-CCA` = c(
    paste0(round(QAUCMayugeCCA["50%", "M_CCA_1"] * 100, 1), "% (",
           round(QAUCMayugeCCA["2.5%", "M_CCA_1"] * 100, 1), "% - ",
           round(QAUCMayugeCCA["97.5%", "M_CCA_1"] * 100, 1), "%)"),
    paste0(round(QAUCMayugeCCA["50%", "M_CCA_2"] * 100, 1), "% (",
           round(QAUCMayugeCCA["2.5%", "M_CCA_2"] * 100, 1), "% - ",
           round(QAUCMayugeCCA["97.5%", "M_CCA_2"] * 100, 1), "%)"),
    paste0(round(QAUCMayugeCCA["50%", "M_CCA_3"] * 100, 1), "% (",
           round(QAUCMayugeCCA["2.5%", "M_CCA_3"] * 100, 1), "% - ",
           round(QAUCMayugeCCA["97.5%", "M_CCA_3"] * 100, 1), "%)")
  ),
  `Mayuge POC-CCA3` = c(
    paste0(round(QAUCMayugeREC["50%", "M_REC_1"] * 100, 1), "% (",
           round(QAUCMayugeREC["2.5%", "M_REC_1"] * 100, 1), "% - ",
           round(QAUCMayugeREC["97.5%", "M_REC_1"] * 100, 1), "%)"),
    paste0(round(QAUCMayugeREC["50%", "M_REC_2"] * 100, 1), "% (",
           round(QAUCMayugeREC["2.5%", "M_REC_2"] * 100, 1), "% - ",
           round(QAUCMayugeREC["97.5%", "M_REC_2"] * 100, 1), "%)"),
    paste0(round(QAUCMayugeREC["50%", "M_REC_3"] * 100, 1), "% (",
           round(QAUCMayugeREC["2.5%", "M_REC_3"] * 100, 1), "% - ",
           round(QAUCMayugeREC["97.5%", "M_REC_3"] * 100, 1), "%)")
  )
)

write.csv(AUC_table, file = "AUC_tabletst.csv", row.names = FALSE)
### Target profile ####
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauCCA_inter"],Results$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tauRECCCA_inter"],Results$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
K <- c(Results$mcmc[[1]][,"k"],Results$mcmc[[2]][,"k"])[IDs]
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]
# Extract the three tauCCA_intra columns
intraCCAsd <- rbind(Results$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                    Results$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
intraCCAsd <- sqrt(1 / intraCCAsd)
intraRECCCAsd <- rbind(Results$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                       Results$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
intraRECCCAsd <- sqrt(1 / intraRECCCAsd)
multim1_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[1]"],Results$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results$mcmc[[1]][,"multiParamCCA[2]"],Results$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[1]"],Results$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results$mcmc[[1]][,"multiParamRECCCA[2]"],Results$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

TPMayugeCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(100,x)})
TPMayugeREC <- sapply(1:length(IDs),function(x){modelSensSpecREC(100,x)})
TPMayugeCCA <- t(TPMayugeCCA)
TPMayugeREC <- t(TPMayugeREC)

interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauCCA_inter"],Results2$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
K <- c(Results2$mcmc[[1]][,"k"],Results2$mcmc[[2]][,"k"])[IDs]
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]
# Extract the three tauCCA_intra columns
intraCCAsd <- rbind(Results2$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                    Results2$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
intraCCAsd <- sqrt(1 / intraCCAsd)
intraRECCCAsd <- rbind(Results2$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                       Results2$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
intraRECCCAsd <- sqrt(1 / intraRECCCAsd)
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
MC_p_2 <- sum(TPMayugeCCA[,1] >= x & TPMayugeCCA[,2] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,5] >= x & TPMayugeCCA[,6] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,7] >= x & TPMayugeCCA[,8] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,1] >= x & TPMayugeREC[,2] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,5] >= x & TPMayugeREC[,6] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,7] >= x & TPMayugeREC[,8] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,1] >= x & TPTororoCCA[,2] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,5] >= x & TPTororoCCA[,6] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,7] >= x & TPTororoCCA[,8] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,1] >= x & TPTororoREC[,2] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,5] >= x & TPTororoREC[,6] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,7] >= x & TPTororoREC[,8] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100

d1_TPP <- data.frame(c(TC_p_2, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_3, TR_p_4),
                     c(MC_p_2, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_3, MR_p_4),
                     c("2 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d1_TPP) <- c("A", "B", "C", "D", "Treshold")


d1_TPP <- melt(d1_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d1_m <- ggplot(d1_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = NULL, 
       y = "Percentage achieving\nminimum TPP", 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

MC_p_2 <- sum(TPMayugeCCA[,9] >= x & TPMayugeCCA[,10] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,13] >= x & TPMayugeCCA[,14] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,15] >= x & TPMayugeCCA[,16] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,9] >= x & TPMayugeREC[,10] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,13] >= x & TPMayugeREC[,14] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,15] >= x & TPMayugeREC[,16] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,9] >= x & TPTororoCCA[,10] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,13] >= x & TPTororoCCA[,14] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,15] >= x & TPTororoCCA[,16] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,9] >= x & TPTororoREC[,10] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,13] >= x & TPTororoREC[,14] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,15] >= x & TPTororoREC[,16] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100

d2_TPP <- data.frame(c(TC_p_2, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_3, TR_p_4),
                     c(MC_p_2, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_3, MR_p_4),
                     c("2 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d2_TPP) <- c("A", "B", "C", "D", "Treshold")

d2_TPP <- melt(d2_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d2_m <- ggplot(d2_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = NULL, 
       y = NULL, 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

MC_p_2 <- sum(TPMayugeCCA[,17] >= x & TPMayugeCCA[,18] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,21] >= x & TPMayugeCCA[,22] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,23] >= x & TPMayugeCCA[,24] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,17] >= x & TPMayugeREC[,18] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,21] >= x & TPMayugeREC[,22] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,23] >= x & TPMayugeREC[,24] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,17] >= x & TPTororoCCA[,18] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,21] >= x & TPTororoCCA[,22] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,23] >= x & TPTororoCCA[,24] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,17] >= x & TPTororoREC[,18] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,21] >= x & TPTororoREC[,22] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,23] >= x & TPTororoREC[,24] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100

d3_TPP <- data.frame(c(TC_p_2, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_3, TR_p_4),
                     c(MC_p_2, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_3, MR_p_4),
                     c("2 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d3_TPP) <- c("A", "B", "C", "D", "Treshold")

d3_TPP <- melt(d3_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d3_m <- ggplot(d3_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = NULL, 
       y = NULL, 
       x = "Setting") +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

# Condition values
x <- 0.75 # Your threshold for column1
y <- 0.965 # Your threshold for column2

# Calculate the percentage
MC_p_2 <- sum(TPMayugeCCA[,1] >= x & TPMayugeCCA[,2] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,5] >= x & TPMayugeCCA[,6] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,7] >= x & TPMayugeCCA[,8] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,1] >= x & TPMayugeREC[,2] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,5] >= x & TPMayugeREC[,6] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,7] >= x & TPMayugeREC[,8] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,1] >= x & TPTororoCCA[,2] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,5] >= x & TPTororoCCA[,6] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,7] >= x & TPTororoCCA[,8] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,1] >= x & TPTororoREC[,2] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,5] >= x & TPTororoREC[,6] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,7] >= x & TPTororoREC[,8] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100

d1_TPP <- data.frame(c(TC_p_2, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_3, TR_p_4),
                     c(MC_p_2, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_3, MR_p_4),
                     c("2 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d1_TPP) <- c("A", "B", "C", "D", "Treshold")


d1_TPP <- melt(d1_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d1_o <- ggplot(d1_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "1 days sampling", 
       y = "Percentage achieving\noptimal TPP", 
       x = NULL) +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

MC_p_2 <- sum(TPMayugeCCA[,9] >= x & TPMayugeCCA[,10] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,13] >= x & TPMayugeCCA[,14] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,15] >= x & TPMayugeCCA[,16] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,9] >= x & TPMayugeREC[,10] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,13] >= x & TPMayugeREC[,14] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,15] >= x & TPMayugeREC[,16] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,9] >= x & TPTororoCCA[,10] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,13] >= x & TPTororoCCA[,14] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,15] >= x & TPTororoCCA[,16] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,9] >= x & TPTororoREC[,10] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,13] >= x & TPTororoREC[,14] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,15] >= x & TPTororoREC[,16] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100

d2_TPP <- data.frame(c(TC_p_2, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_3, TR_p_4),
                     c(MC_p_2, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_3, MR_p_4),
                     c("2 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d2_TPP) <- c("A", "B", "C", "D", "Treshold")

d2_TPP <- melt(d2_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d2_o <- ggplot(d2_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "2 days sampling", 
       y = NULL, 
       x = NULL) +
  theme_bw()+
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0,100))

MC_p_2 <- sum(TPMayugeCCA[,17] >= x & TPMayugeCCA[,18] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_3 <- sum(TPMayugeCCA[,21] >= x & TPMayugeCCA[,22] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100
MC_p_4 <- sum(TPMayugeCCA[,23] >= x & TPMayugeCCA[,24] >= y, na.rm = TRUE) / nrow(TPMayugeCCA) * 100

MR_p_2 <- sum(TPMayugeREC[,17] >= x & TPMayugeREC[,18] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_3 <- sum(TPMayugeREC[,21] >= x & TPMayugeREC[,22] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100
MR_p_4 <- sum(TPMayugeREC[,23] >= x & TPMayugeREC[,24] >= y, na.rm = TRUE) / nrow(TPMayugeREC) * 100

TC_p_2 <- sum(TPTororoCCA[,17] >= x & TPTororoCCA[,18] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_3 <- sum(TPTororoCCA[,21] >= x & TPTororoCCA[,22] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100
TC_p_4 <- sum(TPTororoCCA[,23] >= x & TPTororoCCA[,24] >= y, na.rm = TRUE) / nrow(TPTororoCCA) * 100

TR_p_2 <- sum(TPTororoREC[,17] >= x & TPTororoREC[,18] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_3 <- sum(TPTororoREC[,21] >= x & TPTororoREC[,22] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100
TR_p_4 <- sum(TPTororoREC[,23] >= x & TPTororoREC[,24] >= y, na.rm = TRUE) / nrow(TPTororoREC) * 100

d3_TPP <- data.frame(c(TC_p_2, TC_p_3, TC_p_4),
                     c(TR_p_2, TR_p_3, TR_p_4),
                     c(MC_p_2, MC_p_3, MC_p_4),
                     c(MR_p_2, MR_p_3, MR_p_4),
                     c("2 G score threshold", "3 G score threshold", "4 G score threshold"))
colnames(d3_TPP) <- c("A", "B", "C", "D", "Treshold")

d3_TPP <- melt(d3_TPP, id.vars = "Treshold", variable.name = "Location", value.name = "Value")

d3_o <- ggplot(d3_TPP, aes(x = Location, y = Value, fill = Treshold)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.9) +  
  scale_fill_manual(values = c("2 G score threshold" = "grey0",
                               "3 G score threshold" = "grey66",
                               "4 G score threshold" = "grey100"),
                    labels = c("2 G score threshold", "3 G score threshold", "4 G score threshold")) +
  labs(title = "3 days sampling", 
       x = NULL,
       y = NULL) +
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
  rel_widths = c(1.2, 1, 1, 0.5)
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

#### Simulate 10% prev ####
library(MASS)
interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauCCA_inter"],Results2$mcmc[[2]][,"tauCCA_inter"])[IDs])
interRECCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tauRECCCA_inter"],Results2$mcmc[[2]][,"tauRECCCA_inter"])[IDs])
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
indices <- order(P)[1:100]
IDs <- sample(indices, 10000, replace = TRUE)
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]
K <- c(Results2$mcmc[[1]][,"k"],Results2$mcmc[[2]][,"k"])[IDs]
P <- rep(0.2, 10000)
# Extract the three tauCCA_intra columns
intraCCAsd <- rbind(Results2$mcmc[[1]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")],
                             Results2$mcmc[[2]][, c("tauCCA_batch[1]", "tauCCA_batch[2]", "tauCCA_batch[3]")])
intraCCAsd <- sqrt(1 / intraCCAsd)
intraRECCCAsd <- rbind(Results2$mcmc[[1]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")],
                                Results2$mcmc[[2]][, c("tauRECCCA_batch[1]", "tauRECCCA_batch[2]", "tauRECCCA_batch[3]")])
intraRECCCAsd <- sqrt(1 / intraRECCCAsd)
multim1_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[1]"],Results2$mcmc[[2]][,"multiParamCCA[1]"])[IDs]
multim2_CCA <- c(Results2$mcmc[[1]][,"multiParamCCA[2]"],Results2$mcmc[[2]][,"multiParamCCA[2]"])[IDs]
multim1_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[1]"],Results2$mcmc[[2]][,"multiParamRECCCA[1]"])[IDs]
multim2_REC <- c(Results2$mcmc[[1]][,"multiParamRECCCA[2]"],Results2$mcmc[[2]][,"multiParamRECCCA[2]"])[IDs]

# Function to compute mean and 95% credible interval
compute_summary <- function(data) {
  mean_val <- mean(data)
  lower_95 <- quantile(data, 0.025)
  upper_95 <- quantile(data, 0.975)
  
  return(sprintf("%.1f%% (%.1f%% - %.1f%%)", mean_val * 100, lower_95 * 100, upper_95 * 100))
}

# Compute per G_Score (element-wise means)
prevCCA <- sapply(1:length(IDs), function(x) modelprevCCA(100, x))
prevREC <- sapply(1:length(IDs), function(x) modelprevREC(100, x))
prevKK <- sapply(1:length(IDs), function(x) modelprevKK(100, x))

# Get the indices for the 1000 elements of prevKK closest to 0.1
indices <- order(abs(prevKK - 0.1))[1:1000]

# Subset the lists based on these indices
prevCCA <- prevCCA[,indices]
prevREC <- prevREC[,indices]
prevKK  <- prevKK[indices]

# Take the column means for prevCCA and prevREC
prevCCA_summary <- apply(prevCCA, 1, compute_summary)
prevREC_summary <- apply(prevREC, 1, compute_summary)

# Compute a single mean for prevKK
prevKK_summary <- compute_summary(prevKK)

# Define G_Score values (matching the unique elements)
G_Score <- c(2, 2.5, 3, 4)  # Assuming these are the G-Scores

# Create final summarised data frame
results <- data.frame(
  G_Score = G_Score,
  prevKK = prevKK_summary,  # Repeated single value
  prevREC = prevREC_summary,
  prevCCA = prevCCA_summary
)

# Print the results
print(results)

# Save to CSV
write.csv(results, "diagprev.csv", row.names = FALSE)
