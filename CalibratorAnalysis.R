############ Calibrator Data #####
### J. M. Prada
rm(list=objects())  #Erase workspace
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

dt <- read.csv("Calibrators_Mayuge_recPOC_CCA.csv")

dtf <- data.frame("Batch" = factor(rep(dt$Batch_recPOC_CCA,3)),
                  "Box" = factor(rep(dt$Box_recPOC_CCA,3)),
                  "Calibrator" = factor(c(rep(1,nrow(dt)),rep(2,nrow(dt)),rep(3,nrow(dt)))),
                  "Value" = c(dt$Gsc_Cal_1,dt$Gsc_Cal_2,dt$Gsc_Cal_3) )

## Linear regression
res <- lm(Value ~ Batch + Calibrator, data=dtf)
summary(res)

## Linear regression adding box as a random effect
library(lme4)
res <- lmer(Value ~ Batch + Calibrator + (1|Box), data=dtf)
summary(res)

## Plot batch x Calibrator
pdf("BatchXCal.pdf",width=7,height=4.64)
par(font=2, cex.axis=1, lwd=1.5, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

boxplot(Value ~ Batch*Calibrator , data=dtf, axes=F, xlab="",ylab="",
        outline = F,ylim = c(2,9)) 

abline(v=c(3.5,6.5), lty = "dashed")

for(i in 1:3){ #batch
  for(j in 1:3){ #cal
    index = i + (j-1)*3
    vals = subset(dtf,Batch==i & Calibrator==j)
    jit <- jitter(rep(index,nrow(vals)), amount = .25)
    points(jit,vals$Value, pch =19)
  }
}

axis(1, at = 1:9,labels = paste("Bt",rep(1:3,3)))
mtext(c("Calibrator 1","Calibrator 2","Calibrator 3"), 1, 1.2, at=c(2,5,8), las=1)

axis(2)
mtext("G-Score",side=2,cex=1,line=1.2)

dev.off()