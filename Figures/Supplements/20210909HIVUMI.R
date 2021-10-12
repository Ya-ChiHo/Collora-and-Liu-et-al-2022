library(Seurat)
library(ggplot2)

setwd("~/project/JointAnalyses/")

unstim<-readRDS("Objects/20210519UnstimulatedMaster.rds")

stim<-readRDS("Objects/210623masterstimulated.rds")

unstim<-subset(unstim, HIV_infected=="infected")
stim<-subset(stim, HIVpositive==TRUE)
#normalizing HIVUMI to library depth and taking the log 
unstim$normUMIHIV<-log10((unstim$HIVUMI/unstim$nCount_RNA)*10000)
stim$normUMIHIV<-log10((stim$HIVUMI/stim$nCount_RNA)*10000)
stim$stageag<-paste(stim$stage, stim$antigen2, sep="_")
stim$stageag<-factor(stim$stageag, levels = c("Viremic_pCMVresp","Suppressed_pCMVresp","Viremic_pHIVresp","Suppressed_pHIVresp","Viremic_pmemory","Suppressed_pmemory"))
p1<-VlnPlot(unstim, "normUMIHIV", group.by = "infection", y.max = 3.5)
p2<-VlnPlot(stim, "normUMIHIV", group.by = "stageag",y.max=3.5)
#plotting
pdf("Supp/HIVUMI.pdf")
p1
p2
dev.off()
