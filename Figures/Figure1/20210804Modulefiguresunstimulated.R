library(Seurat)
library(dplyr)
library(ggplot2)

setwd("~/project/JointAnalyses/")

pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}


#unstim scoring on unstim
mods<-readRDS("Objects/20210803UnstimulatedConsensusmods.rds")

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
samples$infection<-factor(samples$infection, c("Viremic", "Suppressed", "Uninfected"))

samples<-AddModuleScore(samples, mods, name="unstimulatedmods")


####mod4, TNF mod, related to second half of figure 1#####
FeaturePlot(samples, "unstimulatedmods4",min.cutoff='q5', max.cutoff='q95')

pubfiglongp("Figures/Modulefigs/UMAPmod4")
FeaturePlot(samples, "unstimulatedmods4",min.cutoff=0.2, max.cutoff=1, split.by = "infection")
dev.off()

VlnPlot(samples, split.by = "infection", features = "unstimulatedmods4", pt.size = 0)

pubfigsqp("Figures/Modulefigs/vlnplot4GZMK")
VlnPlot(samples, split.by = "infection", features = "unstimulatedmods4", pt.size = 0, idents = "CD4-GZMK Th1")
dev.off()
#sig testing in GZMK
{#p<2.2e-16
  wilcox.test(samples$unstimulatedmods4[samples$infection=="Viremic"& samples$cluster_names=="CD4-GZMK Th1"], 
              samples$unstimulatedmods4[samples$infection=="Suppressed"& samples$cluster_names=="CD4-GZMK Th1"])
  #p<2.2e-16
  
  wilcox.test(samples$unstimulatedmods4[samples$infection=="Viremic"& samples$cluster_names=="CD4-GZMK Th1"], 
              samples$unstimulatedmods4[samples$infection=="Uninfected"& samples$cluster_names=="CD4-GZMK Th1"])
  #p<2.2e-16
  wilcox.test(samples$unstimulatedmods4[samples$infection=="Uninfected"& samples$cluster_names=="CD4-GZMK Th1"], 
              samples$unstimulatedmods4[samples$infection=="Suppressed"& samples$cluster_names=="CD4-GZMK Th1"])
  
}


pubfigsqp("Figures/Modulefigs/vlnplot4Proliferating")
VlnPlot(samples, split.by = "infection", features = "unstimulatedmods4", pt.size = 0, idents = "CD4-CXCR5 Memory")
dev.off()

#sig testing in Proliferating
{
  #p=0.034
  wilcox.test(samples$unstimulatedmods4[samples$infection=="Viremic"& samples$cluster_names=="CD4-CXCR5 Memory"], 
              samples$unstimulatedmods4[samples$infection=="Suppressed"& samples$cluster_names=="CD4-CXCR5 Memory"])
  #p<2.2e-16
  wilcox.test(samples$unstimulatedmods4[samples$infection=="Viremic"& samples$cluster_names=="CD4-CXCR5 Memory"], 
              samples$unstimulatedmods4[samples$infection=="Uninfected"& samples$cluster_names=="CD4-CXCR5 Memory"])
  #p<2.2e-16
  wilcox.test(samples$unstimulatedmods4[samples$infection=="Uninfected"& samples$cluster_names=="CD4-CXCR5 Memory"], 
              samples$unstimulatedmods4[samples$infection=="Suppressed"& samples$cluster_names=="CD4-CXCR5 Memory"])
}




