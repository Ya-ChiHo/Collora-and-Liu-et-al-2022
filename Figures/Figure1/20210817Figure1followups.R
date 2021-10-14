library(Seurat)
library(dplyr)
library(ggplot2)

setwd("~/project/JointAnalyses/")

pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

#generating targetted dotplot
samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")

genelist<-c("GZMB","CCL5","TBX21","GZMK","GATA3","RORC","CTSH","CCR6","FOXP3","IL2RA","MKI67","CXCR5","CCR7","SELL","TCF7","percent.mt")

pubfig169p("Figures/Fig1stuff/clusterumap")
DimPlot(samples, group.by = "cluster_names", label=T)
dev.off()

pubfig169p("Figures/Fig1stuff/clusterumap")
DimPlot(samples, group.by = "memory", label=T)
dev.off()



samples$cluster_names<-factor(samples$cluster_names, levels = rev(levels(samples$cluster_names)))
Idents(samples)<-samples$cluster_names

pubfig169p("Figures/Fig1stuff/Dotplotinfectionsplit")
DotPlot(samples, features = genelist, cols = c("#2297E6","#DF536B","grey62"), split.by = "infection")
dev.off()


#checking if there are more "th1" in proliferating during viremia 
samples$TBX21<-samples@assays$RNA@data["TBX21",]>0
samples2<-subset(samples, idents = "CD4-Proliferating")
meta<-samples2@meta.data
meta<-group_by(meta, infection,UID,TBX21)%>%summarise(n=n())%>%mutate(freq=n/sum(n))
meta<-meta[meta$TBX21==TRUE,]
pubfigsqp("Figures/Fig1stuff/TBX21positive")
ggplot(meta, aes(x=infection, y=freq ))+geom_jitter()+theme_classic()
dev.off()

samples$log2fc<-log2(samples$results_freq*10000)

pubfig169p("Figures/Fig1stuff/clonesizeumap")
FeaturePlot(samples, features="log2fc", min.cutoff='q5', max.cutoff='q95')
dev.off()

#things are not normal so we cant T test them
wilcox.test(meta$freq[meta$infection=="Viremic"], meta$freq[meta$infection=="Suppressed"])
wilcox.test(meta$freq[meta$infection=="Viremic"],meta$freq[meta$infection=="Uninfected"])
wilcox.test(meta$freq[meta$infection=="Suppressed"],meta$freq[meta$infection=="Uninfected"])


