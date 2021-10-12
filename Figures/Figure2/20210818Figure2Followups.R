library(Seurat)
library(DescTools)
library(dplyr)
library(ggplot2)

setwd("~/project/JointAnalyses/")

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")

pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

meta<-samples@meta.data

#based on nonbulk 

meta<-samples@meta.data
meta<-meta[!is.na(meta$results_count),]
#counted clones in each cluster to avoid double counting for cells that go across cluisters
meta<-group_by(meta, UID,infection, cluster_names, junction)%>%summarise(n=n())%>%summarise(gini=Gini(n))
#p values relative to naive
for(i in names(table(meta$cluster_names))){
  print(i)
  print(wilcox.test(meta$gini[meta$cluster_names==i],meta$gini[meta$cluster_names=="CD4-Naive"] )$p.value)
}

meta$infection<-factor(meta$infection, levels = c("Viremic","Suppressed","Uninfected"))
pubfig169p("Figures/Clonalfigs/ginisingleunstimulated")
ggplot(meta, aes(x=cluster_names, y=gini, color=infection))+geom_point(position=position_dodge(0.25))+theme_classic()
dev.off()
