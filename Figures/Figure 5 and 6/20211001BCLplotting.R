library(Seurat)
library(dplyr)
setwd("~/project/JointAnalyses/")


pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

samples<-readRDS("Objects/210623masterstimulated.rds")
samples$cluster_annotation<-factor(samples$cluster_annotation, levels =  c("CD4-GZMB Th1","CD4-IL2 Th1","CD4-Lymphotoxin",
                                                                           "CD4-CD45RA","CD4-HLA-DR","CD4-GZMH","CD4-Treg",
                                                                           "CD4-IFN","CD4-Memory-1","CD4-Memory-2","CD4-Memory-3",
                                                                           "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-HSP"))
Idents(samples)<-samples$cluster_annotation
DefaultAssay(samples)<-"RNA"
samples<-ScaleData(samples, c("BCL2","BCL2L1","BCL2A1","BCL2L11","BCL3","BCL7C","BCLAF1","BCL7B","BCL10","BCL9L","BCL11B"))
pubfigsqp("Figures/Fig4stuff/stimHIVBCLdotplot")
DotPlot(samples, cols = c("#B8B8B8","#DF536B"), features = c("BCL2","BCL2L1","BCL2A1","BCL2L11","BCL3","BCL7C","BCLAF1","BCL7B","BCL10","BCL9L","BCL11B"),assay="RNA", split.by = "HIVpositive")
dev.off()
samples$keep<-samples$cluster_annotation %in% c("CD4-Lymphotoxin", "CD4-GZMH","CD4-GZMB Th1","CD4-IL2 Th1","CD4-CD45RA","CD4-HLA-DR")
samples$antigen5<-case_when(grepl("non",samples$antigen) ~"Memory", T~as.character(samples$antigen))

samples2<-subset(samples, keep==TRUE)
pubfigsqp("Figures/Fig4stuff/stimclusterBCLdotplot")
DotPlot(samples2, features = c("BCL2","BCL2L1","BCL2A1","BCL2L11","BCL3","BCL7C","BCLAF1","BCL7B","BCL10","BCL9L","BCL11B"), split.by = "antigen5", cols = c("#6CBAEF","#B8B8B8","#DF536B"))
dev.off()

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
samples<-ScaleData(samples, features=c("BCL2","BCL2L1","BCL2A1","BCL2L11","BCL3","BCL7C","BCLAF1","BCL7B","BCL10","BCL9L","BCL11B"))
pubfigsqp("Figures/Fig4stuff/unstimHIVBCLdotplot")
DotPlot(samples, cols = c("#B8B8B8","#DF536B"), features = c("BCL2","BCL2L1","BCL2A1","BCL2L11","BCL3","BCL7C","BCLAF1","BCL7B","BCL10","BCL9L","BCL11B"),assay="RNA", split.by = "HIV_infected")
dev.off()
pubfigsqp("Figures/Fig4stuff/unstimclusterBCLdotplot")
DotPlot(samples, features = c("BCL2","BCL2L1","BCL2A1","BCL2L11","BCL3","BCL7C","BCLAF1","BCL7B","BCL10","BCL9L","BCL11B"), split.by = "infection", cols = c("#6CBAEF","#DF536B","#B8B8B8"))
dev.off()
