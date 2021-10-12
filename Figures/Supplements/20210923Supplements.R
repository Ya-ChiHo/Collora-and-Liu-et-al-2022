library(Seurat)
library(ggplot2)
library(dplyr)
setwd("~/project/JointAnalyses/")

pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

####heatmaps####

#unstim
samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
DefaultAssay(samples)<-"RNA"
samples<-NormalizeData(samples)
markers<-FindAllMarkers(samples)

dir.create("Supp")

write.table(markers, "Supp/Unstimulatedmarkers.tsv", sep="\t", quote=FALSE)
topmarks<-group_by(markers, cluster)%>%slice_max(order_by=avg_log2FC, n=12)
topmarks<-unique(topmarks$gene)
samples<-ScaleData(samples, features = topmarks)

pubfig169p("Supp/Unstimulatedheatmap")
DoHeatmap(samples, topmarks, angle=90)
dev.off()

#Stim
samples<-readRDS("Objects/210623masterstimulated.rds")
DefaultAssay(samples)<-"RNA"
samples<-DietSeurat(samples, assays="RNA")
samples<-NormalizeData(samples)

samples$cluster_annotation<-factor(samples$cluster_annotation, levels =  c("CD4-GZMB Th1","CD4-IL2 Th1","CD4-Lymphotoxin",
                                                                           "CD4-CD45RA","CD4-HLA-DR","CD4-GZMH","CD4-Treg",
                                                                           "CD4-IFN","CD4-Memory-1","CD4-Memory-2","CD4-Memory-3",
                                                                           "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-HSP"))

Idents(samples)<-samples$cluster_annotation
markers<-FindAllMarkers(samples)

write.table(markers, "Supp/Stimulatedmarkers.tsv", sep="\t", quote=FALSE)
topmarks<-group_by(markers, cluster)%>%slice_max(order_by=avg_log2FC, n=12)
topmarks<-unique(topmarks$gene)
samples<-ScaleData(samples, features = topmarks)


pubfig169p("Supp/Stimulatedheatmap")
DoHeatmap(samples, topmarks, angle=90)
dev.off()


####integration UMAPS####

#unstim
samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
samples$PT<-case_when( samples$PT == "M2" |samples$PT=="M1" ~"HD2023", samples$PT=="M3"~"HD2004", TRUE~as.character(samples$PT))


pubfig169p("Supp/PTcolor_unstim")
DimPlot(samples, group.by = "PT", raster=FALSE)
dev.off()

pubfig169p("Supp/PTsplit_unstim")
DimPlot(samples, group.by="PT",split.by = "PT",ncol=3, raster=FALSE)+
  scale_color_manual(values=c(`236` = "#0E5381",`640`= "#188EDC",`829`= "#91CCF3",
                              `739` = "#29711E",`799` = "#46C133",`910` = "#A8E59E",
                              'HD2004' = "#666666",`HD2023` = "#CCCCCC"))
dev.off()

#stim
samples<-readRDS("Objects/210623masterstimulated.rds")
pubfig169p("Supp/PTcolor_stim")
DimPlot(samples, group.by = "PT", raster=FALSE)
dev.off()

pubfig169p("Supp/PTsplit_stim")
DimPlot(samples, group.by="PT", split.by="PT", ncol=3, raster=FALSE)+
  scale_color_manual(values=c(`236` = "#0E5381",`640`= "#188EDC",`829`= "#91CCF3",
                              `739` = "#29711E",`799` = "#46C133",`910` = "#A8E59E",
                              'HD2004' = "#666666",`HD2022` = "#CCCCCC"))
dev.off()

####HIV Markers####

#both
markers<-readRDS("Objects/20210602markersofHIVinfectedcells.rds")

for(i in names(markers)){
  write.table(markers[[i]], paste("Supp/",i, "_HIVRNAMarkers.tsv", sep=""),sep="\t", quote=FALSE)
}



