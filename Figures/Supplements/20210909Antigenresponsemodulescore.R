library(Seurat)
library(ggplot2)
library(msigdbr)
#reading in files and the modules from msigdb
setwd("~/project/JointAnalyses/")

unstim<-readRDS("Objects/20210519UnstimulatedMaster.rds")

stim<-readRDS("Objects/210623masterstimulated.rds")

setwd("~/project/JointAnalyses/")

m_df_H<- msigdbr(species = "Homo sapiens", category = "C2")
fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
#scoring modules, calculating means based on cluster 
unstim<-AddModuleScore(unstim, list(fgsea_sets$GOLDRATH_ANTIGEN_RESPONSE))
stim<-AddModuleScore(stim, list(fgsea_sets$GOLDRATH_ANTIGEN_RESPONSE))

FeaturePlot(stim, "Cluster1", min.cutoff = "q5",max.cutoff = 'q95')
FeaturePlot(unstim, "Cluster1", min.cutoff = "q5",max.cutoff = 'q95')


metau<-unstim@meta.data

metau2<-group_by(metau, infection, cluster_names)%>%summarise(antigenmod=mean(Cluster1))
metau2$cluster_names<-factor(metau2$cluster_names, levels = rev(levels(metau2$cluster_names)))
#doing plotting for unstimed and then the two stimed plots
pubfig169p("Figures/Fig6stuff/activationmodulescoreunstim")
ggplot(metau2, aes(x=infection, y=cluster_names, fill=antigenmod))+geom_tile()+theme_classic()
dev.off()


stim<-AddModuleScore(stim, list(fgsea_sets$GOLDRATH_ANTIGEN_RESPONSE))
stim$cluster_annotation<-factor(stim$cluster_annotation, levels = rev(c("CD4-GZMB Th1","CD4-IL2 Th1","CD4-Lymphotoxin",
  "CD4-CD45RA","CD4-HLA-DR","CD4-GZMH","CD4-Treg",
  "CD4-IFN","CD4-Memory-1","CD4-Memory-2","CD4-Memory-3",
  "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-HSP")))
Idents(stim)<-stim$cluster_annotation
stim$stage<-case_when(grepl("HD", stim$stage)~"Uninfected", T~as.character(stim$stage))
#memory clusters
stimmem<-subset(stim, idents = c("CD4-GZMH","CD4-HSP","CD4-Memory-1","CD4-IFN","CD4-Memory-2","CD4-Memory-3",
                                 "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-Treg"))
#stimulated clusters
stimag<-subset(stim, idents = c("CD4-CD45RA","CD4-GZMB Th1","CD4-HLA-DR","CD4-IL2 Th1","CD4-Lymphotoxin"))

stimmem<-AddModuleScore(stimmem, list(fgsea_sets$GOLDRATH_ANTIGEN_RESPONSE))
stimag<-AddModuleScore(stimag, list(fgsea_sets$GOLDRATH_ANTIGEN_RESPONSE))

FeaturePlot(stimmem, "Cluster1", min.cutoff='q5', max.cutoff='q95')
FeaturePlot(stimag, "Cluster1", min.cutoff='q5', max.cutoff='q95')
pubfig169p("Figures/Fig6stuff/activationmodulescorestimmem")
group_by(stimmem@meta.data, stage, cluster_annotation)%>%summarise(antigenmod=mean(Cluster1))%>%ggplot(aes(x=stage, y=cluster_annotation, fill=antigenmod))+geom_tile()+theme_classic()
dev.off()
pubfig169p("Figures/Fig6stuff/activationmodulescorestimag")
group_by(stimag@meta.data, stage, cluster_annotation)%>%summarise(antigenmod=mean(Cluster1))%>%ggplot(aes(x=stage, y=cluster_annotation, fill=antigenmod))+geom_tile()+theme_classic()
dev.off()
#plots for memory status 
pubfig169p("Figures/Fig1stuff/memory.pdf")
DimPlot(unstim, group.by = "memory", raster = FALSE)
dev.off()

pubfig169p("Figures/Fig3stuff/memory.pdf")
DimPlot(stim, group.by = "memory", raster = FALSE)
dev.off()
