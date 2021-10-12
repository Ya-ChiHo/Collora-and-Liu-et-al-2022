library(Seurat)
library(SeuratDisk)
library(dplyr)
#goal of script is to make anndata objects for scRFE
setwd("project/JointAnalyses/")
samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
JackClones<-read.csv("sklearn/ClonesHIVJackData.csv", col.names = "cells")
RunxiaClones<-read.csv("sklearn/ClonesHIVRunxiaData.csv", col.names = "cells")
samples$HIVclone<-case_when(colnames(samples)%in%JackClones$cells ~"resting", 
                            colnames(samples)%in%RunxiaClones$cells ~"stimulated",
                            TRUE ~ "nonHIVclone")
#15k cells remain 
samples<-subset(samples, results_count>1)
samples<-subset(samples, infection!="Uninfected")

#subsets are a little different so we'll balance it out
samples2<-SplitObject(samples, split.by = "HIVclone")
props_current<-table(samples2$nonHIVclone$cluster_names)/ncol(samples2$nonHIVclone)
table(samples2$stimulated$cluster_names)/ncol(samples2$stimulated)
props_goal<-table(samples2$resting$cluster_names)/ncol(samples2$resting)
#4 fold enrichment for GZMB Th1 is out most imbalanced so that's going to be the limit
total<-sort(props_goal/props_current, decreasing = TRUE)[1]
total<-as.numeric(1/props_goal[names(total)])*table(samples2$nonHIVclone$cluster_names)[names(total)]
props_goal<-floor(props_goal*total)
props_goal<-props_goal[props_goal>0]


samples2<-subset(samples, HIVclone=="nonHIVclone")
samples2<-SplitObject(samples2, "cluster_names")

for (i in names(props_goal)){
  cells_keep<-sample(colnames(samples2[[i]]), size = props_goal[[i]],replace = FALSE)
  samples2[[i]]<-subset(samples2[[i]], cells=cells_keep)
}
samples3<-subset(samples, HIVclone!="nonHIVclone")
#bring objects back together
for(i in names(props_goal)){
  samples3<-merge(x=samples3, y=samples2[[i]])
}

#writing one copy with stimulated clones, this will be the object used for testing the resulting gene model
SaveH5Seurat(samples3, filename = "sklearn/20210521allclonesforscRFE.h5Seurat")
Convert("sklearn/20210521allclonesforscRFE.h5Seurat", dest = "h5ad")
#writing a second copy without stimulated clones, this will be the training object
samples3<-subset(samples3, HIVclone!="stimulated")
SaveH5Seurat(samples3, filename = "sklearn/20210521restingclonesforscRFE.h5Seurat")
Convert("sklearn/20210521restingclonesforscRFE.h5Seurat", dest = "h5ad")

