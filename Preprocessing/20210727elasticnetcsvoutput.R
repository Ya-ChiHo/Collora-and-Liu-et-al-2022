library(Seurat)

#this is generating csv files of gene expression matrix and meta.data, tried it with Seurat disk and it didnt work due to an unspecified error

setwd("~/project/JointAnalyses/")

#generate for clone size unstimulated, by infection condition

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
meta<-readRDS("Objects/metadata with both clone timepoint quantifications")
samples<-DietSeurat(samples, assays = "RNA")

samples$keep<-!is.na(samples$results_count)
samples<-subset(samples, keep==TRUE)
samples<-SplitObject(samples, split.by = "infection")

for(i in 1:length(samples)){
  write.csv(samples[[i]]@assays$RNA@data, paste("sklearn/clonality/unstimulated_data_",names(samples)[[i]], sep = ""), quote=FALSE)
  write.table(samples[[i]]@meta.data, paste("sklearn/clonality/unstimulated_meta_",names(samples)[[i]], sep = ""),sep="\t", quote=FALSE) 
}
