#generating cor matricies for gsea and matricies for elastic net regression 
library(Seurat)
library(qlcMatrix)
setwd("~/project/JointAnalyses/")
#loading objects, spliting them by participant 
samples<-readRDS("Objects/210623masterstimulated.rds")
samples$keep<-!is.na(samples$bulk_n)
samples<-subset(samples, keep==TRUE)

samples<-SplitObject(samples, split.by = "antigen2")
samples<-lapply(samples, SplitObject, split.by="UID")
#this sort of loop occurs three times, in each, first the clone size is added as a feature, then a correlation matrix is generated and saved

for(i in 1:length(samples)){
  for(j in 1:length(samples[[i]])){
    data<-samples[[i]][[j]]@assays$RNA@data
    data<-rbind(data, samples[[i]][[j]]$bulk_n)
    rownames(data)<-c(rownames(data)[1:nrow(data)-1], "clonesize")
    cor<-corSparse(t(data))
    rownames(cor)<-colnames(cor)<-rownames(data)
    saveRDS(cor,paste("Objects/cor/",names(samples)[[i]],"_", names(samples[[i]])[[j]],"_","stimulated.rds", sep=""))
  }
  }

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
samples$keep<-!is.na(samples$results_count)
samples<-subset(samples, keep==TRUE) 
samples<-SplitObject(samples, split.by = "UID")


for(i in 1:length(samples)){
    data<-samples[[i]]@assays$RNA@data
    data<-rbind(data, samples[[i]]$results_count)
    rownames(data)<-c(rownames(data)[1:nrow(data)-1], "clonesize")
    cor<-corSparse(t(data))
    rownames(cor)<-colnames(cor)<-rownames(data)
    saveRDS(cor,paste("Objects/cor/",names(samples)[[i]],"_", "unstimulated.rds", sep=""))
    
}

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
samples$keep<-!is.na(samples$results_count)
samples<-subset(samples, keep==TRUE)
samples<-SplitObject(samples, split.by = "infection")


for(i in 1:length(samples)){
  data<-samples[[i]]@assays$RNA@data
  data<-rbind(data, samples[[i]]$results_count)
  rownames(data)<-c(rownames(data)[1:nrow(data)-1], "clonesize")
  cor<-corSparse(t(data))
  rownames(cor)<-colnames(cor)<-rownames(data)
  saveRDS(cor,paste("Objects/cor/",names(samples)[[i]],"_", "unstimulated.rds", sep=""))
  
}


