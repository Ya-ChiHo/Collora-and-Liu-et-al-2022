library(Seurat)
library(qlcMatrix)
setwd("~/project/JointAnalyses/")
#this script just generates gene-gene correlation matrices for unstimulated

#infection split

  samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
  samples<-SplitObject(samples, split.by = "infection")
  
  for(i in 1:length(samples)){
    
      data<-samples[[i]]@assays$RNA@data
      cor<-corSparse(t(data))
      rownames(cor)<-colnames(cor)<-rownames(data)
      saveRDS(cor,paste("Objects/cor/",names(samples)[[i]],"_", names(samples)[[i]],"_","unstimulated.rds", sep=""))
    }
#infection2 split
  samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
  samples<-SplitObject(samples, split.by = "infection2")
  
  for(i in 1:length(samples)){
    
    data<-samples[[i]]@assays$RNA@data
    cor<-corSparse(t(data))
    rownames(cor)<-colnames(cor)<-rownames(data)
    saveRDS(cor,paste("Objects/cor/",names(samples)[[i]],"_", names(samples)[[i]],"_","unstimulated.rds", sep=""))
  }
