library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(BiocParallel)
options(MulticoreParam=quote(MulticoreParam(workers=10)))

#startng by reading in my files 
finalobj<-readRDS("~/project/RNA/checkpoint/samplesforDSA.rds")
meta<-readRDS("~/project/RNA/checkpoint/metaforDSA.rds")
#subset to my post QC cells
finalobj$keep<-colnames(finalobj)%in%rownames(meta)
finalobj<-subset(finalobj, keep==TRUE)
meta$UID[meta$UID=="M1_0"]<-"M2_0"
finalobj<-AddMetaData(finalobj, meta)
DefaultAssay(finalobj)<-"RNA"
#splite up to run sctransform and call variable genes, save the intermediate
finalobj<-SplitObject(finalobj, "UID")
finalobj<-lapply(finalobj, SCTransform)
features<-SelectIntegrationFeatures(finalobj)
features<-grep("^TR.V", features, value = T, invert = T)
saveRDS(finalobj,"~/project/RNA/checkpoint/samplesforMNN.rds")

#runfast MNN to integrate and save the result
finalobj <- RunFastMNN(object.list = finalobj, features=features)
finalobj <- RunUMAP(finalobj, reduction = "mnn", dims = 1:30)
finalobj <- FindNeighbors(finalobj, reduction = "mnn", dims = 1:30)
finalobj <- FindClusters(finalobj)

saveRDS(finalobj,"~/project/RNA/checkpoint/FASTMNN.rds")

