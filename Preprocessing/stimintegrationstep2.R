library(Seurat)
library(SeuratData)
library(SeuratWrappers)

finalobj<-readRDS("~/scratch60/2021analysis/checkpoint/210129checkpoint2.rds")


finalobj <- NormalizeData(finalobj)
finalobj <- FindVariableFeatures(finalobj, selection.method = "vst", nfeatures = 2000)
features<-VariableFeatures(finalobj)
features<-grep("^TR.V", features, value = T, invert = T)

finalobj <- RunFastMNN(object.list = SplitObject(finalobj), features=features)
finalobj <- RunUMAP(finalobj, reduction = "mnn", dims = 1:30)
finalobj <- FindNeighbors(finalobj, reduction = "mnn", dims = 1:30)
finalobj <- FindClusters(finalobj)

saveRDS(finalobj,"~/scratch60/2021analysis/checkpoint/210130checkpoint3-integrationfastMNN.rds")





