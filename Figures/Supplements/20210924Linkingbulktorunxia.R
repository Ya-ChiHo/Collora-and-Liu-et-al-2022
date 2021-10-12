library(dplyr)
library(Seurat)
library(ggplot2)

setwd("~/project/JointAnalyses/")

samples<-readRDS("Objects/210623masterstimulated.rds")
samples$UID<-case_when( samples$UID == "HD2004_2" ~"HD2022_2", samples$UID=="HD2022_2"~"HD2004_2", TRUE~as.character(samples$UID))  
RunxiaMeta<-samples@meta.data
BulkMeta<-readRDS("Objects/20210601BulkQuant.rds")




# For bulk data, always use time point + individual
# First based on time point and individual. IE I compared your suppressed time point for patient 799 to 
# my suppressed timepoint for patient 799. This is more accurate in the cluster assignment (since it's going to be 
# the same phenotypically) but will miss things with smaller clone size (since they're less likely to be picked up). 


#first by UID
colnames(BulkMeta)<-paste(colnames(BulkMeta),"test", sep="_")
BulkUIDSplit<-split(BulkMeta, BulkMeta$ID2_test)
RunxiaUIDSplit<-split(RunxiaMeta, RunxiaMeta$UID)

todo<-names(RunxiaUIDSplit)
#check that all are present in both; they are not matching
# bulk doesn't have M2_2; but it has M1_1 and M1_2 which are not in RunxiaMeta
todo==names(RunxiaUIDSplit)

for (i in todo){
  #first generate a lookup table for each sample
  lookuptable<-BulkUIDSplit[[i]]
  colnames(lookuptable)<-paste("bulk",colnames(lookuptable), sep="_")
  RunxiaUIDSplit[[i]]<-merge(RunxiaUIDSplit[[i]],lookuptable, by.y = "bulk_JUNCTION_test",by.x = "TRB_junction",all.x=TRUE)
}

# bind function when the number of columns are not equal
RunxiaUIDSplit<-bind_rows(RunxiaUIDSplit)

# add rownames
rownames(RunxiaUIDSplit)<-RunxiaUIDSplit$cell
meta<-group_by(RunxiaUIDSplit, UID, antigen2)%>%summarise(cells=n(),mean_genes=mean(nFeature_RNA), mean_umi=mean(nCount_RNA),TCRsingle=sum(!is.na(TRB_junction)), clonal=sum(bulk_n_test
                                                                                                                                                                          >1& !is.na(bulk_n_test)), HIVRNApos=sum(HIVpositive))

saveRDS(RunxiaUIDSplit,"Objects/20210924Runxiameta.rds")

x<-readRDS("Objects/20210924Runxiameta.rds")
x<-group_by(x, TRB_junction, PT, antigen2)%>%slice_max(n=1, order_by=bulk_n_test, with_ties=FALSE)
x<-x[x$bulk_n_test>1& !is.na(x$bulk_n_test
                             ),]
table(x$antigen2)
