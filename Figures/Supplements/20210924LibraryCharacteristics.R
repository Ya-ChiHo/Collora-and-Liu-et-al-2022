library(Seurat)
library(dplyr)
library(tidyr)

setwd("~/project/JointAnalyses/")

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
#fixing a labeling issue
samples$PT<-case_when( samples$PT == "M2" |samples$PT=="M1" ~"HD2023", samples$PT=="M3"~"HD2004", TRUE~as.character(samples$PT))
samples$UID<-paste(samples$PT, samples$timepoint)

meta<-samples@meta.data
#grouping by library and writing to file
meta<-group_by(meta, UID)%>%summarise(cells=n(),mean_genes=mean(nFeature_RNA), mean_umi=mean(nCount_RNA),TCRsingle=sum(!is.na(junction)), clonal=sum(results_count>1&!is.na(junction)), HIVRNApos=sum(HIV_infected=="infected"))

write.table(meta, "Supp/unstimulated_library_chars.tsv", sep="\t", row.names = FALSE, quote=FALSE)

#updated stimulated meta with right bulkn
samples<-readRDS("Objects/20210924Runxiameta.rds")
#simplify annotation
samples$antigen2<-case_when(samples$antigen2=="hdCMVresp"|samples$antigen2=="pCMVresp"~"CMV", samples$antigen2=="hdmemory" | samples$antigen2=="pmemory"~"Memory",samples$antigen2=="pHIVresp"~"HIV", TRUE~"resting")
#grouping by library and writing to file

meta<-group_by(samples, UID, antigen2)%>%summarise(cells=n(),mean_genes=mean(nFeature_RNA), mean_umi=mean(nCount_RNA),TCRsingle=sum(!is.na(TRB_junction)), clonal=sum(bulk_n_test>1& !is.na(bulk_n_test)), HIVRNApos=sum(HIVpositive))

write.table(meta, "Supp/stimulated_library_chars.tsv", sep="\t", row.names = FALSE, quote=FALSE)


