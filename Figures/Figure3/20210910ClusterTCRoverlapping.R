library(Seurat)
library(dplyr)
library(ggplot2)
pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}
setwd("~/project/JointAnalyses/")

#we want to strengthen the connection between the two datasets more generally, but specifically with regard to the GZMB and GZMB/H clusters 

#the first goal of this analysis is to identify 1) the cluster pairs based on TCR overlap in the antigen specific, and memory clusters for stimulated and unstimulated datasets 

#the second goal of this analysis is to identify the antigen specificity pair between stim and unstim datasets 

####data loading, cleaning, formating####

unstim<-readRDS("Objects/20210519UnstimulatedMaster.rds")
stim<-readRDS("Objects/210623masterstimulated.rds")

#subset both to not have uninfected 

unstim<-subset(unstim, infection!="Uninfected")
stim<-subset(stim, stage=="Suppressed"|stage=="Viremia")

####goal 1 - cluster overlaps####
#first things first we're going to subset the metadata to only include the TCRs present in both datasets
stimmeta<-stim@meta.data[stim$TRB_junction%in%unstim$junction,]
unstimmeta<-unstim@meta.data[unstim$junction%in%stim$TRB_junction,]
#next remove those with NA
stimmeta<-stimmeta[!is.na(stimmeta$TRB_junction),]
unstimmeta<-unstimmeta[!is.na(unstimmeta$junction),]

#all as one
#split only by cluster
stimmeta<-split(stimmeta, stimmeta$cluster_annotation)
unstimmeta<-split(unstimmeta, unstimmeta$cluster_names)

#what cluster in unstim meta has the most overlap 
results<-list()
for(i in names(unstimmeta)){
  cur<-c()
  for(j in names(stimmeta)){
    cur<-c(cur, sum(unstimmeta[[i]]$junction%in%stimmeta[[j]]$TRB_junction)/nrow(unstimmeta[[i]]))
  }
  results[[i]]<-cur
}
names(results)<-names(unstimmeta)
results<-unlist(results, recursive = FALSE)

results<-as.data.frame(matrix(data = results, ncol  = 10))
rownames(results)<-names(stimmeta)
colnames(results)<-names(unstimmeta)
pubfigsqp("Figures/Fig3stuff/overlapunstimtostim")
heatmap(as.matrix(results),scale="col")
dev.off()

#first things first we're going to subset the metadata to only include the TCRs present in both datasets
stimmeta<-stim@meta.data[stim$TRB_junction%in%unstim$junction,]
unstimmeta<-unstim@meta.data[unstim$junction%in%stim$TRB_junction,]
#next remove those with NA
stimmeta<-stimmeta[!is.na(stimmeta$TRB_junction),]
unstimmeta<-unstimmeta[!is.na(unstimmeta$junction),]

#all as one
#split only by cluster
stimmeta<-split(stimmeta, stimmeta$cluster_annotation)
unstimmeta<-split(unstimmeta, unstimmeta$cluster_names)

#what cluster in unstim meta has the most overlap 
results<-list()
for(i in names(stimmeta)){
  cur<-c()
  for(j in names(unstimmeta)){
    cur<-c(cur, sum(stimmeta[[i]]$TRB_junction%in%unstimmeta[[j]]$junction)/nrow(stimmeta[[i]]))
  }
  results[[i]]<-cur
}
names(results)<-names(stimmeta)
results<-unlist(results, recursive = FALSE)
results<-as.data.frame(matrix(data = results, ncol  = 15))
results<-t(results)
rownames(results)<-names(stimmeta)
colnames(results)<-names(unstimmeta)
results<-t(results)
pubfigsqp("Figures/Fig3stuff/overlapstimtounstim")
heatmap(as.matrix(results),scale="col")
dev.off()
#with each individual 

#first things first we're going to subset the metadata to only include the TCRs present in both datasets
stimmeta<-stim@meta.data[stim$TRB_junction%in%unstim$junction,]
unstimmeta<-unstim@meta.data[unstim$junction%in%stim$TRB_junction,]
#next remove those with NA
stimmeta<-stimmeta[!is.na(stimmeta$TRB_junction),]
unstimmeta<-unstimmeta[!is.na(unstimmeta$junction),]

#all as one
#split by a combo of cluster and individual
stimmeta<-split(stimmeta, paste(stimmeta$cluster_annotation, stimmeta$PT))
unstimmeta<-split(unstimmeta, paste(unstimmeta$cluster_names, unstimmeta$PT))

#what cluster in unstim meta has the most overlap 
results<-list()
for(i in names(unstimmeta)){
  cur<-c()
  for(j in names(stimmeta)){
    cur<-c(cur, sum(unstimmeta[[i]]$junction%in%stimmeta[[j]]$TRB_junction)/nrow(unstimmeta[[i]]))
  }
  results[[i]]<-cur
}
names(results)<-names(unstimmeta)
results<-unlist(results, recursive = FALSE)

results<-as.data.frame(matrix(data = results, ncol  = 60))
rownames(results)<-names(stimmeta)
colnames(results)<-names(unstimmeta)
pubfigsqp("Figures/Fig3stuff/overlapunstimtostimPT")
heatmap(as.matrix(results), scale = "col")
dev.off()

####goal 2 - antigen overlaps####

#first things first we're going to subset the metadata to only include the TCRs present in both datasets
stimmeta<-stim@meta.data[stim$TRB_junction%in%unstim$junction,]
unstimmeta<-unstim@meta.data[unstim$junction%in%stim$TRB_junction,]
#next remove those with NA
stimmeta<-stimmeta[!is.na(stimmeta$TRB_junction),]
unstimmeta<-unstimmeta[!is.na(unstimmeta$junction),]

#split by cluster in unstim but antigen specificity in stim
stimmeta<-split(stimmeta, stimmeta$antigen2)
unstimmeta<-split(unstimmeta, unstimmeta$cluster_names)

#what cluster in unstim meta has the most overlap 
results<-list()
for(i in names(unstimmeta)){
  cur<-c()
  for(j in names(stimmeta)){
    cur<-c(cur, sum(unstimmeta[[i]]$junction%in%stimmeta[[j]]$TRB_junction)/nrow(unstimmeta[[i]]))
  }
  results[[i]]<-cur
}
names(results)<-names(unstimmeta)
results<-unlist(results, recursive = FALSE)

results<-as.data.frame(matrix(data = results, ncol  = 10))
rownames(results)<-names(stimmeta)
colnames(results)<-names(unstimmeta)
pdf("Figures/Fig3stuff/antigentocluster")
heatmap(as.matrix(results), scale="row")
dev.off()
