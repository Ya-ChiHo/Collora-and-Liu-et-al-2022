library(Seurat)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tibble)
#this script takes the clonal correlation and runs the GSEA ranking, no figures are generated
setwd("~/project/JointAnalyses/")

#genesets to test
m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C2"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C6"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C7"), m_df_H)


fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}


GSEA<-function(corobj, genesets){
  genes<-sort(corobj["clonesize",],decreasing=TRUE)
  genes<-fgsea(pathways=genesets, stats=genes, eps=FALSE)
  genes<-list(genes,sort(corobj["clonesize",],decreasing=TRUE))
  names(genes)<-c("results","rnk")
  return(genes)}

todo<-list.files("Objects/cor/")
todo<-grep("HD", todo, value=TRUE, invert = TRUE)
res<-list()
for(i in todo){
  res2<-try(GSEA(readRDS(paste("Objects/cor/", i, sep="")),fgsea_sets))
  if(inherits(res2, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    print(paste(i, "failed"))
    next
  }
  res[[i]]<-res2
}
names(res)<-gsub(".rds","",todo)
saveRDS(res, "Objects/cor/finalGSEA.rds")

todo<-list.files("Objects/cor/")
todo<-grep("HD", todo, value=TRUE, invert = FALSE)
res<-list()
for(i in todo){
  res2<-try(GSEA(readRDS(paste("Objects/cor/", i, sep="")),fgsea_sets))
  if(inherits(res2, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    print(paste(i, "failed"))
    next
  }
  res[[i]]<-res2
}
names(res)<-gsub(".rds","",names(res))
saveRDS(res, "Objects/cor/HDGSEA.rds")
