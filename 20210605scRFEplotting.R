library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(SeuratDisk)
#utilities
{pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

VolPlot<-function(results,top=TRUE, adthresh=TRUE, thresh=0.2, threshn=30, Title=NA , pvalcutoff=300){
  results$value<-(-log10(results$p_val_adj))
  results$value[results$value==Inf | results$value>pvalcutoff]<-pvalcutoff
  results$gene<-rownames(results)
  results$color<-case_when(results$avg_log2FC>0.25 & results$value>(-log10(0.05))~"upregulated",results$avg_log2FC<(-0.25) & results$value>(-log10(0.05))~"downregulated", T~"notsignificant" )
  p<-ggplot(results,aes(x=avg_log2FC, y=value, color=color))+geom_point()+theme_classic()+scale_color_manual(values = c("notsignificant"="grey62","upregulated"="red", "downregulated"="blue"))
  if(adthresh==TRUE){thresh<-AdaptiveThreshold(results, Tstart=thresh,n=threshn)}
  if(top==TRUE){
    p<- p + geom_text_repel(aes(label=ifelse((avg_log2FC>thresh|avg_log2FC<(-thresh))&value>(-log10(0.05)) ,as.character(gene),'')),hjust=1,vjust=1, max.overlaps = 30000)
  }
  
  if(!is.na(Title)){
    p<-p+ggtitle(Title)+xlab("Average Log Fold Change")+ylab("-log10 adjusted p value")+theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}

AdaptiveThreshold<-function(results, n=30,Tstart=0.25){
  m=Inf
  while(n<=m){
    Tstart=Tstart+0.05
    m<-length(rownames(results)[(results$avg_log2FC>Tstart|results$avg_log2FC<(-Tstart))&results$value>10])
  }
  return(Tstart)
}}

#read in samples

setwd("~/project/JointAnalyses/")

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
samples<-subset(samples, infection!="Uninfected")

#add annotations for the "type" of clone
HIVclone<-readRDS("Objects/20210601CellsincloneswithHIV.rds")
#these are the stimulated clones
Runxiaclone<-read.csv("sklearn/ClonesHIVRunxiaData.csv", header = FALSE)
#these are the unstimulated clones
Jackclone<-read.csv("sklearn/ClonesHIVJackData.csv", header = FALSE)
#these are the results of the model validation (correct or incorrect) for stimulated clones
Runxiares<-read.csv("sklearn/runxiaresultsnew.txt", header = FALSE)
Runxiares<-gsub("rue|alse","",Runxiares[[1]])
Runxiares<-strsplit(Runxiares, "")[[1]]
Runxiaclone<-as.data.frame(Runxiaclone)
Runxiaclone$res<-Runxiares
rownames(Runxiaclone)<-Runxiaclone$V1
Runxiaclone["res"]
#these are the results of the model validation (correct or incorrect) for unstimulated clones
Jackres<-read.csv("sklearn/Jackresultsnew.txt",header = FALSE)
Jackres<-gsub("rue|alse","",Jackres[[1]])
Jackres<-strsplit(Jackres, "")[[1]]

Jacktrue<-read.csv("sklearn/JackTruthnew.txt")
Jacktrue$prediction<-Jackres
rownames(Jacktrue)<-Jacktrue$X
Jacktrue<-Jacktrue[,c("HIVclone", "prediction")]

samples$HIV_clone<-colnames(samples)%in%HIVclone
samples<-AddMetaData(samples, Runxiaclone["res"])
samples<-AddMetaData(samples, Jacktrue)
####umap of all clones####

#first subset to just clones

samples2<-subset(samples, results_count>1)


#umap of correctly vs incorrectly identified as "HIV infected"
samples2$prediction[is.na(samples2$prediction)]<-0
samples2$correct<-case_when(samples2$res=="T"~"correct",samples2$res=="F"~"incorrect", samples2$prediction=="T"&samples2$HIVclone=="True"~"correct", samples2$prediction=="F"&samples2$HIV_clone==TRUE~"incorrect", samples2$HIV_clone==TRUE~"trainingset",T~"notclone")
data$correct<-samples2$correct
pubfig169p("Figures/scRFEfigs/HIVclonesaccuracysplit")

ggplot(data, aes(x=UMAP_1, y=UMAP_2, color="NotHIVclone"))+geom_point()+geom_point(data = data[data$HIV==TRUE,],
                                                                                   mapping=aes(x=UMAP_1,y=UMAP_2, color=correct))+
  scale_color_manual(values = c("NotHIVclone"="grey62","FALSE"="grey62", "incorrect"="red", "correct"="blue", "trainingset"="green"))+theme_classic()
dev.off()

#differential expression volcano
samples2$HIV_clone<-factor(samples2$HIV_clone, c(TRUE,FALSE))
Idents(samples2)<-samples2$HIV_clone
markers<-FindMarkers(samples2, TRUE,FALSE)
pubfig169p("Figures/scRFEfigs/Alldifferentiallyexpressednoclustermatching")
VolPlot(markers, threshn = 100, Title = "Alldifferentiallyexpressednoclustermatching")
dev.off()
#Select key genes to show via violins 
pubfig169p("Figures/scRFEfigs/Differentiallyexpressedviolinsnoclustermatching")
VlnPlot(samples2, c("CST7","CCL4","CCL5","NKG7","GZMH","ZEB2","GZMB","GZMA","PRF1","GZMK","IL2RG","HLA-DRB5"))
dev.off()

####same but we just do the testing for the top 200 genes that we validated the model on
#this is all of the scRFE genes ranked by proportion in top 100 genes of the model
x<-read.csv("sklearn/20210605top100scRFE.csv")

#selecting the top 200 genes 
genes<-x$genes[order(x$Freq, decreasing = TRUE)[1:201]]
genes<-genes[!genes==""]
samples2$HIV_clone<-case_when(samples2$HIVclone=="nonHIVclone"~FALSE, T~TRUE)
samples2$HIV_clone<-factor(samples2$HIV_clone, c(TRUE,FALSE))
Idents(samples2)<-samples2$HIV_clone
#DE testing
markers<-FindMarkers(samples2, TRUE,FALSE, features =genes, logfc.threshold = -Inf, min.pct = -Inf, min.diff.pct = -Inf)
genes<-x[order(x$Freq, decreasing = TRUE)[1:201],]
genes<-genes[!genes$genes=="",c(2,3)]
rownames(genes)<-genes$genes
markers$genes<-rownames(markers)
cur<-merge(genes, markers, by="genes")

#saving and plotting
write.table(cur, "Supp/scRFEgenetable.tsv", quote=FALSE, row.names=FALSE, sep="\t")
pubfig169p("Figures/scRFEfigs/keydifferentiallyexpressednoclustermatching")
VolPlot(markers, threshn = 1000, Title = "keydifferentiallyexpressednoclustermatching")
dev.off()
