library(Seurat)
library(cvequality)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)
#reading in files and clones  

setwd("~/project/JointAnalyses/")

HIV<-readRDS("Objects/20210601CellsincloneswithHIV.rds")
stim<-readRDS("Objects/210623masterstimulated.rds")
unstim<-readRDS("Objects/20210519UnstimulatedMaster.rds")


#HIV RNA+ clones
stim$HIV<-colnames(stim)%in%HIV
unstim$HIV<-colnames(unstim)%in%HIV

stim$keep<-stim$bulk_n>1&stim$supprstage!="Uninfected"
stim<-subset(stim, keep==TRUE)
unstim$keep<-unstim$results_count>1&unstim$infection!="Uninfected"
unstim<-subset(unstim, keep==TRUE)

pubfig169p("Figures/Fig5stuff/stimRNAclones")
DimPlot(stim,cells.highlight = grep(TRUE,stim$HIV))+NoLegend()
dev.off()

pubfig169p("Figures/Fig5stuff/unstimRNAclones")
DimPlot(unstim,cells.highlight = grep(TRUE,unstim$HIV), raste=FALSE)+NoLegend()
dev.off()




meta<-readRDS("Objects/metadata with both clone timepoint quantifications")
#using the updated two timpoint frequency, removing non clonal cells
meta$HIV<-rownames(meta)%in%HIV
meta$HIV<-factor(meta$HIV, levels = c(TRUE, FALSE))
meta<-meta[meta$bulk>1,]
minsize<-min(unlist(c(meta$bulkfreq_1),meta$bulkfreq_2)[unlist(c(meta$bulkfreq_1),meta$bulkfreq_2)>0])
maxsize<-max(unlist(c(meta$bulkfreq_1),meta$bulkfreq_2)[unlist(c(meta$bulkfreq_1),meta$bulkfreq_2)>0])
#plotting for viremic and suppresed timepoint each condition to make 6A-D
CMV<-meta[meta$antigen2=="CMV",]
CMV<-group_by(CMV, HIV,UID,junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
p1<-ggplot(CMV, aes(x=HIV, fill=HIV, y=bulkfreq_1))+geom_violin()+scale_y_log10(limits=c(minsize, maxsize))+theme_classic()
p2<-ggplot(CMV, aes(x=HIV, fill=HIV, y=bulkfreq_2))+geom_violin()+scale_y_log10(limits=c(minsize, maxsize))+theme_classic()
pubfigsqp("Figures/Fig5stuff/CMVclonesize")
p1+p2
dev.off()

HIV<-meta[meta$antigen2=="HIV",]
HIV<-group_by(HIV, HIV,UID,junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
p1<-ggplot(HIV, aes(x=HIV, fill=HIV, y=bulkfreq_1))+geom_violin()+scale_y_log10(limits=c(minsize, maxsize))+theme_classic()
p2<-ggplot(HIV, aes(x=HIV, fill=HIV, y=bulkfreq_2))+geom_violin()+theme_classic()+scale_y_log10(limits=c(minsize, maxsize))
pubfigsqp("Figures/Fig5stuff/HIVclonesize")
p1+p2
dev.off()


Memory<-meta[meta$antigen2=="Memory",]
Memory<-group_by(Memory, HIV,UID,junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
p1<-ggplot(Memory, aes(x=HIV, fill=HIV, y=bulkfreq_1))+geom_violin()+scale_y_log10(limits=c(minsize, maxsize))+theme_classic()
p2<-ggplot(Memory, aes(x=HIV, fill=HIV, y=bulkfreq_2))+geom_violin()+scale_y_log10(limits=c(minsize, maxsize))+theme_classic()
pubfigsqp("Figures/Fig5stuff/Memoryclonesize")
p1+p2
dev.off()

resting<-meta[meta$antigen2=="resting",]
resting<-group_by(resting, HIV,UID,junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
p1<-ggplot(resting, aes(x=HIV, fill=HIV, y=bulkfreq_1))+geom_violin()+scale_y_log10(limits=c(minsize, maxsize))+theme_classic()
p2<-ggplot(resting, aes(x=HIV, fill=HIV, y=bulkfreq_2))+geom_violin()+theme_classic()+scale_y_log10(limits=c(minsize, maxsize))
pubfigsqp("Figures/Fig5stuff/Unstimulatedclonesize")
p1+p2
dev.off()

