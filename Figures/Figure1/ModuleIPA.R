library(reshape2)
library(dplyr)
library(ggplot2)

setwd("~/Documents/HoLab/JoinedAnalyses/Modules/")
#read in pathways, wrangle, pick the top n, removing header
x<-read.table(file = "Mod4_pathways.txt", sep = "\t")
colnames(x)<-c("pathways","pval","ratio","crap","crap","crap")
x<-x[2:11,1:2]
x$pval<-as.numeric(x$pval)
x$pathways<-factor(x$pathways, levels=rev(x$pathways))
pdf("module4pathways.pdf", width = 5, height = 8)
ggplot(x, aes(x = 1, y=pathways, fill=pval))+geom_tile()+theme_classic()+theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
dev.off()
#read in upstream regulators, wrangle, pick the top n, removing header, remove chemicals

x<-read.table(file = "Mod4_upstream.txt", sep = "\t")
colnames(x)<-c("regulator","mole","pval","crap","crap")
x<-x[2:100,1:3]
x$pval<- (-log10(as.numeric(x$pval)))
x$pathways<-factor(x$regulator, levels=rev(x$regulator))
x<-x[!grepl("chemical",x$mole),]
keep<-c("GPER1","TGFB1", "CSF2","TNF","RELA","TCR","IL1B","Jnk","IFNG","IGF1","IRAK4","F2","FADD","NR3C1")
pdf("module4upstream.pdf", width = 5, height = 8)
ggplot(x[x$regulator%in%keep,], aes(x = 1, y=pathways, fill=pval))+geom_tile()+theme_classic()+theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())+scale_color_continuous(limits=)
dev.off()






