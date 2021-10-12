library(Seurat)
library(ggplot2)
library(dplyr)
setwd("~/project/JointAnalyses/")
#generate some key genes violins in HIV RNA+ vs negative in figure 7
pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

#reading in files 
samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
HIV<-readRDS("Objects/20210601CellsincloneswithHIV.rds")
samples$HIV<-colnames(samples)%in%HIV
table(samples$HIV)
samples$keep<-samples$infection!="Uninfected" & samples$results_count>1
samples<-subset(samples, keep==TRUE)
#selected some genes 
genes<-c("GZMB", "GZMH", "GZMA", "PRF1", "CCL5", "CST7", "NKG7", "IL2RG")
samples$HIV<-factor(samples$HIV, levels=c(TRUE,FALSE))
Idents(samples)<-samples$HIV

#plotted
pubfig169p("Figures/Fig7newfollowups/vlnplotwithoutdotswithboxandwiskar")
cur<-VlnPlot(samples, features = genes, pt.size = 0,combine = FALSE, same.y.lims = TRUE)
for(i in 1:length(cur)){
  cur[[i]]<-cur[[i]]+geom_boxplot(width=0.1)+NoLegend()
}
patchwork::wrap_plots(cur, ncol=4)
dev.off()
