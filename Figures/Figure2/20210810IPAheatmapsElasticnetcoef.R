library(ggplot2)
library(reshape2)
library(dplyr)

setwd("~/Documents/HoLab/JoinedAnalyses/IPAElasticnetresults/")

#unstim pathways
data<-read.table("unstimulatedcomparisonpathways.txt", sep="\t", header = 1)
colnames(data)<-gsub("unstimulated_meta_","",colnames(data))
colnames(data)<-gsub("results.{1,}","",colnames(data))
data<-data[,c("Canonical.Pathways", "Viremic","Suppressed","Uninfected")]

data<-melt(data)
orderpaths<-data$Canonical.Pathways[data$variable=="Suppressed"][order(data$value[data$variable=="Suppressed"], decreasing = TRUE)]
data$Canonical.Pathways<-factor(data$Canonical.Pathways, levels = rev(orderpaths))
pdf("Pathways_unstimulated.pdf")
ggplot(data[data$Canonical.Pathways%in%orderpaths[1:20],], aes(x=variable, y=Canonical.Pathways, fill=value))+geom_tile()+theme_classic()+theme(title = element_blank(),axis.ticks = element_blank())
dev.off()



#unstim upstream
data<-read.table("unstimulatedcomparisonupstream.txt", sep="\t", header = 1)
colnames(data)<-gsub("unstimulated_meta_","",colnames(data))
colnames(data)<-gsub("results.{1,}","",colnames(data))
data<-data[,c("Upstream.Regulators", "Viremic","Suppressed","Uninfected")]

data<-melt(data)
orderpaths<-data$Upstream.Regulators[data$variable=="Suppressed"][order(data$value[data$variable=="Suppressed"], decreasing = TRUE)]
data$Upstream.Regulators<-factor(data$Upstream.Regulators, levels = rev(orderpaths))
pdf("upstream_unstimulated.pdf")
ggplot(data[data$Upstream.Regulators%in%orderpaths[1:20],], aes(x=variable, y=Upstream.Regulators, fill=value))+geom_tile()+theme_classic()+theme(title = element_blank(),axis.ticks = element_blank())
dev.off()

