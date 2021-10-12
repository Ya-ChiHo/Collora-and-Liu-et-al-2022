library(Seurat)
library(dplyr)
library(ggplot2)

setwd("~/project/JointAnalyses/")

pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

#reading in files, annotating HIV+ clones, cleaning and wrangling

unstim<-readRDS("Objects/20210519UnstimulatedMaster.rds")
stim<-readRDS("Objects/210623masterstimulated.rds")

X<-readRDS("Objects/20210601CellsincloneswithHIV.rds")
stim$HIV<-colnames(stim)%in%X
unstim$HIV<-colnames(unstim)%in%X

metau<-unstim@meta.data
metas<-stim@meta.data

metau$HIV<-factor(metau$HIV, levels = c(TRUE,FALSE))
metas$HIV<-factor(metas$HIV, levels = c(TRUE,FALSE))
metas$cluster_annotation<-factor(metas$cluster_annotation, levels=c("CD4-GZMB Th1","CD4-IL2 Th1","CD4-Lymphotoxin",
                                                                    "CD4-CD45RA","CD4-HLA-DR","CD4-GZMH","CD4-Treg",
                                                                    "CD4-IFN","CD4-Memory-1","CD4-Memory-2","CD4-Memory-3",
                                                                    "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-HSP"))
metas<-metas[metas$stage=="Viremic"|metas$stage=="Suppressed",]
metau<-metau[metau$infection!="Uninfected",]

#plotting antigen specificity differences
pubfig169p("Figures/Fig6stuff/specificity")
metas%>%group_by( stage,HIV,antigen2)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIV,y=n, fill=antigen2)) + geom_bar(position="fill", stat = "identity") +
  theme_classic() +facet_wrap(facets=.~ stage)+
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Category",
       x = NULL,
       y = NULL,
       title = "Proportion of cluster") +
  coord_polar(theta="y")
dev.off()
#plotting memory differences
pubfig169p("Figures/Fig6stuff/mem_stimgrouped")
metas%>%group_by(stage,HIV,memory)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIV,y=n, fill=memory))+
  facet_wrap(facets=.~stage ) + geom_bar(position="fill", stat = "identity") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Category",
       x = NULL,
       y = NULL,
       title = "Proportion of cluster") +
  coord_polar(theta="y")
dev.off()

pubfig169p("Figures/Fig6stuff/mem_unstimgrouped")
metau%>%group_by(infection, HIV,memory)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIV,y=n, fill=memory))+
  facet_wrap(facets=.~ infection) + geom_bar(position="fill", stat = "identity") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Category",
       x = NULL,
       y = NULL,
       title = "Proportion of cluster") +
  coord_polar(theta="y")

dev.off()

#plotting cluster differences
pubfig169p("Figures/Fig6stuff/clusters_stimgrouped")
metas%>%group_by(stage, HIV,cluster_annotation)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIV,y=n, fill=cluster_annotation))+
  facet_wrap(facets=.~ stage) + geom_bar(position="fill", stat = "identity") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Category",
       x = NULL,
       y = NULL,
       title = "Proportion of cluster") +
  coord_polar(theta="y")

dev.off()
pubfig169p("Figures/Fig6stuff/clusters_unstimgrouped")
metau%>%group_by(infection, HIV,cluster_names)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIV,y=n, fill=cluster_names))+
  facet_wrap(facets=.~ infection) + geom_bar(position="fill", stat = "identity") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Category",
       x = NULL,
       y = NULL,
       title = "Proportion of cluster") +
  coord_polar(theta="y")
dev.off()

