library(Seurat)
library(dplyr)
library(ggplot2)

setwd("~/project/JointAnalyses/")

pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}


unstim<-readRDS("Objects/20210519UnstimulatedMaster.rds")
stim<-readRDS("Objects/210623masterstimulated.rds")

metau<-unstim@meta.data
metas<-stim@meta.data

metau$HIV_infected<-factor(metau$HIV_infected, levels=c("infected","uninfected"))
metas$HIVpositive<-factor(metas$HIVpositive, levels=c(TRUE,FALSE))
metas$cluster_annotation<-factor(metas$cluster_annotation, levels=c("CD4-GZMB Th1","CD4-IL2 Th1","CD4-Lymphotoxin",
                                                                   "CD4-CD45RA","CD4-HLA-DR","CD4-GZMH","CD4-Treg",
                                                                   "CD4-IFN","CD4-Memory-1","CD4-Memory-2","CD4-Memory-3",
                                                                   "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-HSP"))

metas$stage<-case_when(metas$stage=="Suppressed"~"Suppressed", metas$supprstage=="Viremic"~"Viremic", T~"Uninfected")
#cluster distribution of HIV-1 RNA+ and RNA- cells in unstim and stim
pubfig169p("Figures/Fig4stuff/allPTClusterunstim")
metau%>%group_by(infection, HIV_infected,cluster_names)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIV_infected,y=n, fill=cluster_names))+facet_wrap(facets = .~infection)+
   geom_bar(position="fill", stat = "identity") +
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
pubfig169p("Figures/Fig4stuff/allPTClusterstim")
metas%>%group_by( stage,HIVpositive,cluster_annotation)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIVpositive,y=n, fill=cluster_annotation))+ facet_wrap(facets = .~stage)+
  geom_bar(position="fill", stat = "identity") +
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
#memory distribution of HIV-1 RNA+ and RNA- cells in unstim and stim
pubfig169p("Figures/Fig4stuff/allPTmemoryunstim")
metau%>%group_by(infection,HIV_infected,memory)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIV_infected,y=n, fill=memory))+facet_wrap(facets = .~infection)+
  geom_bar(position="fill", stat = "identity") +
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
pubfig169p("Figures/Fig4stuff/allPTMemorystim")
metas%>%group_by( stage,HIVpositive,memory)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIVpositive,y=n, fill=memory))+ facet_wrap(facets = .~stage)+
  geom_bar(position="fill", stat = "identity") +
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
#antigen distribution of HIV-1 RNA+ and RNA- cells in   stim

pubfig169p("Figures/Fig4stuff/allPTantigenstim")
metas%>%group_by( stage,HIVpositive,antigen2)%>%summarise(n=n())%>%mutate(n=n/sum(n))%>%
  ggplot(aes(x=HIVpositive,y=n, fill=antigen2))+ facet_wrap(facets = .~stage)+
  geom_bar(position="fill", stat = "identity") +
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

#cluster stats 


metau<-split(metau, metau$infection)
metas<-split(metas, metas$stage)
#testing unstim data
for(i in c("Viremic", "Suppressed")){
  totalHIV<-sum(metau[[i]]$HIV_infected=="infected")
  totalnotHIV<-sum(!metau[[i]]$HIV_infected=="infected")
  
  if (totalHIV==0){print("skip")
    next()}
  for(j in names(table(metau[[i]]$cluster_names))){
    print(i)
    print(j)
    HIVcluster<-sum(metau[[i]]$HIV_infected=="infected"&metau[[i]]$cluster_names==j)
    totalcluster<-sum(metau[[i]]$HIV_infected=="uninfected"&metau[[i]]$cluster_names==j)
    print(HIVcluster)
    print(totalcluster)
    print(fisher.test(rbind(c(HIVcluster, totalcluster),c(totalHIV-HIVcluster, totalnotHIV-HIVcluster))))
}
}
#testing stim data
for(i in c("Viremic", "Suppressed")){
  totalHIV<-sum(metas[[i]]$HIVpositive==TRUE)
  totalnotHIV<-sum(!metas[[i]]$HIVpositive==TRUE)
  
  if (totalHIV==0){print("skip")
    next()}
  for(j in names(table(metas[[i]]$cluster_annotation))){
    print(i)
    print(j)
    HIVcluster<-sum(metas[[i]]$HIVpositive==TRUE&metas[[i]]$cluster_annotation==j)
    totalcluster<-sum(metas[[i]]$HIVpositive==FALSE&metas[[i]]$cluster_annotation==j)
    print(HIVcluster)
    print(totalcluster)
    print(fisher.test(rbind(c(HIVcluster, totalcluster),c(totalHIV-HIVcluster, totalnotHIV-HIVcluster))))
  }
}

#memory distribution 
#unstim
for(i in c("Viremic", "Suppressed")){
  totalHIV<-sum(metau[[i]]$HIV_infected=="infected")
  totalnotHIV<-sum(!metau[[i]]$HIV_infected=="infected")
  
  if (totalHIV==0){print("skip")
    next()}
  for(j in names(table(metau[[i]]$memory))){
    print(i)
    print(j)
    HIVcluster<-sum(metau[[i]]$HIV_infected=="infected"&metau[[i]]$memory==j)
    totalcluster<-sum(metau[[i]]$HIV_infected=="uninfected"&metau[[i]]$memory==j)
    print(HIVcluster)
    print(totalcluster)
    print(fisher.test(rbind(c(HIVcluster, totalcluster),c(totalHIV-HIVcluster, totalnotHIV-HIVcluster))))
  }
}
#stim
for(i in c("Viremic", "Suppressed")){
  totalHIV<-sum(metas[[i]]$HIVpositive==TRUE)
  totalnotHIV<-sum(!metas[[i]]$HIVpositive==TRUE)
  
  if (totalHIV==0){print("skip")
    next()}
  for(j in names(table(metas[[i]]$memory))){
    print(i)
    print(j)
    HIVcluster<-sum(metas[[i]]$HIVpositive==TRUE&metas[[i]]$memory==j)
    totalcluster<-sum(metas[[i]]$HIVpositive==FALSE&metas[[i]]$memory==j)
    print(HIVcluster)
    print(totalcluster)
    print(fisher.test(rbind(c(HIVcluster, totalcluster),c(totalHIV-HIVcluster, totalnotHIV-HIVcluster))))
  }
}

#getting positive markers for enrichr plot
res<-readRDS("Objects/20210602markersofHIVinfectedcells.rds")
names(res)
res<-res$infected_resting_resting_Viremic_vs_uninfected_resting_resting_Viremic
res<-res[res$p_val_adj<0.05 & res$avg_log2FC>0,]
write.table(rownames(res), "Figures/Fig4stuff/HIVRNAviremicmarkers.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
