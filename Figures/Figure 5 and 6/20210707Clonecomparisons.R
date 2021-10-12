library(Seurat)
library(cvequality)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)

#goal of this script is to generate clonal comparisons of HIV infected and uninfected. This script is mainly related to figures 5 and 6. 

setwd("~/project/JointAnalyses/")
#reading in files, fixing some metadata for simplicity
samples<-readRDS("Objects/20210521BothDatasetsLite.rds")

samples$antigen2[is.na(samples$antigen2)]<-""

samples$antigen2<-case_when(samples$antigen2=="hdCMVresp"|samples$antigen2=="pCMVresp"~"CMV", samples$antigen2=="hdmemory" | samples$antigen2=="pmemory"~"Memory",samples$antigen2=="pHIVresp"~"HIV", TRUE~"resting")
samples$infection<-case_when(samples$infection2=="Viremic"|samples$supprstage=="Viremic" ~"Viremic", samples$infection2=="Suppressed_Delayed"|samples$supprstage=="DelayedSuppr"~"Delayed_Suppression",samples$infection2=="Suppressed_Immediate"|samples$supprstage=="ImmedSuppr"~"Immediate_Suppression", TRUE~"Uninfected")
samples$HIV_infected<-case_when(samples$HIV_infected=="infected" | samples$HIVpositive==TRUE ~ "infected", TRUE~"uninfected")
samples$infection2<-case_when(samples$infection2=="Viremic"|samples$supprstage=="Viremic" ~"Viremic", samples$infection2=="Suppressed_Delayed"|samples$supprstage=="DelayedSuppr"|samples$infection2=="Suppressed_Immediate"|samples$supprstage=="ImmedSuppr"~"Suppressed", TRUE~"Uninfected")
samples$HIVcondtion2<-paste(samples$HIV_infected, samples$condition, samples$antigen2, samples$infection2,sep="_")
samples<-subset(samples, infection2!="Uninfected")
samples$junction2<-samples$junction
samples$junction2[is.na(samples$junction2)]<-samples$TRB_junction[is.na(samples$junction2)]
samples$bulk<-samples$results_count
samples$bulk[is.na(samples$bulk)]<- samples$bulk_n[is.na(samples$bulk)]
#now doing the actual testing
meta<-samples@meta.data
meta$results_total<-meta$results_count/meta$results_freq 
meta$results_total[meta$results_total==meta$results_total[[1]]]<-0

#first need to make sure everyone gets their timepoint two annotated 

meta<-meta[!is.na(meta$junction),]

meta<-split(meta, meta$PT)

for(i in 1:length(meta)){
  current<-split(meta[[i]], meta[[i]]$UID)
  #need to identify matching clones in both and add in the n, total, and calculate frequency for the alternate time point
  #querry tp2 with tp1
  tp1<-current[[grep("_1", names(current))]]
  tp2<-current[[grep("_2", names(current))]]
  tp1res<-list()
  tp1total<-max(tp1$results_total,na.rm = TRUE)
  tp2total<-max(tp2$results_total, na.rm=TRUE)
  
  for(j in 1:nrow(tp1)){
    if(tp1$junction2[[j]] %in% tp2$junction2){
      hit<-grep(tp1$junction2[[j]], tp2$junction2)[[1]]
      tp1res[[j]]<-c(tp2$bulk[[hit]], tp2total,tp2$bulk[[hit]]/tp2total )
    }else{tp1res[[j]]<-c(0,tp2total,0)}
  }
  tp2res<-list()
  #querry tp 1 with tp2
  for(j in 1:nrow(tp2)){
    if(tp2$junction2[[j]] %in% tp1$junction2){
      hit<-grep(tp2$junction2[[j]], tp1$junction2)[[1]]
      tp2res[[j]]<-c(tp1$bulk[[hit]], tp1total,tp1$bulk[[hit]]/tp1total )
    }else{tp2res[[j]]<-c(0,tp1total, 0)}
  }
  #reshape, annotate and append to existing data 
  tp1res<-as.data.frame(matrix(unlist(tp1res), ncol=3, byrow = TRUE))
  tp2res<-as.data.frame(matrix(unlist(tp2res), ncol=3, byrow = TRUE))
  
  names(tp1res)<-c("bulkn_2","bulktotal_2","bulkfreq_2")
  names(tp2res)<-c("bulkn_1","bulktotal_1","bulkfreq_1")
  tp1<-cbind( tp1,tp1res)
  tp2<-cbind( tp2,tp2res)
  
  tp1$bulkn_1<-tp1$bulk
  tp1$bulktotal_1<-tp1total
  tp1$bulkfreq_1<-tp1$bulk/tp1total
  tp2$bulkn_2<-tp2$bulk
  tp2$bulktotal_2<-tp2total
  tp2$bulkfreq_2<-tp2$bulk/tp2total
  
  current<-rbind(tp1,tp2)
  meta[[i]]<-current
}

meta<-rbind(meta[[1]],meta[[2]],meta[[3]],meta[[4]],meta[[5]],meta[[6]])

saveRDS(meta, "Objects/metadata with both clone timepoint quantifications")
meta<-readRDS("Objects/metadata with both clone timepoint quantifications")
  #first complete a fisher test to see if there is a timepoint difference

fisherp<-c()
for(i in 1:nrow(meta)){
  fisherp[[i]]<-fisher.test(rbind(c(meta$bulkn_1[[i]], meta$bulktotal_1[[i]]),c(meta$bulkn_2[[i]],meta$bulkn_2[[i]])))$p.value
}
meta$fisherres<-fisherp
#90k are significant at this p cutoff, cells not clones
table(fisherp<0.05)
meta$change<-case_when(meta$fisherres<0.05 &(meta$bulkfreq_1>meta$bulkfreq_2)~"contracting",meta$fisherres<0.05 &(meta$bulkfreq_1<meta$bulkfreq_2)~"expanding", TRUE~"stable"  )

meta2<-meta%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
meta2<-meta2[meta2$bulk>1,]
#8.4k clones are different at this cutoff
table(meta2$fisherres<0.05)
#these remain signficant even with a BH correction
table(p.adjust(meta2$fisherres,"BH")<0.05)

#repeating the above with HIV as an additional axis of interest 
HIVcells<-readRDS("Objects/20210601CellsincloneswithHIV.rds")
meta$HIV <- rownames(meta)%in%HIVcells

#first we'll check if HIV infecrted clones are different
metaHIV<-meta[meta$HIV==TRUE,]

metanotHIV<-meta[meta$HIV==FALSE,]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
meta2<-rbind(metaHIV, metanotHIV)
meta2<-meta2[meta2$bulk>1,]
#55/68 clones are different at this cutoff (80%)
table(meta2$fisherres[meta2$HIV==TRUE]<0.05)
#35 (64%) of those 55 expand, 20 (36%) contract, and the last 13 are stable
table(meta2$change[meta2$HIV==TRUE])

#of 11841 clones 8.371k are different (70%)
table(meta2$fisherres[meta2$HIV==FALSE]<0.05)
#6792 (81%) of those 11841 expand, 1579 (19%) contract, and the last 3470 are stable
table(meta2$change[meta2$HIV==FALSE])

#now for each condition we're going to look at the clone fold change (Figure 6E)
#just resting
metaHIV<-meta[meta$HIV==TRUE&meta$antigen2=="resting",]

metanotHIV<-meta[meta$HIV==FALSE&meta$antigen2=="resting",]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
resting<-rbind(metaHIV, metanotHIV)
resting<-resting[resting$bulk>1,]

#just CMV
metaHIV<-meta[meta$HIV==TRUE&meta$antigen2=="CMV",]

metanotHIV<-meta[meta$HIV==FALSE&meta$antigen2=="CMV",]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
CMV<-rbind(metaHIV, metanotHIV)
CMV<-CMV[CMV$bulk>1,]

#just HIV
metaHIV<-meta[meta$HIV==TRUE&meta$antigen2=="HIV",]

metanotHIV<-meta[meta$HIV==FALSE&meta$antigen2=="HIV",]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
HIV<-rbind(metaHIV, metanotHIV)
HIV<-HIV[HIV$bulk>1,]

#just Memory
metaHIV<-meta[meta$HIV==TRUE&meta$antigen2=="Memory",]

metanotHIV<-meta[meta$HIV==FALSE&meta$antigen2=="Memory",]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
Memory<-rbind(metaHIV, metanotHIV)
Memory<-Memory[Memory$bulk>1,]

#need to make an adjustment for 0s in the log
resting$FC<-log2((resting$bulkfreq_2+min(resting$bulkfreq_2[resting$bulkfreq_2>0]))/(resting$bulkfreq_1+min(resting$bulkfreq_1[resting$bulkfreq_1>0])))
HIV$FC<-log2((HIV$bulkfreq_2+min(HIV$bulkfreq_2[HIV$bulkfreq_2>0]))/(HIV$bulkfreq_1+min(HIV$bulkfreq_1[HIV$bulkfreq_1>0])))
CMV$FC<-log2((CMV$bulkfreq_2+min(CMV$bulkfreq_2[CMV$bulkfreq_2>0]))/(CMV$bulkfreq_1+min(CMV$bulkfreq_1[CMV$bulkfreq_1>0])))
Memory$FC<-log2((Memory$bulkfreq_2+min(Memory$bulkfreq_2[Memory$bulkfreq_2>0]))/(Memory$bulkfreq_1+min(Memory$bulkfreq_1[Memory$bulkfreq_1>0])))


#wilcox tests to look for changes in clone size, every case is significantly different , Figure 6A-D

wilcox.test(resting$bulkfreq_1[resting$HIV==TRUE],resting$bulkfreq_1[resting$HIV!=TRUE])
wilcox.test(HIV$bulkfreq_1[HIV$HIV==TRUE],HIV$bulkfreq_1[HIV$HIV!=TRUE])
wilcox.test(CMV$bulkfreq_1[CMV$HIV==TRUE],CMV$bulkfreq_1[CMV$HIV!=TRUE])
wilcox.test(Memory$bulkfreq_1[Memory$HIV==TRUE],Memory$bulkfreq_1[Memory$HIV!=TRUE])

wilcox.test(resting$bulkfreq_2[resting$HIV==TRUE],resting$bulkfreq_2[resting$HIV!=TRUE])
wilcox.test(HIV$bulkfreq_2[HIV$HIV==TRUE],HIV$bulkfreq_2[HIV$HIV!=TRUE])
wilcox.test(CMV$bulkfreq_2[CMV$HIV==TRUE],CMV$bulkfreq_2[CMV$HIV!=TRUE])
wilcox.test(Memory$bulkfreq_2[Memory$HIV==TRUE],Memory$bulkfreq_2[Memory$HIV!=TRUE])

####clonesize HIV vs not#### 

resting<-resting[,c("HIV","FC")]
resting$ID<-"Resting"
Memory<-Memory[,c("HIV","FC")]
Memory$ID<-"Memory"
CMV<-CMV[,c("HIV","FC")]
CMV$ID<-"CMV"
HIV<-HIV[,c("HIV","FC")]
HIV$ID<-"HIV"
data<-rbind(resting, Memory, CMV, HIV)

#figure 6E
p1<-ggplot(data[data$HIV==TRUE,], aes(x=FC, color=ID))+geom_density()+theme_classic()+ylim(0,0.7)+xlim(-10,10)
p2<-ggplot(data[data$HIV==FALSE,], aes(x=FC, color=ID))+geom_density()+theme_classic()+ylim(0,0.7)+xlim(-10,10)
####clonaldynamics hiv vs not####
pubfig169p("Figures/Clonalfigs/Clonaldynamicscomparison")
p1+p2
dev.off()



####are HIV containing clones more persistent####
#related to figure 6F-J
{
  #first ensure we're just to clones, then categorize them
  meta2<-meta[!is.na(meta$bulk),]
  meta2<-meta2[meta2$bulk>1,]
  meta2$lost<-case_when(meta2$bulkn_1==0 & meta2$bulkn_2>0 ~"Suppression only",meta2$bulkn_1>0&meta2$bulkn_2==0~"Viremia only",T~"both")
  table(meta2$lost[meta2$bulk>1], meta2$HIV[meta2$bulk>1])
  table(meta2$lost[meta2$bulk>1], meta2$antigen2[meta2$bulk>1])
#just CMV
metaHIV<-meta2[meta2$HIV==TRUE&meta2$antigen2=="CMV",]

metanotHIV<-meta2[meta2$HIV==FALSE&meta2$antigen2=="CMV",]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
CMV<-rbind(metaHIV, metanotHIV)
CMV<-cbind(table(CMV$lost[CMV$HIV=="TRUE"]), table(CMV$lost[CMV$HIV!="TRUE"]))
CMV<-melt(CMV)
CMV$ID<-"CMV"

#just HIV
metaHIV<-meta2[meta2$HIV==TRUE&meta2$antigen2=="HIV",]

metanotHIV<-meta2[meta2$HIV==FALSE&meta2$antigen2=="HIV",]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
HIV<-rbind(metaHIV, metanotHIV)
HIV<-cbind(table(HIV$lost[HIV$HIV=="TRUE"]), table(HIV$lost[HIV$HIV!="TRUE"]))
HIV<-melt(HIV)
HIV$ID<-"HIV"

#just Memory
metaHIV<-meta2[meta2$HIV==TRUE&meta2$antigen2=="Memory",]

metanotHIV<-meta2[meta2$HIV==FALSE&meta2$antigen2=="Memory",]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
Memory<-rbind(metaHIV, metanotHIV)
Memory<-cbind(table(Memory$lost[Memory$HIV=="TRUE"]), table(Memory$lost[Memory$HIV!="TRUE"]))
Memory<-melt(Memory)
Memory$ID<-"Memory"

#just resting
metaHIV<-meta2[meta2$HIV==TRUE&meta2$antigen2=="resting",]

metanotHIV<-meta2[meta2$HIV==FALSE&meta2$antigen2=="resting",]
metaHIV<-metaHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
metanotHIV<-metanotHIV%>%group_by( PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
resting<-rbind(metaHIV, metanotHIV)
resting<-cbind(table(resting$lost[resting$HIV=="TRUE"]), table(resting$lost[resting$HIV!="TRUE"]))
resting<-melt(resting)
resting$ID<-"resting"

#testing if the props are different between HIV RNA+ and RNA- clones 
fisher.test(rbind(resting$value[resting$Var2==1],resting$value[resting$Var2==2]))
fisher.test(rbind(CMV$value[CMV$Var2==1],CMV$value[CMV$Var2==2]))
fisher.test(rbind(HIV$value[HIV$Var2==1],HIV$value[HIV$Var2==2]))
fisher.test(rbind(Memory$value[Memory$Var2==1],Memory$value[Memory$Var2==2]))
data<-rbind(resting, CMV, HIV, Memory)

p1<-ggplot(data[data$ID=="CMV",], aes(x=Var2, y=value, fill=Var1))+geom_bar(stat="identity", position="fill")+theme_classic()+ylab("Proportion")+xlab("")+ggtitle("Clones Harboring HIV-1")
p2<-ggplot(data[data$ID=="HIV",], aes(x=Var2, y=value, fill=Var1))+geom_bar(stat="identity", position="fill")+theme_classic()+ylab("Proportion")+xlab("")+ggtitle("Clones Harboring HIV-1")
p3<-ggplot(data[data$ID=="Memory",], aes(x=Var2, y=value, fill=Var1))+geom_bar(stat="identity", position="fill")+theme_classic()+ylab("Proportion")+xlab("")+ggtitle("Clones Harboring HIV-1")
p4<-ggplot(data[data$ID=="resting",], aes(x=Var2, y=value, fill=Var1))+geom_bar(stat="identity", position="fill")+theme_classic()+ylab("Proportion")+xlab("")+ggtitle("Clones Harboring HIV-1")
library(ggpubr)
#figure 6F
pubfig169p("Figures/Clonalfigs/ClonescontainingHIVaremorepersistent")
ggarrange(p1,p2,p3,p4, nrow=1)
dev.off()
}

#starting the top 10/100 clones and highlighting them based on their origin

HIVcells<-readRDS("Objects/20210601CellsincloneswithHIV.rds")
meta$HIV <- rownames(meta)%in%HIVcells

agpromo <- function(cat){
  print(cat)
  if("CMV"%in%names(table(cat))|"HIV"%in%names(table(cat))){
    
    cat<-cat[cat=="HIV"|cat=="CMV"]
    print(names(table(cat))[which.max(table(cat))])

  }else{
  print(names(table(cat))[which.max(table(cat))])
  }}


#first just find the top 100 clones
#unique only
meta2<-group_by(meta, PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
#clones only
meta2<-meta2[meta2$bulk>1,]
#top 100
meta2<-ungroup(meta2)%>%group_by(PT)%>%slice_max(order_by=bulkfreq_2,n=100, with_ties=FALSE)
#generate the PT junction pairs 
meta2<-meta2[,c("PT","junction2")]
meta2<-split(meta2, meta2$PT)
metaPTsplit<-split(meta,meta$PT )

#for each patient grab the junctions related 
for(i in 1:length(meta2)){
  cur<-meta2[[i]]
  res<-metaPTsplit[[i]][metaPTsplit[[i]]$junction2==cur$junction2[[1]],]
  for(j in 2:nrow(cur)){
    res<-rbind(res,metaPTsplit[[i]][metaPTsplit[[i]]$junction2==cur$junction2[[j]],] )
  }
  metaPTsplit[[i]]<-res
}

meta2<-bind_rows(metaPTsplit)
meta2<-group_by(meta2, PT, junction2)%>%mutate(cond=agpromo(antigen2))%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)%>%ungroup()%>%group_by(PT)%>%mutate(freq=bulkfreq_2/sum(bulkfreq_2))
meta2<-meta2[order(meta2$freq, decreasing = TRUE),]
meta2$HIV[meta2$HIV==FALSE]<-NA
meta2$HIV[meta2$HIV==TRUE]<-"*"

ggplot(meta2, aes(x=PT, y=freq,fill=cond, group=freq,color="grey"))+geom_bar(stat="identity",size=.1)+scale_color_manual(values = c("grey"="black"))+geom_text(aes(label=HIV), color="black", size=3, position = position_stack(vjust=0.5))+
  theme_classic()+scale_fill_manual(values=c( "CMV"="#DF536B","HIV"="#2297E6","Memory"="#61D04F", "resting"="#F5C710"))

meta2$cond[meta2$cond=="Memory"|meta2$cond=="resting"]<-"antigen unknown"

metas<-meta2
metas$PT<-paste(metas$PT, "_2",sep="")

#viremia
meta2<-group_by(meta, PT, junction2)%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)
#clones only
meta2<-meta2[meta2$bulk>1,]
#top 100
meta2<-ungroup(meta2)%>%group_by(PT)%>%slice_max(order_by=bulkfreq_1,n=100, with_ties=FALSE)
#generate the PT junction pairs 
meta2<-meta2[,c("PT","junction2")]
meta2<-split(meta2, meta2$PT)
metaPTsplit<-split(meta,meta$PT )

#for each patient grab the junctions related 
for(i in 1:length(meta2)){
  cur<-meta2[[i]]
  res<-metaPTsplit[[i]][metaPTsplit[[i]]$junction2==cur$junction2[[1]],]
  for(j in 2:nrow(cur)){
    res<-rbind(res,metaPTsplit[[i]][metaPTsplit[[i]]$junction2==cur$junction2[[j]],] )
  }
  metaPTsplit[[i]]<-res
}

meta2<-bind_rows(metaPTsplit)
meta2<-group_by(meta2, PT, junction2)%>%mutate(cond=agpromo(antigen2))%>%slice_max(order_by=bulk, n=1, with_ties=FALSE)%>%ungroup()%>%group_by(PT)%>%mutate(freq=bulkfreq_1/sum(bulkfreq_1))
meta2<-meta2[order(meta2$freq, decreasing = TRUE),]
meta2$HIV[meta2$HIV==FALSE]<-NA
meta2$HIV[meta2$HIV==TRUE]<-"*"
meta2$cond[meta2$cond=="Memory"|meta2$cond=="resting"]<-"antigen unknown"


####combo plot, Figure 5D####
meta2$PT<-paste(meta2$PT, "_1", sep="")
meta2<-rbind(meta2, metas)

pubfig169p("Figures/Clonalfigs/Top100clones_both")
ggplot(meta2, aes(x=PT, y=freq,fill=cond, group=freq,color="grey"))+geom_bar(stat="identity",size=.1)+scale_color_manual(values = c("grey"="black"))+geom_text(aes(label=HIV), color="black", size=3, position = position_stack(vjust=0.5))+
  theme_classic()+scale_fill_manual(values=c( "CMV"="#2297E6","HIV"="#DF536B","antigen unknown"="#61D04F"))
dev.off()

####start of visualizing where certain things are, Figure 6 G-J####


#first with the unstimulated dataset, I want to know where the clones that change where we can find them vs the clones we do not are

#first need to rewrite things so that all of those clones get identified
meta2<-meta[!is.na(meta$bulk),]
meta2<-meta2[meta2$bulk>1,]
meta2$lost<-case_when(meta2$bulkn_1==0 & meta2$bulkn_2>0 ~"Suppression only",meta2$bulkn_1>0&meta2$bulkn_2==0~"Viremia only",T~"both")

#next we'll add that metadata to my master seurat object

samples<-readRDS("Objects/20210519UnstimulatedMaster.rds")
order_of_clusters<-levels(samples$cluster_names)
samples<-AddMetaData(samples, meta2)

#we will subset to just those that have an annotation ie are clonal

samples$keep<-!is.na(samples$lost)
samples<-subset(samples, keep==TRUE)

#12k cells remain, 6G
pubfig169p("Figures/Clonalfigs/Unstimulatedperistentcells_groupislost_splitisHIV")
DimPlot(samples, group.by = "lost", split.by = "HIV")
dev.off()

#Bars of peristence
samples$cluster_names<-factor(samples$cluster_names, levels=order_of_clusters)
samples$HIV<-case_when(samples$HIV==TRUE ~ "HIV RNA+", T~"HIV RNA-")%>%factor(levels = c("HIV RNA+","HIV RNA-"))
meta<-group_by(samples@meta.data, lost,HIV, cluster_names )%>%summarise(n=n(), .drop=FALSE)

p1<-ggplot(meta[meta$lost=="both",], aes(x=HIV, y=n, fill=cluster_names))+geom_bar(position="fill",stat="identity")+theme_classic()+ggtitle("cells in clones that persist")
p2<-ggplot(meta[meta$lost=="Suppression only",], aes(x=HIV, y=n, fill=cluster_names))+geom_bar(position="fill",stat="identity")+theme_classic()+ggtitle("cells in clones that are only seen during suppression")
p3<-ggplot(meta[meta$lost=="Viremia only",], aes(x=HIV, y=n, fill=cluster_names))+geom_bar(position="fill",stat="identity")+theme_classic()+ggtitle("cells in clones that are only seen during viremia")
#6H
pubfig169p("Figures/Clonalfigs/Phenotypeofcellsinclonesthatperist")
p1+p3+p2
dev.off()

#now were going todo the same with runxia's data
samples<-readRDS("Objects/210623masterstimulated.rds")

samples<-AddMetaData(samples, meta2)

samples$keep<-!is.na(samples$lost)
samples<-subset(samples, keep==TRUE)
#6I
pubfig169p("Figures/Clonalfigs/stimulatedperistentcells_groupislost_splitisHIV")
DimPlot(samples, group.by = "lost", split.by = "HIV")
dev.off()
meta<-group_by(samples@meta.data, lost,HIV, cluster_annotation )%>%summarise(n=n(), .drop=FALSE)
meta$cluster_annotation<-factor(meta$cluster_annotation, levels=c("CD4-GZMB Th1","CD4-IL2 Th1","CD4-Lymphotoxin",
                                                                  "CD4-CD45RA","CD4-HLA-DR","CD4-GZMH","CD4-Treg",
                                                                  "CD4-IFN","CD4-Memory-1","CD4-Memory-2","CD4-Memory-3",
                                                                  "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-HSP"))
meta$HIV<-case_when(meta$HIV==TRUE~"HIV RNA+",T~"HIV RNA-")%>%factor(levels=c("HIV RNA+","HIV RNA-"))
p1<-ggplot(meta[meta$lost=="both",], aes(x=HIV, y=n, fill=cluster_annotation))+geom_bar(position="fill",stat="identity")+theme_classic()+ggtitle("cells in clones that persist")
p2<-ggplot(meta[meta$lost=="Suppression only",], aes(x=HIV, y=n, fill=cluster_annotation))+geom_bar(position="fill",stat="identity")+theme_classic()+ggtitle("cells in clones that are only seen during suppression")
p3<-ggplot(meta[meta$lost=="Viremia only",], aes(x=HIV, y=n, fill=cluster_annotation))+geom_bar(position="fill",stat="identity")+theme_classic()+ggtitle("cells in clones that are only seen during viremia")
#6J
pubfig169p("Figures/Clonalfigs/Phenotypeofcellsinclonesthatperiststimulated")
p1+p3+p2
dev.off()

