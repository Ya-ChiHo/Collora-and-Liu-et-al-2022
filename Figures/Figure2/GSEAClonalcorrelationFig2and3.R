library(scales)
library(ggplot2)
library(fgsea)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(Seurat)
library(ggrepel)
library(tidyr)
library(reshape2)
library(msigdbr)
library(tibble)

#goal of this script is to run the GSEA on clonal-gene correlation ranking, starting with the given genes and then going through a bunch of helper functions for plotting
m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C2"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C6"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C7"), m_df_H)
fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)

#Function 1 augments the existing result table with pathway information, including the entire rank order list of gene hits and NES scores 
GSEATable<-function(GSEAwrap_out,gmt, gseaParam=1, name){
  #doing constants first since those are easy
  finalresults<-GSEAwrap_out[1][[1]]
  #first get the order of genes for each pathway, modified from the fGSEA plot function  
  rnk <- rank(GSEAwrap_out[2][[1]])
  ord <- order(rnk, decreasing = T)
  statsAdj <- GSEAwrap_out[2][[1]][ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  pathwaysetup<-list()
  NES<-list()
  print("start ranking")
  #basically here I just want to get everything into a character vector
  for (i in 1:length(finalresults$pathway)){
    pathway <- unname(as.vector(na.omit(match(gmt[[finalresults$pathway[[i]]]], names(statsAdj)))))
    pathway <- sort(pathway)
    NES[[i]]<-calcGseaStat(statsAdj, selectedStats = pathway,returnAllExtremes = TRUE)
    NES[[i]]<-NES[[i]]$tops
    pathwaysetup[[i]]<-pathway}
  print("done ranking")
  finalresults$rnkorder<-pathwaysetup
  finalresults$NESorder<-NES
  finalresults$ID<-name
  return(finalresults)
}


GSEAbig<-function(listofGSEAtables){
  GSEA<-bind_rows(listofGSEAtables)
  return(GSEA)
}



GSEAEnrichmentPlotComparison<-function(PathwayName, GSEACompOut, returnplot="none"){
  top<-max(unlist(GSEACompOut$rnkorder))
  GSEACompOut<-GSEACompOut[GSEACompOut$pathway==PathwayName,]
  curves<-list()
  for (i in 1:length(rownames(GSEACompOut))){curves[[i]]<-list(c(0,GSEACompOut$rnkorder[[i]],top),c(0,GSEACompOut$NESorder[[i]],0))
  curves[[i]]<-as.data.frame(curves[[i]])
  curves[[i]]$name<-GSEACompOut$ID[[i]]
  names(curves[[i]])<-c("Rank","ES","Name")}
  curves<-bind_rows(curves)
  if(returnplot=="ES"){p<-EnrichmentScorePlot(curves)
  return(p)}
  if(returnplot=="BC"){p<-BarcodePlot(curves)
  return(p)}
  if(returnplot=="BOTH"){p<-ComboPlot(curves, Title = PathwayName, Table=GSEACompOut[GSEACompOut$pathway==PathwayName,c("ID","padj")])
  return(p)}
  return(curves)
}

BarcodePlot<-function(GSEACompPathway,stacked=FALSE){
  if(stacked==TRUE){GSEACompPathway<-split(GSEACompPathway,GSEACompPathway$Name)
  for (i in 1:length(GSEACompPathway)){GSEACompPathway[[i]]$y<-i}
  GSEACompPathway<-bind_rows(GSEACompPathway)
  GSEACompPathway$y<-((-GSEACompPathway$y)+1)
  p<-ggplot(GSEACompPathway,aes(x=Rank, y=y,xend=Rank, yend=y-1,color=Name))+geom_segment()+theme_classic()+theme(axis.title.y=element_blank(),
                                                                                                                  axis.text.y=element_blank(),
                                                                                                                  axis.ticks.y=element_blank())+xlab("Gene Rank") + xlim(-5,max(GSEACompPathway$Rank+100) )
  }else{
    p<-ggplot(GSEACompPathway,aes(x=Rank, y=-1,xend=Rank, yend=1,color=Name))+geom_segment()+theme_classic()
  }
  return(p)}


EnrichmentScorePlot<-function(GSEACompPathway){
  p<-ggplot(GSEACompPathway, aes(x=Rank, y=ES, color=Name))+geom_line(size=1.5)+geom_hline(yintercept = 0)+theme_classic()+theme(axis.title.x=element_blank(),
                                                                                                                                 axis.text.x=element_blank(),
                                                                                                                                 axis.ticks.x=element_blank())+ylab("Enrichment Score")+xlim(-5,max(GSEACompPathway$Rank+100) )
  return(p)}
ComboPlot<-function(GSEACompPathway, Title="Enrichment Plot",Table=NA){
  EP<-EnrichmentScorePlot(GSEACompPathway)
  if(!is.na(Table)){EP<-EP+annotation_custom(tableGrob(Table, theme = ttheme_minimal(), rows = NULL, cols = c("Name","q Value")), xmin = 10000, ymin = 0.60)}
  BC<-BarcodePlot(GSEACompPathway, stacked=TRUE)
  p<-ggarrange(EP,BC, common.legend = TRUE, ncol = 1, heights = c(0.6,0.25), align = "v", labels = Title)
  
  return(p)}


UniquePaths<-function(GSEAres, IDgrep, degrees_of_freedom="max"){
  todo<-names(table(grep(IDgrep, GSEAres$ID, value=TRUE)))
  res<-c()
  for(i in 1:length(todo)){
    res<-c(res,GSEAres$pathway[GSEAres$ID==todo[[i]]&GSEAres$padj<0.05&GSEAres$NES>0])
  }
  res<-table(res)
  if(degrees_of_freedom=="max"){return(names(res)[res==i])}else{return(names(res)[res>=(i-degrees_of_freedom)])}
}


Meanenrichmentplot<-function(GSEAmain,greplist=c(),pathwaytodo=""){
    final<-list()
    IDssmooth<-paste(greplist, "smooth")
    IDsgray<-paste(greplist, "rough")
    #for each group, go through and make a smooth track 
    for(i in 1:length(greplist)){
    GSEAres<-GSEAmain[grepl(greplist[[i]], GSEAmain$ID) & grepl(pathwaytodo, GSEAmain$pathway),]
    GSEAres$split<-IDsgray[[i]]
    maximum<-max(unlist(GSEAres$rnkorder))
    #generate bins of size 5 (helps with smoothing)
    bins<-seq(1,maximum ,5)
    bins<-c(bins, maximum,maximum)
    res<-list()
    #average anypoint within those bins
    for (j in 1:(length(bins)-1)){
      res[[j]]<-mean(unlist(GSEAres$NESorder)[unlist(GSEAres$rnkorder)>bins[[j]]& unlist(GSEAres$rnkorder)<=bins[[j+1]]])
    }
    res<-c(res, 0)
    #remove any bin that does not have a point 
    bins<-bins[!is.na(res)]
    res<-res[!is.na(res)]
    #make a new track, add it to the overall 
    res<-list(GSEAres$pathway[[1]], NA, NA, NA, NA, NA, NA, NA, list(bins), list(res), "smooth", IDssmooth[[i]])
    GSEAres<-rbind(GSEAres, res)
    GSEAres$NESorder<-lapply(GSEAres$NESorder, as.double)
    curves<-list()
    #generate coordinate curves for each track including the smooth 
    for (j in 1:length(rownames(GSEAres))){curves[[j]]<-list(c(0,GSEAres$rnkorder[[j]],maximum),c(0,GSEAres$NESorder[[j]],0))
    curves[[j]]<-as.data.frame(curves[[j]])
    curves[[j]]$name<-GSEAres$ID[[j]]
    curves[[j]]$split<-GSEAres$split[[j]]
    names(curves[[j]])<-c("Rank","ES","Name","split")}
    final[[i]]<-bind_rows(curves)
    }
    #merge all the groups together
  final<-bind_rows(final)
  #plot things
  final<-ggplot(final[!grepl("smooth", final$split),],aes(x=Rank, y=ES, group=Name,color=split))+geom_line(size=0.25)+geom_hline(yintercept = 0) + 
    geom_smooth(data = final[grepl("smooth", final$split),], aes(x=Rank, y=ES, group=split, color=split), size=1.5, method = "gam")+theme_classic()+
    ggtitle(pathwaytodo)+ylab("Enrichment Score")+xlim(-5,max(maximum+100))
  
  return(final)}
#read in files, make a big table for plotting
results<-list()
results[[1]]<-readRDS("~/Documents/HoLab/JoinedAnalyses/clonalcorrelationfinalGSEA.rds")
results[[2]]<-readRDS("~/Documents/HoLab/JoinedAnalyses/ClonalcorHDGSEA.rds")
results<-unlist(results, recursive = FALSE)
for (i in names(results)){
  results[[i]]<-GSEATable(results[[i]], gmt = fgsea_sets, name = i)
}

results<-GSEAbig(results)


  p1<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._1_stimulated","pHIVresp_..._2_stimulated", "pCMVresp_..._1_stimulated", "pCMVresp_..._2_stimulated","hdCMVresp_HD"), pathwaytodo = "GOLDRATH_ANTIGEN_RESPONSE")
  p2<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._1_stimulated","pHIVresp_..._2_stimulated", "pCMVresp_..._1_stimulated", "pCMVresp_..._2_stimulated","hdCMVresp_HD"), pathwaytodo = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  p3<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._1_stimulated","pHIVresp_..._2_stimulated", "pCMVresp_..._1_stimulated", "pCMVresp_..._2_stimulated","hdCMVresp_HD"), pathwaytodo = "BOSCO_TH1_CYTOTOXIC_MODULE")
  pdf("ClonalCorGSEA/antigenresponsepathways.pdf", width=48, height = 9)
  
  p1+p2+p3
  dev.off()
  
  p1<-Meanenrichmentplot(GSEAmain = results, greplist=c("0_unstimulated","1_unstimulated", "2_unstimulated"), pathwaytodo = "GOLDRATH_ANTIGEN_RESPONSE")
  p2<-Meanenrichmentplot(GSEAmain = results, greplist=c("0_unstimulated","1_unstimulated", "2_unstimulated"), pathwaytodo = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  p3<-Meanenrichmentplot(GSEAmain = results, greplist=c("0_unstimulated","1_unstimulated", "2_unstimulated"), pathwaytodo = "BOSCO_TH1_CYTOTOXIC_MODULE")
  pdf("ClonalCorGSEA/meanunstimulatedpathways.pdf", width=48, height = 9)
  
  p1+p2+p3
  dev.off()
  
  p1<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._1_stimulated","pCMVresp_..._1_stimulated", "pmemory_..._1_stimulated"), pathwaytodo = "GOLDRATH_ANTIGEN_RESPONSE")
  p2<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._1_stimulated","pCMVresp_..._1_stimulated", "pmemory_..._1_stimulated"), pathwaytodo = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  p3<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._1_stimulated","pCMVresp_..._1_stimulated", "pmemory_..._1_stimulated"), pathwaytodo = "BOSCO_TH1_CYTOTOXIC_MODULE")
  pdf("ClonalCorGSEA/Viremicmeanstimulatedpathways.pdf", width=48, height = 9)
  
  p1+p2+p3
  dev.off()
  
  p1<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._2_stimulated","pCMVresp_..._2_stimulated", "pmemory_..._2_stimulated"), pathwaytodo = "GOLDRATH_ANTIGEN_RESPONSE")
  p2<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._2_stimulated","pCMVresp_..._2_stimulated", "pmemory_..._2_stimulated"), pathwaytodo = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  p3<-Meanenrichmentplot(GSEAmain = results, greplist=c("pHIVresp_..._2_stimulated","pCMVresp_..._2_stimulated", "pmemory_..._2_stimulated"), pathwaytodo = "BOSCO_TH1_CYTOTOXIC_MODULE")
  pdf("ClonalCorGSEA/Suppressedmeanstimulatedpathways.pdf", width=48, height = 9)
  
  p1+p2+p3
  dev.off()

  
