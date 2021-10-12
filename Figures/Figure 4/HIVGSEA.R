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

m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C2"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C6"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C7"), m_df_H)
fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)

pubfig169p<-function(filename){pdf(paste(filename, "pdf", sep="."), width = 16, height = 16)}

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


results<-readRDS("~/Documents/HoLab/JoinedAnalyses/HIVinfectedcells_vs_eachother_infection_blind_GSEAresults.rds")

for (i in names(results)){
  for(j in names(results[[i]])){if(i==j){next()}
    results[[i]][[j]]<-GSEATable(results[[i]][[j]], gmt = fgsea_sets, name = paste(i,"vs",j, sep = "_"))
  }
  
  
}

results<-unlist(results, recursive = FALSE)
results<-results[!unlist(lapply(results, is.null))]
results<-GSEAbig(results)

GSEAEnrichmentPlotComparison<-function(PathwayName, GSEACompOut, returnplot="none", ylim=NULL){
  GSEACompOut<-GSEACompOut[GSEACompOut$pathway==PathwayName,]
  curves<-list()
  for (i in 1:length(rownames(GSEACompOut))){curves[[i]]<-list(c(0,GSEACompOut$rnkorder[[i]],33539),c(0,GSEACompOut$NESorder[[i]],0))
  curves[[i]]<-as.data.frame(curves[[i]])
  curves[[i]]$name<-GSEACompOut$ID[[i]]
  names(curves[[i]])<-c("Rank","ES","Name")}
  curves<-bind_rows(curves)
  if(returnplot=="ES"){p<-EnrichmentScorePlot(curves, ylim=ylim)
  return(p)}
  if(returnplot=="BC"){p<-BarcodePlot(curves)
  return(p)}
  if(returnplot=="BOTH"){p<-ComboPlot(curves, Title = PathwayName, Table=GSEACompOut[GSEACompOut$pathway==PathwayName,c("ID","padj")], ylim=ylim)
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


EnrichmentScorePlot<-function(GSEACompPathway, ylim=NULL){
  if(!is.na(ylim)){p<-ggplot(GSEACompPathway, aes(x=Rank, y=ES, color=Name))+geom_line(size=1.5)+geom_hline(yintercept = 0)+theme_classic()+theme(axis.title.x=element_blank(),
                                                                                                                                 axis.text.x=element_blank(),
                                                                                                                                 axis.ticks.x=element_blank())+ylab("Enrichment Score")+xlim(-5,max(GSEACompPathway$Rank+100) )+ylim(ylim[1],ylim[2])
  }else{p<-ggplot(GSEACompPathway, aes(x=Rank, y=ES, color=Name))+geom_line(size=1.5)+geom_hline(yintercept = 0)+theme_classic()+theme(axis.title.x=element_blank(),
                                                                                                                                       axis.text.x=element_blank(),
                                                                                                                        axis.ticks.x=element_blank())+ylab("Enrichment Score")+xlim(-5,max(GSEACompPathway$Rank+100) )
  }
  return(p)}

ComboPlot<-function(GSEACompPathway, Title="Enrichment Plot",Table=NA, ylim=NULL){
  EP<-EnrichmentScorePlot(GSEACompPathway, ylim=ylim)
  BC<-BarcodePlot(GSEACompPathway, stacked=TRUE)
  p<-ggarrange(EP,BC ,common.legend = TRUE, ncol = 1, heights = c(0.6,0.25), align = "v", labels = Title)
  p<-p/tableGrob(Table)
  return(p)}
GSEAEnrichmentPlotComparison(results[grep("^infected resting",results$ID),], PathwayName = "GSE36476_CTRL_VS_TSST_ACT_16H_MEMORY_CD4_TCELL_YOUNG_UP",returnplot = "BOTH")


UniquePaths<-function(GSEAres, IDgrep, degrees_of_freedom="max"){
  todo<-names(table(grep(IDgrep, GSEAres$ID, value=TRUE)))
  res<-c()
  for(i in 1:length(todo)){
    res<-c(res,GSEAres$pathway[GSEAres$ID==todo[[i]]&GSEAres$padj<0.05&GSEAres$NES>0])
  }
  res<-table(res)
  if(degrees_of_freedom=="max"){return(names(res)[res==i])}else{return(names(res)[res>=(i-degrees_of_freedom)])}
}

#first comparing infected to uninfected counterparts

results<-readRDS("~/Documents/HoLab/JoinedAnalyses/HIVinfectedcells_vs_uninfectedcells_GSEAresults.rds")

for (i in names(results)){
  
  results[[i]]<-GSEATable(results[[i]], gmt = fgsea_sets, name = i)
}

results2<-GSEAbig(results)

pubfig169p("HIVGSEAfigs/resting_HBV")
GSEAEnrichmentPlotComparison(results2[grep("resting", results2$ID),], PathwayName = "WIELAND_UP_BY_HBV_INFECTION",returnplot = "BOTH")
dev.off()
pubfig169p("HIVGSEAfigs/resting_IFNG")
GSEAEnrichmentPlotComparison(results2[grep("resting", results2$ID),], PathwayName = "HALLMARK_INTERFERON_GAMMA_RESPONSE",returnplot = "BOTH")
dev.off()


#next comparing across infection condtions 

results<-readRDS("~/Documents/HoLab/JoinedAnalyses/HIVinfectedcells_vs_eachother_GSEAresults.rds")

for (i in names(results)){
  for(j in names(results[[i]])){if(i==j){next()}
    results[[i]][[j]]<-GSEATable(results[[i]][[j]], gmt = fgsea_sets, name = paste(i,"vs",j, sep = "_"))
  }
  
  
}

results<-unlist(results, recursive = FALSE)
results<-results[!unlist(lapply(results, is.null))]
results<-GSEAbig(results)

#KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION
pubfig169p("HIVGSEAacross/CytorecptHIVCMV")
GSEAEnrichmentPlotComparison(results[grep("infected_stimulated_CMV_Viremic_vs_infected_stimulated_HIV_Viremic|infected_stimulated_CMV_Suppressed_vs_infected_stimulated_HIV_Suppressed", results$ID),], PathwayName = "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",returnplot = "BOTH")
dev.off()


####HIV RNA+ clones vs their non HIV RNA+ counterparts
results2<-readRDS("~/Documents/HoLab/JoinedAnalyses/HIVinfectedclones_vs_uninfectedclones_GSEAresults2.rds")

for (i in names(results2)){
  
  results2[[i]]<-GSEATable(results2[[i]], gmt = fgsea_sets, name = i)
}

results22<-GSEAbig(results2)

#figure 7 
pubfig169p("Documents/HoLab/JoinedAnalyses/HIVclonalGSEAfigsvsself/restmemact")
GSEAEnrichmentPlotComparison(results22[grep("resting|Memory", results22$ID),], PathwayName = "GSE45739_UNSTIM_VS_ACD3_ACD28_STIM_WT_CD4_TCELL_DN",returnplot = "BOTH")
dev.off()
