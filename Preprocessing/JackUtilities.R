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

#these functions are for generating consistently sized pdf figures
pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

#this function is for writing a list of genes to a file for the java version of GSEA
write.GSEA<-function(seuratoutput, fileprefix){
  write.table(seuratoutput["avg_logFC"], file = paste(fileprefix, "rnk", sep="."), quote=F, row.names = T, col.names = F, sep = "\t")
}
#prints the most correlated and anticorrelated genes with a gene of interest
print.cor<-function(cormatrix, gene, ngene=50){
  for (i in 1:length(cormatrix)){
  print(names(cormatrix)[[i]])
    print(paste("top",ngene,"genes correlated with",gene))
  print(sort(cormatrix[[i]][gene,], decreasing = T)[1:ngene])
  print(paste("bottom",ngene,"genes correlated with",gene))
  print(sort(cormatrix[[i]][gene,], decreasing = F)[1:ngene])
  }
}

#grabs the most correlated genes with gene of interest
extract.cor<-function(cormatrix, gene, ngene=50){
  test<-list()
  for (i in 1:length(cormatrix)){test[[i]]<-names(sort(cormatrix[[i]][gene,],decreasing = T)[1:ngene])}
  return(test)
}

#generates a heatmap based on metadata
mean_meta_heatmap<-function(seuratobj,feature,grouping_x, grouping_y, limits=NA){
  meta<-seuratobj@meta.data
  data<-meta%>%group_by(!!sym(grouping_x), !!sym(grouping_y))%>%summarise(vr=mean(!!sym(feature)))
  if(is.na(limits[[1]])){
    print(ggplot(data, aes(fill=vr, x=!!sym(grouping_x), y=!!sym(grouping_y)))+geom_tile()+scale_fill_gradient(low="gray64", high="gray30")+labs(title=feature)+theme_classic())
  }else{
  print(ggplot(data, aes(fill=vr, x=!!sym(grouping_x), y=!!sym(grouping_y)))+geom_tile()+scale_fill_gradient(low="gray64", high="gray30", limits=limits, oob=squish)+labs(title=feature)+theme_classic())}
}
#used to integrating AUCell data 
Add_AUC_calls<-function(seuratobj, AUC_asssignments, labels=NA){
  if (is.na(labels)){
    labels<-names(AUC_asssignments)
  }
  for (i in 1:length(labels)){
    currentset<-as.data.frame(AUC_asssignments[[labels[i]]]$assignment)
    if(length(rownames(currentset))<1){next()}
    rownames(currentset)<-currentset[,1]
    currentset[[labels[[i]]]]<-labels[[i]]
    names(currentset)[1]<-"Dump"
    seuratobj<-AddMetaData(seuratobj,currentset)
  }
  return(seuratobj)
}
#column normalize for dataframe
{col_normalize <- function(x){
  (x / max(x, na.rm=TRUE))
}

#custom heatmap for generating a by gene normalized dataframe split by a given variable 

genes_heatmap<-function(seuratobj, featurelist, grouping_var, level=TRUE){
  
  #step 0 grab the actual data, both the matrix of features and the metadata, then bind them together and flip so its columns - features, rows observations 
  meta<-seuratobj@meta.data
  genes<-featurelist[featurelist%in%rownames(seuratobj@assays$SCT)]
  meta_features<-c(featurelist[featurelist%in%colnames(meta)],"orig.ident")
  matrix<-seuratobj@assays$SCT@data[genes,]
  matrix<-t(matrix)
  matrix<-cbind(matrix, meta[,meta_features])
  
  #step 1 is to split things by celltype 
  matrix_group<-group_by(matrix,!!sym(grouping_var), as_index=FALSE)
  #step 2 is to calculate the mean expression of a given gene with that cell type
  matrix_group<-summarise_if(matrix_group,is.numeric,"mean")
  #return(matrix_group)
  #step 3 is to scale/row normalize
  matrix_group2<-summarise_if(matrix_group,is.numeric,"mean")%>%as.data.frame()%>%transmute_if(is.numeric, col_normalize)
  matrix_group2$group<-matrix_group[[grouping_var]]
  matrix_group<-matrix_group2
  #step4 is to colapse all the numeric features into a single column with the corresponding grouping var in all the remaining 
  matrix_group<-gather(matrix_group,Gene, Normalized_Expression,  1:length(matrix_group)-1)
  if(level==TRUE){matrix_group$Gene<-factor(matrix_group$Gene, levels = featurelist)}
  #step 5 is to return dataframe for plotting 
  return(matrix_group)
}
}

#GSEA Ploting - many of these functions were refactored/tooled for GSEA in their respective GSEA scripts. 

#function 0 takes all my inputs and figures it out (reduce human error)
GSEAWrap<-function(gmtfile,seuratoutput,name="undefined",metric="avg_logFC"){
  #read in files
  gmt<- gmtPathways(gmtfile)
  rnk<- setNames(seuratoutput[[metric]],as.character(rownames(seuratoutput)))
  #run GSEA
  results <- fgsea(pathways=gmt, stats=rnk,minSize = 15,maxSize = 500, nperm=100000)
  results$name<-name
  #construct object with GSEA results and rnk file
  results<-list(results, rnk, gmt)
  return(results)
}

#Function 1 augments the existing result table with pathway information, including the entire rank order list of gene hits and NES scores 
GSEATable<-function(GSEAwrap_out, gseaParam=1){
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
 #basically here I just want to get everything into a character vector
 for (i in 1:length(finalresults$pathway)){
  pathway <- unname(as.vector(na.omit(match(GSEAwrap_out[3][[1]][[finalresults$pathway[[i]]]], names(statsAdj)))))
  pathway <- sort(pathway)
  NES[[i]]<-calcGseaStat(statsAdj, selectedStats = pathway,returnAllExtremes = TRUE)
  NES[[i]]<-NES[[i]]$tops
  pathwaysetup[[i]]<-pathway}
 finalresults$rnkorder<-pathwaysetup
 finalresults$NESorder<-NES
 return(finalresults)
}

#Function 2 uses function 1 results to wrangle several of the same pathway for different conditions 
GSEAComparisons<-function(gmtfile,seuratobj_lists, IDs, metric="avg_logFC"){
  results<-list()
  for (i in 1:length(seuratobj_lists)){
    results[[i]]<-GSEAWrap(gmtfile, seuratobj_lists[[i]], IDs[[i]], metric=metric)
    results[[i]]<-GSEATable(results[[i]])
  }
  results<-bind_rows(results)
  return(results)
}
#does the actual plotting, can plot multiple conditions at once
GSEAEnrichmentPlotComparison<-function(PathwayName, GSEACompOut, returnplot="none"){
  GSEACompOut<-GSEACompOut[GSEACompOut$pathway==PathwayName,]
  curves<-list()
  for (i in 1:length(rownames(GSEACompOut))){curves[[i]]<-list(c(0,GSEACompOut$rnkorder[[i]],17062),c(0,GSEACompOut$NESorder[[i]],0))
  curves[[i]]<-as.data.frame(curves[[i]])
  curves[[i]]$name<-GSEACompOut$name[[i]]
  names(curves[[i]])<-c("Rank","ES","Name")}
  curves<-bind_rows(curves)
  if(returnplot=="ES"){p<-EnrichmentScorePlot(curves)
    return(p)}
  if(returnplot=="BC"){p<-BarcodePlot(curves)
  return(p)}
  if(returnplot=="BOTH"){p<-ComboPlot(curves, Title = PathwayName, Table=GSEACompOut[GSEACompOut$pathway==PathwayName,c("name","padj")])
  return(p)}
  return(curves)
}



#various minor ploting functions

#barcodeplot version for if I can get two separate graphs working right
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


#modified star command and plot command from: https://github.com/andrewheiss/Attitudes-in-the-Arab-World/blob/master/figure12.R



# Build dataframe for plotting

#modified from: https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/ until line 249

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90, hjust =1, vjust = 0.2), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


mode <- function(cat){
  names(table(cat))[which.max(table(cat))]
}

#takes in a list of cor matrix and genes and saves one plot per cor matrix 
corheatmap<-function(cor, genes, prefix, n_ext=50,exp=TRUE){
  samples<-names(cor)
  plot<-c()
  if(exp==TRUE){for(i in 1:length(genes)){
    plot<-unique(c(plot,extract.cor(cor,genes[[i]], ngene = n_ext)))
  }
  plot2<-c()
  for(i in 1:length(plot)){
    plot2<-c(plot2, plot[[i]])
  }
  plot2<-unique(plot2)}else{plot2<-genes}
  
  for (z in 1:length(cor)){
    corx<-cor[[z]][plot2,plot2]
    x<-setNames(melt(corx), c("x","y","weight"))
    x$x<-factor(x$x, levels = plot2)
    x$y<-factor(x$y, levels = plot2)
    x[is.na(x)]<-0
    p<-ggplot(x, aes(x=x,y, fill=weight))+geom_tile()+
      scale_fill_gradient2(low = "blue",high = "red", limits=c(-0.3, 0.3), oob=squish)+
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.4))
    savename<-paste(prefix,"_",samples[z],".pdf", sep = "")
    ggsave(filename = savename, plot = p, device = "pdf", width = 10, height = 9)
  }
}
#makes a combo of the two inputmatrix
splitcorheatmap<-function(cor, genes, prefix, n_ext=50,exp=TRUE){
  samples<-names(cor)
  plot<-c()
  if(exp==TRUE){for(i in 1:length(genes)){
    plot<-unique(c(plot,extract.cor(cor,genes[[i]], ngene = n_ext)))
  }
    plot2<-c()
    for(i in 1:length(plot)){
      plot2<-c(plot2, plot[[i]])
    }
    plot2<-unique(plot2)}else{plot2<-genes}
    
    corx<-cor[[1]][plot2,plot2]
    cory<-cor[[2]][plot2,plot2]
    corx[lower.tri(corx)]<-cory[lower.tri(cory)]
    x<-setNames(melt(corx), c("x","y","weight"))
    x[is.na(x)]<-0
    x$x<-factor(x$x, levels = plot2)
    x$y<-factor(x$y, levels = plot2)
    p<-ggplot(x, aes(x=x,y, fill=weight))+geom_tile()+
      scale_fill_gradient2(low = "blue",high = "red", limits=c(-0.3, 0.3), oob=squish)+
      theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.4))+theme(axis.ticks = element_blank())
    savename<-paste(prefix,"_",samples[1],"_",samples[2],"split.pdf", sep = "")
    ggsave(filename = savename, plot = p, device = "pdf", width = 10, height = 9)
  
}
#makes all hte comparisons in each correlation plot
splitwrap<-function(cor, genes, n_ext=50, exp=FALSE, prefix="matrix"){
  genes<-unique(genes)
  plot<-c()
  if(exp==TRUE){for(i in 1:length(genes)){
    plot<-unique(c(plot,extract.cor(cor,genes[[i]], ngene = n_ext)))
  }
    plot<-unique(unlist(plot))}else(plot<-genes)
  for (i in 1:length(cor)){
    for (j in 1:length(cor)){
      splitcorheatmap(cor[c(i,j)], genes = plot, prefix=prefix, exp=FALSE)
    }
  }
}
#basic volcano function for seurat output
VolPlot<-function(results,top=TRUE, adthresh=TRUE, thresh=0.2, threshn=30, Title=NA){
  results$value<-(-log10(results$p_val_adj))
  results$value[results$value==Inf]<-300
  results$gene<-rownames(results)
  p<-ggplot(results,aes(x=avg_logFC, y=value, label=gene))+geom_point()+theme_classic()
  if(adthresh==TRUE){thresh<-AdaptiveThreshold(results, Tstart=thresh,n=threshn)}
  if(top==TRUE){
    p<- p + geom_text_repel(aes(label=ifelse((avg_logFC>thresh|avg_logFC<(-thresh))&value>(-log10(0.01)) ,as.character(gene),'')),hjust=1,vjust=1)
  }
  if(!is.na(Title)){
    p<-p+ggtitle(Title)+xlab("Average Log Fold Change")+ylab("-log10 adjusted p value")+theme(plot.title = element_text(hjust = 0.5))
  }
  return(p)
}
#quick selection of top n genes basically
AdaptiveThreshold<-function(results, n=30,Tstart=0.25){
  m=Inf
  while(n<=m){
    Tstart=Tstart+0.05
    m<-length(rownames(results)[(results$avg_logFC>Tstart|results$avg_logFC<(-Tstart))&results$value>10])
  }
  return(Tstart)
}

#as above but better for certain circumstances. 
splitwrap<-function(cor, genes, n_ext=50, exp=FALSE, prefix="matrix"){
  genes<-unique(genes)
  plot<-c()
  if(exp==TRUE){for(i in 1:length(genes)){
    plot<-unique(c(plot,extract.cor(cor,genes[[i]], ngene = n_ext)))
  }
    plot<-unique(unlist(plot))}else(plot<-genes)
  for (i in 1:length(cor)){
    for (j in 1:length(cor)){
      splitcorheatmap(cor[c(i,j)], genes = plot, prefix=prefix, exp=FALSE)
    }
  }
}
