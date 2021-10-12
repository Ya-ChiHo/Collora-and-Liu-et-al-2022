library(Seurat)
library(ggplot2)
library(dplyr)

#dependent upon Jackutilities.R
#Figure 3 mainly
setwd("~/project/JointAnalyses/")

samples<-readRDS("Objects/210623masterstimulated.rds")
#need to factor cluster levels and annotation orders 
samples$cluster_annotation<-factor(samples$cluster_annotation, levels =  c("CD4-GZMB Th1","CD4-IL2 Th1","CD4-Lymphotoxin",
                                         "CD4-CD45RA","CD4-HLA-DR","CD4-GZMH","CD4-Treg",
                                         "CD4-IFN","CD4-Memory-1","CD4-Memory-2","CD4-Memory-3",
                                         "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-HSP"))


Idents(samples)<-samples$cluster_annotation
samples$antigen5<-case_when(grepl("non",samples$antigen) ~"Memory", T~as.character(samples$antigen))
#ensure all testing and work is done in RNA
DefaultAssay(samples)<-"RNA"
samples<-DietSeurat(samples, assays = "RNA", dimreducs = "umap")

pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

samples<-NormalizeData(samples)
samples<-ScaleData(samples)
#basic plots for figure 3
pubfig169p("Figures/stimulatedHW/clusterumap")
DimPlot(samples, group.by = "cluster_annotation", label=T)
dev.off()

samples$PT<-factor(samples$PT, levels=c("236","640","829","739","799","910","HD2004","HD2022"))
pubfig169p("Figures/stimulatedHW/PTumap")
DimPlot(samples, group.by="PT")+
  scale_color_manual(values=c(`236` = "#0E5381",`640`= "#188EDC",`829`= "#91CCF3",
                              `739` = "#29711E",`799` = "#46C133",`910` = "#A8E59E",
                              'HD2004' = "#666666",`HD2022` = "#CCCCCC"))
dev.off()


samples$PT<-factor(samples$PT, levels=c("236","640","829","739","799","910","HD2004","HD2022"))
pubfig169p("Figures/stimulatedHW/PTumapsplit")
DimPlot(samples, group.by="PT", split.by="PT", ncol=3)+
  scale_color_manual(values=c(`236` = "#0E5381",`640`= "#188EDC",`829`= "#91CCF3",
                              `739` = "#29711E",`799` = "#46C133",`910` = "#A8E59E",
                              'HD2004' = "#666666",`HD2022` = "#CCCCCC"))
dev.off()


# passing colors for antigen specific cells
samples$antigen2<-factor(samples$antigen2, levels=c("pHIVresp","pCMVresp","hdCMVresp","bystander","pmemory","hdmemory"))
pubfig169p("Figures/stimulatedHW/antigenumap")
DimPlot(samples, group.by="antigen2")+
  scale_color_manual(values=c(`hdCMVresp`= "#6CBAEF",`hdmemory`= "#B8B8B8",
                              `pCMVresp` = "#126AA5",`pHIVresp` = "#DF536B",`pmemory` = "#707070"))
dev.off()

pubfig169p("Figures/stimulatedHW/antigenumapsplit")
DimPlot(samples, group.by="antigen2", split.by="antigen2")+
  scale_color_manual(values=c(`hdCMVresp`= "#6CBAEF",`hdmemory`= "#B8B8B8",
                              `pCMVresp` = "#126AA5",`pHIVresp` = "#DF536B",`pmemory` = "#707070"))
dev.off()

### do heatmap for the top10 genes using  cluster_annotation and the correct order
DefaultAssay(samples)<-"RNA"
Idents(samples)<-"cluster_annotation"

samples$cluster_annotation<-factor(samples$cluster_annotation, levels=c("CD4-GZMB Th1","CD4-IL2 Th1","CD4-Lymphotoxin",
                                                                        "CD4-CD45RA","CD4-HLA-DR","CD4-GZMH","CD4-Treg",
                                                                        "CD4-IFN","CD4-Memory-1","CD4-Memory-2","CD4-Memory-3",
                                                                        "CD4-Memory-4","CD4-Memory-5","CD4-MT","CD4-HSP"))


samples.markers <- FindAllMarkers(samples, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- samples.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pubfiglongp("Figures/stimulatedHW/top10heatmap")
DoHeatmap(samples, features = unique(top10$gene)) + theme(text = element_text(size = 4))
dev.off()


# plotting HIV infected cells
pubfig169p("Figures/stimulatedHW/HIVpos")
DimPlot(samples, cells.highlight=HIV_vec)+NoLegend() 
dev.off()


data<-as.data.frame(samples@reductions$umap@cell.embeddings)
data$HIV<-samples$HIVpositive
data$infection<-samples$stage2
colnames(data)<-c("UMAP_1","UMAP_2", "HIV", "infection")


samples$stage2<-factor(samples$stage2, levels=c("Viremic","Suppressed","Uninfected"))
pubfig169p("Figures/stimulatedHW/HIVinfectioncolor")
ggplot(data, aes(x=UMAP_1, y=UMAP_2,color="Uninfected"))+geom_point(size=0.25)+
  geom_point(data=data[data$HIV=="TRUE",], aes(x=UMAP_1, y=UMAP_2, color=infection),size=1)+
  scale_color_manual(values = c("Viremic"="#DF536B", "Suppressed"="#2297E6","Uninfected"="#C2C2C2"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


samples$supprstage<-factor(samples$supprstage, levels=c("Viremic","ImmedSuppr","DelayedSuppr","Uninfected")) 
samples$antigen5<-factor(samples$antigen5, levels=c("HIVresp","CMVresp","Memory"))

#dotplot of key genes, will be edited to only important clusters, across antigen responses 
genelist<-c("TNFRSF9","CD40LG","TBX21","IFNG","IL2","TNF","LTA","LTB","CCL3","CCL4","GZMB","GNLY","NKG7","CTSW","CCR7","SELL","HLA-DRB1","IL2RA")
samples$keep<-samples$cluster_annotation %in% c("CD4-Lymphotoxin", "CD4-GZMH","CD4-GZMB Th1","CD4-IL2 Th1","CD4-CD45RA","CD4-HLA-DR")
samples2<-subset(samples, keep==TRUE)
pubfigsqp("Figures/stimulatedHW/targetteddotplot")
DotPlot(samples2, features = genelist, split.by = "antigen5", cols = c("#6CBAEF","#B8B8B8","#DF536B"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.4), axis.ticks = element_blank())
dev.off()

#pie charts of cluster composition by antigen 

clusterprops<-group_by(samples@meta.data, cluster_annotation, antigen5)%>%summarise(n=n())%>%mutate(prop=n/sum(n))

pubfiglongp("Figures/stimulatedHW/clusterproportions")
ggplot(clusterprops, aes(x="",y=prop, fill=antigen5))+geom_bar(stat="identity")+
  facet_wrap(facets=. ~ cluster_annotation, ncol = 5) + geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(fill = "Category",
       x = NULL,
       y = NULL,
       title = "Proportion of cluster") + 
  coord_polar("y")
dev.off()

#clone size plot
pubfig169p("Figures/stimulatedHW/clonalsize")
FeaturePlot(samples, "log2bulk_freq")
dev.off()
