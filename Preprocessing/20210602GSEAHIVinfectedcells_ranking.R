library(Seurat)
library(presto)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tibble)

setwd("~/project/JointAnalyses/")

m_df_H<- msigdbr(species = "Homo sapiens", category = "H")
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C2"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C6"), m_df_H)
m_df_H<- rbind(msigdbr(species = "Homo sapiens", category = "C7"), m_df_H)


fgsea_sets<- m_df_H %>% split(x = .$gene_symbol, f = .$gs_name)
pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}
pubfigsqp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 9)}
pubfiglongp<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 48)}

GSEA<-function(seuratobj, groupcol, groups, groupofinterest, genesets){
  genes<- wilcoxauc(seuratobj, group_by=groupcol, groups_use=groups)
  genes<-genes %>%
    dplyr::filter(group == groupofinterest) %>%
    arrange(desc(logFC)) %>% 
    dplyr::select(feature, logFC)
  rnk<- deframe(genes)
  genes<- deframe(genes)
  genes<- fgsea(pathways=genesets, stats = genes,eps=FALSE)
  genes <- genes %>%
    as_tibble() %>%
    arrange(desc(NES))
  genes<-list(genes, rnk)
  names(genes)<-c("results","rnk")
  return(genes)
}

plottop<-function(res, topn=10){
  p<-ggplot(res %>% filter(padj < 0.008) %>% head(n= topn), aes(reorder(pathway, NES), NES)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") + 
    theme_minimal()
  return(p)
}

####comparisons within HIV infected cell population to others 
#this object is just a merge of the two masters subset to only HIV-infected cells
samples<-readRDS("Objects/20210526AllHIVInfectedCells.rds")
#basically just get all the conditions and then iterate through all comparisons
todo<-names(table(samples$HIVcondtion2))

results<-list()
for(i in 1:length(todo)){
  curres2<-list()
  for(j in 1:length(todo)){if(i==j){next()}
    curres<-GSEA(seuratobj = samples, groupcol = "HIVcondtion2",groups = c(todo[[i]], todo[[j]]),
                 groupofinterest = todo[[i]],genesets = fgsea_sets)
    curres2[[j]]<-curres
  }
  if(i==8){names(curres2)<-todo[1:7]}else{names(curres2)<-todo}
  results[[i]]<-curres2
}

names(results)<-todo

saveRDS(results, "Objects/HIVinfectedcells_vs_eachother_GSEAresults.rds")

#this dataset is just the unstimulated and stimulated objects merged, and extra data removed (SCT, ADT, etc.) due to RAM issues
samples<-readRDS("Objects/20210521BothDatasetsLite.rds")
samples$antigen2[is.na(samples$antigen2)]<-""

#simplifying some columns
samples$antigen2<-case_when(samples$antigen2=="hdCMVresp"|samples$antigen2=="pCMVresp"~"CMV", samples$antigen2=="hdmemory" | samples$antigen2=="pmemory"~"Memory",samples$antigen2=="pHIVresp"~"HIV", TRUE~"resting")
samples$infection<-case_when(samples$infection2=="Viremic"|samples$supprstage=="Viremic" ~"Viremic", samples$infection2=="Suppressed_Delayed"|samples$supprstage=="DelayedSuppr"~"Delayed_Suppression",samples$infection2=="Suppressed_Immediate"|samples$supprstage=="ImmedSuppr"~"Immediate_Suppression", TRUE~"Uninfected")
samples$HIV_infected<-case_when(samples$HIV_infected=="infected" | samples$HIVpositive==TRUE ~ "infected", TRUE~"uninfected")
samples$infection2<-case_when(samples$infection2=="Viremic"|samples$supprstage=="Viremic" ~"Viremic", samples$infection2=="Suppressed_Delayed"|samples$supprstage=="DelayedSuppr"|samples$infection2=="Suppressed_Immediate"|samples$supprstage=="ImmedSuppr"~"Suppressed", TRUE~"Uninfected")
#generating condtions column
samples$HIVcondtion2<-paste(samples$HIV_infected, samples$condition, samples$antigen2, samples$infection2,sep="_")
samples<-subset(samples, infection2!="Uninfected")
todo<-names(table(samples$HIVcondtion2))
todo1<-grep("^infected",todo, value = TRUE)
todo2<-grep("^uninfected",todo, value = TRUE)
#iterating through each infected vs uninfected section 
results<-list()
for(i in 1:length(todo1)){
  results[[i]]<-GSEA(seuratobj = samples, groupcol = "HIVcondtion2",groups = c(todo1[[i]], todo2[[i]]),
               groupofinterest = todo1[[i]],genesets = fgsea_sets)
  
}
names(results)<-paste(todo1,"_vs_", todo2, sep="")
saveRDS(results, "Objects/HIVinfectedcells_vs_uninfectedcells_GSEAresults.rds")

