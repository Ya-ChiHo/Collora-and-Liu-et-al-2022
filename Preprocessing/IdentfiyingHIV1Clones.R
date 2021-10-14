library(Seurat)
library(dplyr)

setwd("project/JointAnalyses/")

Jack<-readRDS("Objects/FASTMNN.rds")
Runxia<-readRDS("Objects/210506phateCITEbina.rds")

Runxia[["HTO"]]<-NULL
Runxia[["SCT"]]<-NULL
Runxia[["phate0"]]<-NULL
Runxia[["phatea"]]<-NULL
Runxia<-DietSeurat(Runxia, dimreducs="umap")

saveRDS(Runxia, "Objects/210512RunxiaLite.rds")

Jack[["integrated"]]<-NULL
DefaultAssay(Jack)<-"RNA"
Jack[["SCT"]]<-NULL
Jack<-DietSeurat(Jack, dimreducs="umap")
Jack<-NormalizeData(Jack)

saveRDS(Jack, "Objects/210512JackLite.rds")

clonalmeta<-readRDS("Metadata/clonalmetadata.rds")
Jack<-AddMetaData(Jack,clonalmeta)

JackMeta<-Jack@meta.data
RunxiaMeta<-Runxia@meta.data

#make a master list of patients and TCRBs that are captured as HIV positive
todo<-JackMeta[JackMeta$HIV_infected=="infected",c("junction","PT")]
todo<-todo[!is.na(todo$junction),]
todo2<-RunxiaMeta[RunxiaMeta$HIVpositive==TRUE, c("TRB_junction", "PT")]
todo2<-todo2[!todo2$TRB_junction=="",]
colnames(todo2)<-c("junction","PT")
todo<-rbind(todo,todo2)
#for each of my TCRs, lets find what cells in Runxia's data are present

results<-c()
for(i in 1:nrow(RunxiaMeta)){
  for(j in 1:nrow(todo))
    if(RunxiaMeta$PT[[i]]==todo$PT[[j]] & RunxiaMeta$TRB_junction[[i]]==todo$junction[[j]]){
      results<-c(results, rownames(RunxiaMeta)[[i]])
      break
    }
}
JackMeta<-JackMeta[!is.na(JackMeta$orig.ident),]
Jackcells<-c()
#now check mine
for(i in 1:nrow(JackMeta)){
  for(j in 1:nrow(todo))
    if(JackMeta$PT[[i]]==todo$PT[[j]] & JackMeta$junction[[i]]==todo$junction[[j]]){
      Jackcells<-c(Jackcells, rownames(JackMeta)[[i]])
      break
    }
}



Jackcells<-JackMeta[Jackcells,]
#no overlap PT needed at this stage since each junction is mututally exclusive
test<-split(Jackcells, Jackcells$PT)
for (i in 1:length(test)){
  print(names(test)[[i]])
  for (j in 1:length(test)){
    print(names(test)[[j]])
    print(table(test[[i]]$junction %in%test[[j]]$junction))
  }}
  
  

#40 clones containing HIV-1 infected cells in just Jack's datasets 
table(table(Jackcells$junction)>1)

#751 cells, up from 537
sum(table(Jackcells$junction[Jackcells$junction %in% names(table(Jackcells$junction))[table(Jackcells$junction)>1]]))

Jackcellnames<-rownames(Jackcells[Jackcells$junction %in% names(table(Jackcells$junction))[table(Jackcells$junction)>1],])

#save marked cell names 
saveRDS(Jackcellnames, "Metadata/210512JackdataHIVClones.rds")

Runxiacells<-RunxiaMeta[results,]

#no overlap PT needed at this stage since each junction is mututally exclusive
test<-split(Runxiacells, Runxiacells$PT)

for (i in 1:length(test)){
  print(names(test)[[i]])
  for (j in 1:length(test)){
    print(names(test)[[j]])
    print(table(test[[i]]$TRB_junction %in%test[[j]]$TRB_junction))
  }}
  
  
Runxiacellnames<-rownames(Runxiacells[Runxiacells$TRB_junction %in% names(table(Runxiacells$TRB_junction))[table(Runxiacells$TRB_junction)>1],])

saveRDS(Runxiacellnames, "Metadata/210512RunxiadataHIVClones.rds")


#now with the merged data 
Jackcellsubset<-Jackcells[,"junction"]
Runxiacellsubset<-Runxiacells[,"TRB_junction"]

masterclonelist<-c(rownames(Jackcells)[Jackcells$junction%in%bothclones],rownames(Runxiacells)[Runxiacells$TRB_junction%in%bothclones])
saveRDS(masterclonelist, "Metadata/2120512AllHIVclones.rds")

#which subset of clones are only apparent with the cells that I have 

todo<-JackMeta[JackMeta$HIV_infected=="infected",]
todo<-todo[!is.na(todo$junction),]
#24 clones identified by bulk in my data
junctions<-todo$junction[todo$results_count>1]

#514 cells in those clones 
junctions<-JackMeta[JackMeta$junction%in%junctions, ]

#saving those 514 as the training for scRFE round 2 
write.table(rownames(junctions),"sklearn/ClonesHIVJackData.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

#saving the normalized data and metadata for scRFE round 2 
Jack$keep<-!is.na(Jack$junction) & Jack$infection !="Uninfected"
Jack<-subset(Jack, keep ==TRUE)
#78736 cells remain
Jack$HIV<- colnames(Jack)%in%rownames(junctions)
library(Matrix)
writeMM(Jack@assays$RNA@data,"sklearn/NormalizedJackData.mtx")
write.table(colnames(Jack), "sklearn/CellnamesJack.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(Jack), "sklearn/GenenamesJack.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

#also writing the clones from Runxia's data as a "held out" set 

todo<-RunxiaMeta[RunxiaMeta$HIVpositive==TRUE,]
todo<-todo[!todo$TRB_junction=="",]
oldjuncts<-junctions
junctions<-todo$TRB_junction[todo$bulk_n>1]

junctions<-JackMeta[JackMeta$junction%in%junctions, ]

#also no shared junctions across individuals
table(junctions$PT, junctions$junction)

#just need to remove one 
remove<-intersect(unique(todo$TRB_junction), unique(oldjuncts$junction))
#248 cells remain in these other clones (as determined by bulk level)
junctions<-junctions[junctions$junction!=remove,]

#saving those cell nmaes
write.table(rownames(junctions),"sklearn/ClonesHIVRunxiaData.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
