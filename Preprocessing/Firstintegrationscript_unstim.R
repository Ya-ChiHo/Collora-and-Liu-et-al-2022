library(dplyr)
library(Seurat)
library(future)
library(SeuratWrappers)
plan(multicore, workers=10)
options(future.globals.maxSize = 50000 * 1024^2)

pre_filtering_tcr<- function(tcrfile){
  tcr <- read.csv(tcrfile)
  tcr<-tcr[tcr$productive=="True",]
  return(tcr)
}
#this function will then find barcodes that occur more often than a threshold
tcr_mark_multi<-function(tcrobj, threshold=3){
  tcr_multi<-sort(table(tcrobj$barcode), decreasing = T)
  tcr_multi<-tcr_multi[tcr_multi>threshold]
  #remove those with 3 of the same A/B too
  tcrA <- tcrobj[tcrobj$chain=="TRA",c(1,7:10,13)]
  tcrB <- tcrobj[tcrobj$chain=="TRB",c(1,7:10,13)]
  tcr_multiA<-sort(table(tcrA$barcode), decreasing = T)
  tcr_multiB<-sort(table(tcrB$barcode), decreasing = T)
  #assemble a list
  tcr_multi<-union(names(tcr_multi),names(tcr_multiB[tcr_multiB>(threshold-1)]))
  tcr_multi<-union(tcr_multi,names(tcr_multiA[tcr_multiA>(threshold-1)]))
  return(tcr_multi)
}
#removes the multiplit cells from the seurat object
tcr_doublet_removal<- function(tcrmultinames, seuratobj){
  tcrmultinames <- gsub("-1", "", tcrmultinames)
  seuratobj<-seuratobj[,tcrmultinames, invert=T]
  return(seuratobj)
}
#removes the multiplit cells from the TCR object
tcr_removal<- function(tcrobj, tcrmultinames){
  tcrobj<-tcrobj[!grepl(paste(tcrmultinames,collapse="|"), tcrobj$barcode),]
  return(tcrobj)
}
#adds the alpha and beta chains to the suerat object, must be filtered down to just 2 at a max for prioritizing 
add_tcr_info<-function(tcrobj, seuratobj){
  tcrobj$barcode<- gsub("-1", "", tcrobj$barcode)
  #split by TCR A and TCR B, keep only columns with info we care about, find out which of each is clonal in this sample ie present in greater than 2 barcodes.  
  tcrA <- tcrobj[tcrobj$chain=="TRA",c(1,7:10,13)]
  #first identify and pull out those that have two for prioritization 
  tcrA_multi<- sort(table(tcrA$barcode), decreasing = T)
  tcrA_multi <-names(tcrA_multi[tcrA_multi>1])
  tcrA_to_prioritize<- tcrA[tcrA$barcode%in%tcrA_multi,]
  tcrA_priority_free <-tcrA[!tcrA$barcode%in%tcrA_multi,]
  #second identify what receptors are clonal in this sample, and prioritize the remaining based on clonality. 
  tcrA_clonal<- sort(table(tcrA$cdr3), decreasing = T)
  tcrA_clonal <- tcrA_clonal[tcrA_clonal>1]
  #get obs per cdr3/tcr  
  tcrA_to_prioritize$clonalityscore <- tcrA_clonal[tcrA_to_prioritize$cdr3]
  
  #reindex
  rownames(tcrA_to_prioritize)<- 1:nrow(tcrA_to_prioritize)
  #removeNA
  tcrA_to_prioritize$clonalityscore[is.na(tcrA_to_prioritize$clonalityscore)]<-0
  
  #find the one that is more clonal, defaulting to the one with more reads if equal
  priorityvector<-c()
  for(i in seq(1,length(tcrA_to_prioritize$clonalityscore)-1,2)){
    if(tcrA_to_prioritize$clonalityscore[i]<tcrA_to_prioritize$clonalityscore[i+1]){priorityvector<-c(priorityvector,2,1)
    }
    else{priorityvector<-c(priorityvector,1,2)}
  }
  tcrA_to_prioritize$priority<-priorityvector
  #split priority vectors
  tcrA_p1<- tcrA_to_prioritize[tcrA_to_prioritize$priority==1,1:(length(colnames(tcrA_to_prioritize))-2)]
  tcrA_p2<- tcrA_to_prioritize[tcrA_to_prioritize$priority==2,1:(length(colnames(tcrA_to_prioritize))-2)]
  #add the priority 1s togeter
  tcrA<-rbind(tcrA_priority_free, tcrA_p1)
  rownames(tcrA)<-tcrA$barcode
  rownames(tcrA_p2)<-tcrA_p2$barcode
  #add in priorty 2
  tcrA<-merge(tcrA, tcrA_p2, by=0, all=T)
  #cleanup the dataframe
  Anames<- colnames(tcrA)
  Anames<- gsub("y", Anames, replacement = "Priority_2")
  Anames<- gsub("x", Anames, replacement = "Priority_1")
  colnames(tcrA)<- Anames
  rownames(tcrA)<-tcrA$Row.names
  tcrA<- tcrA[,c(3:7,9:13)]
  ####do same as above but for B
  tcrB <- tcrobj[tcrobj$chain=="TRB",c(1,7:10,13)]
  
  tcrB_multi<- sort(table(tcrB$barcode), decreasing = T)
  tcrB_multi <-names(tcrB_multi[tcrB_multi>1])
  tcrB_to_prioritize<- tcrB[tcrB$barcode%in%tcrB_multi,]
  tcrB_priority_free <-tcrB[!tcrB$barcode%in%tcrB_multi,]
  #second identify what receptors are clonal in this sample, and prioritize the remaining based on clonality. 
  tcrB_clonal<- sort(table(tcrB$cdr3), decreasing = T)
  tcrB_clonal <- tcrB_clonal[tcrB_clonal>1]
  #get obs per cdr3/tcr  
  tcrB_to_prioritize$clonalityscore <- tcrB_clonal[tcrB_to_prioritize$cdr3]
  
  #reindex
  rownames(tcrB_to_prioritize)<- 1:nrow(tcrB_to_prioritize)
  #removeNA
  tcrB_to_prioritize$clonalityscore[is.na(tcrB_to_prioritize$clonalityscore)]<-0
  
  #find the one that is more clonal, defaulting to the one with more reads if equal
  priorityvector<-c()
  for(i in seq(1,length(tcrB_to_prioritize$clonalityscore)-1,2)){
    if(tcrB_to_prioritize$clonalityscore[i]<tcrB_to_prioritize$clonalityscore[i+1]){priorityvector<-c(priorityvector,2,1)
    }
    else{priorityvector<-c(priorityvector,1,2)}
  }
  tcrB_to_prioritize$priority<-priorityvector
  #split priority vectors
  tcrB_p1<- tcrB_to_prioritize[tcrB_to_prioritize$priority==1,1:(length(colnames(tcrB_to_prioritize))-2)]
  tcrB_p2<- tcrB_to_prioritize[tcrB_to_prioritize$priority==2,1:(length(colnames(tcrB_to_prioritize))-2)]
  #add the priority 1s togeter
  tcrB<-rbind(tcrB_priority_free, tcrB_p1)
  rownames(tcrB)<-tcrB$barcode
  rownames(tcrB_p2)<-tcrB_p2$barcode
  #add in priorty 2
  tcrB<-merge(tcrB, tcrB_p2, by=0, all=T)
  #cleanup the dataframe
  Bnames<- colnames(tcrB)
  Bnames<- gsub("y", Bnames, replacement = "Priority_2")
  Bnames<- gsub("x", Bnames, replacement = "Priority_1")
  colnames(tcrB)<- Bnames
  rownames(tcrB)<-tcrB$Row.names
  tcrB<- tcrB[,c(3:7,9:13)]
  #prep both to merge 
  colnames(tcrB)<- paste("B",colnames(tcrB),sep="_")
  colnames(tcrA)<- paste("A",colnames(tcrA),sep="_")
  #merge and cleanup
  tcrobj<-merge(tcrB, tcrA, by=0, all=T)
  rownames(tcrobj)<-tcrobj$Row.names
  tcrobj<-tcrobj[,2:length(colnames(tcrobj))]
  #add it to seuratobj and return
  seuratobj<-AddMetaData(object=seuratobj, metadata=tcrobj)
  return(seuratobj)
}
#calls all the above in an automated way
automated_TCR<-function(tcrfile, seuratobj){
  tcrobj<-pre_filtering_tcr(tcrfile)
  tcrmultinames<-tcr_mark_multi(tcrobj)
  seuratobj<-tcr_doublet_removal(tcrmultinames, seuratobj)
  tcrobj<-tcr_removal(tcrobj, tcrmultinames)
  seuratobj<-add_tcr_info(tcrobj, seuratobj)
  return(seuratobj)
}

#will change things to binary, thresholding positive cells as those above  0.999 percentile 
aggregate_TCR<-function(matrixobj){
  #reflist to remove
  TCR_gene_list<-c("TRGV","TRDV","TRAV", "TRBV")
  #first get the matrix going
  matrix<- matrixobj$`Gene Expression`
  sumofgene=list()
  #get the sum of each of those
  for(i in 1:length(TCR_gene_list)){
    gene<-grep(TCR_gene_list[i],rownames(matrixobj$`Gene Expression`))
    sumofgene[[i]]<-as.data.frame(Matrix::colSums(matrix[gene,]))
    colnames(sumofgene[[i]])<-TCR_gene_list[i]
    rownames(sumofgene[[i]])<-colnames(matrix)
  }
  #restack it into a matrix
  two<-bind_cols(sumofgene)
  rownames(two)<-colnames(matrix)
  two<-t(two)
  #remove the old values
  matrix<-matrix[!grepl(paste(TCR_gene_list,collapse="|"), rownames(matrix)),]
  matrix<-rbind(matrix,two)
  matrixobj$`Gene Expression`<-matrix
  #giveitback
  return(matrixobj)
}

add_HIV<-function(Seuratobj, HIVfile){
  #HIV file generated by first mapping via STAR, extracting corresponding BC/UMI read and running this command:
  #awk '{if(NR%4==2){print substr($1,1,26)}}' positivreads.fq.gz |sort >positivebc.txt
  #first read an clean the HIV data
  HIV2<-read.table(HIVfile)
  #pulloff readcounts for each cell barcode for second list
  HIV3<-table(HIV2$V1)
  #remove duplicate reads
  HIV2<-unique(HIV2)
  #pull barcodes only
  HIV2<-separate(HIV2, V1, into=c("BC","UMI"), sep=16)
  #count barcodes (which are really UMI now)
  HIV2<-table(HIV2$BC)
  #removing all those that only have 1 umi, likely not real, but retaining counts to also add
  HIVUMI<-HIV2
  HIV2<-HIV2[HIV2>1]
  HIV2<-names(HIV2)
  #select only those form HIV3 which are captured 5 or more times
  HIV3<-HIV3[HIV3>4]
  #split to only cell BC
  for (i in 1:length(HIV3)){HIV3[[i]]<-substr(names(HIV3)[[i]], 1,16)}
  HIV3<-unique(HIV3)
  HIV<-as.data.frame(HIVUMI)
  rownames(HIV)<-HIV$Var1
  #now need to figure out intersections
  method1=c()
  method2=c()
  both=c()
  for (i in 1:length(rownames(HIV))){
    if (rownames(HIV)[[i]]%in%HIV2&rownames(HIV)[[i]]%in%HIV3){
      method1[[i]]="yes"
      method2[[i]]="yes"
      both[[i]]="yes"
    }
    else if(rownames(HIV)[[i]]%in%HIV2){
      method1[[i]]="yes"
      method2[[i]]="no"
      both[[i]]="no"
    }
    else if(rownames(HIV)[[i]]%in%HIV3){
      method1[[i]]="no"
      method2[[i]]="yes"
      both[[i]]="no"
    } 
    else{
      method1[[i]]="no"
      method2[[i]]="no"
      both[[i]]="no"
    }
    }
HIV$UMImethodpositive<-method1
HIV$readmethodpositive<-method2
HIV$bothmethodpositive<-both
HIV<-HIV[,2:5]
colnames(HIV)<-c("HIVUMI",colnames(HIV)[2:4])
Seuratobj<-AddMetaData(Seuratobj, HIV)
  #giveitback
  return(Seuratobj)
}

matrixfile<-list.files("scratch60/RNA/unfilteredmatrix")
TCRanno<- list.files("scratch60/RNA/TCRdata/annotations")
identifiers<-gsub("matrix","", matrixfile)
identifiers2<-gsub("annotations.csv", "",TCRanno)
metadata<-read.table("scratch60/RNA/metadata.csv",sep=",", header=T) 
identifiers==identifiers2
rawdata=list()
processeddata=list()
for(i in 1:length(matrixfile)){
  #this section loads the data and puts it into a seurat object
  rawdata[[i]] <- Read10X(data.dir = paste("~/scratch60/RNA/unfilteredmatrix",matrixfile[[i]], sep = "/"))
  rawdata[[i]] <- aggregate_TCR(rawdata[[i]])
  processeddata[[i]] <- CreateSeuratObject(counts = rawdata[[i]]$`Gene Expression`, project = identifiers[[i]])
  processeddata[[i]]$delay<-metadata[[identifiers[[i]]]][[1]]
  processeddata[[i]]$infection<-metadata[[identifiers[[i]]]][[2]]
  processeddata[[i]]$HIV<-metadata[[identifiers[[i]]]][[3]]
  processeddata[[i]]$batch<-metadata[[identifiers[[i]]]][[4]]
  processeddata[[i]][["ADT"]] <- CreateAssayObject(counts = rawdata[[i]]$`Antibody Capture`)
  processeddata[[i]]<-automated_TCR(paste("~/scratch60/RNA/TCRdata/annotations",TCRanno[[i]], sep = "/"),processeddata[[i]])
  
  
  #perform the processing on each patient, this is variable based on which patient population is used 
  processeddata[[i]] <- PercentageFeatureSet(processeddata[[i]], pattern = "^MT-", col.name = "percent.mt")
  processeddata[[i]]<- subset(processeddata[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
  #need this so we can ensure taht we can add on HIV later 
  processeddata[[i]]<-RenameCells(object = processeddata[[i]], add.cell.id = identifiers[[i]])
  #sctransform to start, binarize ADT library and then normalize it. Stopped regressing out percent.mt, it probably doesnt make too much of a difference since its already filtered 
  processeddata[[i]]<-NormalizeData(processeddata[[i]], assay="ADT", normalization.method = "CLR")
  
}

saveRDS(processeddata, "scratch60/RNA/checkpoint/checkpoint1.rds")
#perform SCT specific methods
for(i in 1:length(processeddata)){processeddata[[i]]<-SCTransform(processeddata[[i]])}
features<-SelectIntegrationFeatures(processeddata)
processeddata<-PrepSCTIntegration(processeddata, anchor.features = features)
processeddata <- RunFastMNN(processeddata, normalization.method = "SCT")
saveRDS(processeddata, "scratch60/RNA/checkpoint/checkpoint2.rds")
processeddata<- FindNeighbors(processeddata, reduction="mnn", dims=1:30)
processeddata<-FindClusters(processeddata, resolution = 2)
processeddata <- RunUMAP(processeddata, reduction="mnn", dims=1:30)
saveRDS(processeddata,"scratch60/RNA/checkpoint3.rds")
pdf("scratch60/RNA2/plots/UMAP1h.pdf")
DimPlot(processeddata, label = T) + NoLegend()
dev.off()
pdf("scratch60/RNA2/plots/UMAPbyPTh.pdf")
DimPlot(processeddata,group.by = "orig.ident")
dev.off()
