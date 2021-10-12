# need to set working directory?

library(dplyr)
library(Seurat)
library(future)
# -p bigmem -c 10 -N 1 --mem 750g -t 3- 
plan(multicore, workers=10)
options(future.globals.maxSize = 500000 * 1024^2)

### #put in new version of TCR functions 

library(dplyr)
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


matrixfile<-list.files("~/scratch60/2021analysis/raw_feature_bc_matrix")
TCRanno<- list.files("~/scratch60/2021analysis/TCRannotations")
# if the name of matrix files is hdm1matrix for example, use
# identifiers<-gsub("matrix","", matrixfile)
# In this case, the name of matrix file is hdm1 for example, use
# identifiers<-matrixfile
identifiers<-matrixfile

# if the annotation name is hdm1annotations.csv for example, use
# identifiers2<-gsub("annotations.csv","",TCRanno)
# In this case, the annotation name is hdm1_all_contig_annotations for example, use
# identifiers2<-gsub("_all_contig_annotations.csv", "",TCRanno)
identifiers2<-gsub("_all_contig_annotations.csv", "",TCRanno)
#can remove, or just type it in 
metadata<-read.table("~/scratch60/2021analysis/metadata.txt",sep="\t", header=T) 

# now just grab those rows with our hashes 
# seruatobject[["HTO"]]<- seuratobject@assays$ADT@counts[currenthashes,]

#make sure in same order 
identifiers==identifiers2

# The rawdata list and processeddata list are meant to be empty at that stage. You're initializing them prior to using them in the loop. After the loop they will have your libraries in them. 
rawdata=list()
processeddata=list()

for(i in 1:length(matrixfile)){
  #this section loads the data and puts it into a seurat object
  rawdata[[i]] <- Read10X(data.dir = paste("~/scratch60/2021analysis/raw_feature_bc_matrix",matrixfile[[i]], sep = "/"))
  #skip rawdata[[i]] <- aggregate_TCR(rawdata[[i]])
  processeddata[[i]] <- CreateSeuratObject(counts = rawdata[[i]]$`Gene Expression`, project = identifiers[[i]])
  processeddata[[i]][["ADT"]] <- CreateAssayObject(counts = rawdata[[i]]$`Antibody Capture`)
  processeddata[[i]]<-automated_TCR(paste("~/scratch60/2021analysis/TCRannotations",TCRanno[[i]], sep = "/"),processeddata[[i]])
  
  #perform the processing on each patient, this is variable based on which patient population is used 
  processeddata[[i]] <- PercentageFeatureSet(processeddata[[i]], pattern = "^MT-", col.name = "percent.mt")
  processeddata[[i]]<- subset(processeddata[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
  
  #perform hashtag filtering
  #for individual sample
  #hashtags<-grep("Hashtag", rownames(processeddata[[i]]@assays$ADT), value=FALSE)
  #hashtags<-hashtags[5:10]
  
  ## check rownames of ADT
  # rownames(processeddata[[1]]@assays$ADT@counts)
  
  ### make hashtag numbers in metadata.txt match ADT
  currenthashes<-as.character(metadata[[identifiers[[i]]]][[1]])
  currenthashes<-as.numeric(unlist(strsplit(currenthashes,",")))
  currenthashes<-paste("Hashtag", currenthashes, sep="")
  currenthashes<-paste(currenthashes,"-", sep="")
  
  currenthashes<-paste(currenthashes, collapse="|")
  currenthashes<-grep(currenthashes, rownames(processeddata[[i]]@assays$ADT@counts))
  
  ## This line will first collapse your list of hashtags into a regex expression with or. So its written as Hashtag-1|Hashtag-2 which is read as find something that matches Hashtag 1 or hashtag 2. We're searching the row names of the ADT to find them. 
  # currenthashes<-grep(paste(currenthashes, collapse="|"), rownames(processeddata[[i]]@assays$ADT@counts))
  
  
  # At this point, currenthashes is the numbers associated with your hashtags. 
  # This line will make the hashtag names so that we can search for them. 
  
  processeddata[[i]][["HTO"]]<-CreateAssayObject(processeddata[[i]]@assays$ADT@counts[currenthashes,])
  processeddata[[i]]<-subset(processeddata[[i]], subset = nCount_HTO>2)
  processeddata[[i]]<-NormalizeData(processeddata[[i]], assay = "HTO", normalization.method = "CLR")
  processeddata[[i]]<-HTODemux(processeddata[[i]],init=nrow(processeddata[[i]]@assays$HTO@counts))
  processeddata[[i]]<-subset(processeddata[[i]], subset = hash.ID!="Doublet" & hash.ID!="Negative")
  
  #need this so we can ensure that we can add on HIV later 
  processeddata[[i]]<-RenameCells(object = processeddata[[i]], add.cell.id = identifiers[[i]])
  #sctransform to start, binarize ADT library and then normalize it. Stopped regressing out percent.mt, it probably doesnt make too much of a difference since its already filtered 
  processeddata[[i]]<-NormalizeData(processeddata[[i]], assay="ADT", normalization.method = "CLR")
  
}

saveRDS(processeddata, "~/scratch60/2021analysis/checkpoint/210128checkpoint1.rds")

# check surface antibody names
# rownames(processeddata[[11]]@assays$ADT@counts)

#### no integration method

# SCT transform
processeddata<-lapply(processeddata, SCTransform) 

finalobj<-merge(processeddata[[1]], processeddata[[2]])
finalobj<-merge(finalobj,processeddata[[3]])
finalobj<-merge(finalobj,processeddata[[4]])
finalobj<-merge(finalobj,processeddata[[5]])
finalobj<-merge(finalobj,processeddata[[6]])
finalobj<-merge(finalobj,processeddata[[7]])
finalobj<-merge(finalobj,processeddata[[8]])
finalobj<-merge(finalobj,processeddata[[9]])
finalobj<-merge(finalobj,processeddata[[10]])
finalobj<-merge(finalobj,processeddata[[11]])
finalobj<-merge(finalobj,processeddata[[12]])

# finalobj<-merge(processeddata[[1]], processeddata[[2]],processeddata[[3]],processeddata[[4]],processeddata[[5]], processeddata[[6]], processeddata[[7]], processeddata[[8]],processeddata[[9]],processeddata[[10]],processeddata[[11]],processeddata[[12]])

saveRDS(finalobj, "~/scratch60/2021analysis/checkpoint/210129checkpoint2.rds")

