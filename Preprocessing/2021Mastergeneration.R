library(Seurat)
library(dplyr)
####unstimulated####
#load in integrated object 
samples<-readRDS("FASTMNN.rds")
#add bulk TCR quantification
BulkMeta<-readRDS("BulkQuantificationForRunxia.rds")
JackMeta<-samples@meta.data
#first by UID
BulkUIDSplit<-split(BulkMeta, BulkMeta$ID2)
JackUIDSplit<-split(JackMeta, RunxiaMeta$UID)

todo<-names(JackUIDSplit)

for (i in todo){
  #first generate a lookup table for each sample
  lookuptable<-BulkUIDSplit[[i]]
  colnames(lookuptable)<-paste("bulk",colnames(lookuptable), sep="_")
  JackUIDSplit[[i]]<-merge(JackUIDSplit[[i]],lookuptable, by.y = "bulk_JUNCTION",by.x = "TRB_junction",all.x=TRUE)
}

JackUIDSplit<-bind_rows(JackUIDSplit)

# add rownames
rownames(JackUIDSplit)<-JackUIDSplit$cell
samples<-AddMetaData(samples, JackUIDSplit)

#binarize CITE data

#get antibody names
antibodies<-rownames(samples@assays$ADT@data)
#need to make a list of isotypes, so antibody one has isotype X which is antibody 23 in the list
isotypes<-c(antibodies[22], antibodies[22], antibodies[22], antibodies[24],antibodies[23],
            antibodies[22], antibodies[22], antibodies[23], antibodies[22],antibodies[22],
            antibodies[24], antibodies[22], antibodies[22], antibodies[22],antibodies[22],
            antibodies[22], antibodies[23], antibodies[23], antibodies[22],antibodies[22],
            antibodies[23], antibodies[22], antibodies[23], antibodies[24],antibodies[25],
            antibodies[22], antibodies[22], antibodies[22], antibodies[22],antibodies[22],
            antibodies[22], antibodies[22], antibodies[22], antibodies[22],antibodies[22],
            antibodies[22],antibodies[22])

#done by patient 
splitsamples<-SplitObject(samples, split.by="orig.ident")
metadata<-list()
#for each antibody
for (z in 1:length(antibodies)){
  #for each patient 
  for (i in 1:length(splitsamples)){
    #get the antibody cutoff (isotype 90th percentile)
    cutoff<-quantile(splitsamples[[i]]@assays$ADT@data[isotypes[[i]],],c(0.90))
    #if the value is greater than that cutoff, it's 1
    res<-case_when(splitsamples[[i]]@assays$ADT@data[antibodies[z],]>=cutoff~1, T~0)
    #name that one
    names(res)<-colnames(splitsamples[[i]])
    #append it to the rest
    metadata[[antibodies[z]]]<-c(metadata[[antibodies[z]]], res)
  }}
#for each of these you then add the metadata 
for (i in 1:length(metadata)){
  samples<-AddMetaData(samples, metadata[[i]], col.name = gsub("TotalSeqC","binary",names(metadata)[[i]]))
}

#here's how you then make classificaitons

samples$memory<-case_when(samples$CD45RA.binary==1 & samples$CCR7.binary==1 ~ "Naive", 
                          samples$CD45RA.binary==0 & samples$CCR7.binary==1 ~ "Central Memory",
                          samples$CD45RA.binary==1 & samples$CCR7.binary==0 ~ "Effector",
                          samples$CD45RA.binary==0 & samples$CCR7.binary==0 ~ "Effector Memory")
samples$exhaustion<-case_when(samples$PD.1.binary==1 & samples$TIGIT.binary==1 ~ "Double Positive", 
                              samples$PD.1.binary==0 & samples$TIGIT.binary==1 ~ "TIGIT Positive",
                              samples$PD.1.binary==1 & samples$TIGIT.binary==0 ~ "PD1 Positive",
                              samples$PD.1.binary==0 & samples$TIGIT.binary==0 ~ "Double Negative")
samples$activation<-case_when(samples$HLA.DR.binary==1 & samples$CD25.binary==1 ~ "Double Positive", 
                              samples$HLA.DR.binary==0 & samples$CD25.binary==1 ~ "CD25 positive",
                              samples$HLA.DR.binary==1 & samples$CD25.binary==0 ~ "HLA DR Positive",
                              samples$HLA.DR.binary==0 & samples$CD25.binary==0 ~ "Double Negative")


#save as the new "master file" for downstream analysis 

saveRDS(samples, "20210519UnstimulatedMaster.rds")

####stimulated#####

#load in integrated object 
samples<-readRDS("210130checkpoint3-integrationfastMNN.rds")
#eject SARSdata
new<-subset(samples, subset=antigen3!="hdSARSmemory")

#add bulk TCR quantification
#making the participant and time point naming consistent 

RunxiaMeta$PT<-gsub("[a-zA-Z]{1,2}","",RunxiaMeta$participant.ID)
RunxiaMeta$PT<-gsub("2004","M3",RunxiaMeta$PT)
RunxiaMeta$PT<-gsub("2022","M2",RunxiaMeta$PT)

RunxiaMeta$TP<-case_when(RunxiaMeta$stage=="HDTime1" |RunxiaMeta$stage=="Viremic"~1, T~2 )

RunxiaMeta$UID<-paste(RunxiaMeta$PT, RunxiaMeta$TP, sep="_")
table(RunxiaMeta$UID)


# For bulk data, always use time point + individual

BulkMeta<-readRDS("BulkQuantificationForRunxia.rds")
RunxiaMeta<-samples@meta.data
#first by UID
BulkUIDSplit<-split(BulkMeta, BulkMeta$ID2)
RunxiaUIDSplit<-split(RunxiaMeta, RunxiaMeta$UID)

todo<-names(RunxiaUIDSplit)
#check that all are present in both; they are not matching
# bulk doesn't have M2_2; but it has M1_1 and M1_2 which are not in RunxiaMeta
todo==names(RunxiaUIDSplit)

for (i in todo){
  #first generate a lookup table for each sample
  lookuptable<-BulkUIDSplit[[i]]
  colnames(lookuptable)<-paste("bulk",colnames(lookuptable), sep="_")
  RunxiaUIDSplit[[i]]<-merge(RunxiaUIDSplit[[i]],lookuptable, by.y = "bulk_JUNCTION",by.x = "TRB_junction",all.x=TRUE)
}

# M2_2 is missing in BulkUIDSplit. So any samples after that (M2_2, M3_1 and M3_2) are not added the bulk information.

colnames(RunxiaUIDSplit$`236_1`)
lapply(RunxiaUIDSplit, colnames)

# manually add M3_1 (the 16th in BulkUIDSplit and the 15th in RunxiaUIDSplit)
lookuptable<-BulkUIDSplit[[16]]
colnames(lookuptable)<-paste("bulk",colnames(lookuptable), sep="_")
RunxiaUIDSplit[[15]]<-merge(RunxiaUIDSplit[[15]],lookuptable, by.y = "bulk_JUNCTION",by.x = "TRB_junction",all.x=TRUE)

# manually add M3_2 (the 17th in BulkUIDSplit and the 16th in RunxiaUIDSplit)
lookuptable<-BulkUIDSplit[[17]]
colnames(lookuptable)<-paste("bulk",colnames(lookuptable), sep="_")
RunxiaUIDSplit[[16]]<-merge(RunxiaUIDSplit[[16]],lookuptable, by.y = "bulk_JUNCTION",by.x = "TRB_junction",all.x=TRUE)

# bind function when the number of columns are not equal
RunxiaUIDSplit<-bind_rows(RunxiaUIDSplit)

# check if added bulk information for M2_2 is all NA
M2_2<-subset(RunxiaUIDSplit, RunxiaUIDSplit$UID=='M2_2')
table(is.na(M2_2$bulk_ID2))

# add rownames
rownames(RunxiaUIDSplit)<-RunxiaUIDSplit$cell
samples<-AddMetaData(samples, RunxiaUIDSplit)
#binarize CITE data

#get antibody names
antibodies<-rownames(samples@assays$ADT@data)
#need to make a list of isotypes, so antibody one has isotype X which is antibody 23 in the list
isotypes<-c(antibodies[22], antibodies[22], antibodies[22], antibodies[24],antibodies[23],
            antibodies[22], antibodies[22], antibodies[23], antibodies[22],antibodies[22],
            antibodies[24], antibodies[22], antibodies[22], antibodies[22],antibodies[22],
            antibodies[22], antibodies[23], antibodies[23], antibodies[22],antibodies[22],
            antibodies[23], antibodies[22], antibodies[23], antibodies[24],antibodies[25],
            antibodies[22], antibodies[22], antibodies[22], antibodies[22],antibodies[22],
            antibodies[22], antibodies[22], antibodies[22], antibodies[22],antibodies[22],
            antibodies[22],antibodies[22])

#done by patient 
splitsamples<-SplitObject(samples, split.by="orig.ident")
metadata<-list()
#for each antibody
for (z in 1:length(antibodies)){
  #for each patient 
  for (i in 1:length(splitsamples)){
    #get the antibody cutoff (isotype 90th percentile)
    cutoff<-quantile(splitsamples[[i]]@assays$ADT@data[isotypes[[i]],],c(0.90))
    #if the value is greater than that cutoff, it's 1
    res<-case_when(splitsamples[[i]]@assays$ADT@data[antibodies[z],]>=cutoff~1, T~0)
    #name that one
    names(res)<-colnames(splitsamples[[i]])
    #append it to the rest
    metadata[[antibodies[z]]]<-c(metadata[[antibodies[z]]], res)
  }}
#for each of these you then add the metadata 
for (i in 1:length(metadata)){
  samples<-AddMetaData(samples, metadata[[i]], col.name = gsub("TotalSeqC","binary",names(metadata)[[i]]))
}

#here's how you then make classificaitons

samples$memory<-case_when(samples$CD45RA.binary==1 & samples$CCR7.binary==1 ~ "Naive", 
                          samples$CD45RA.binary==0 & samples$CCR7.binary==1 ~ "Central Memory",
                          samples$CD45RA.binary==1 & samples$CCR7.binary==0 ~ "Effector",
                          samples$CD45RA.binary==0 & samples$CCR7.binary==0 ~ "Effector Memory")
samples$exhaustion<-case_when(samples$PD.1.binary==1 & samples$TIGIT.binary==1 ~ "Double Positive", 
                              samples$PD.1.binary==0 & samples$TIGIT.binary==1 ~ "TIGIT Positive",
                              samples$PD.1.binary==1 & samples$TIGIT.binary==0 ~ "PD1 Positive",
                              samples$PD.1.binary==0 & samples$TIGIT.binary==0 ~ "Double Negative")
samples$activation<-case_when(samples$HLA.DR.binary==1 & samples$CD25.binary==1 ~ "Double Positive", 
                              samples$HLA.DR.binary==0 & samples$CD25.binary==1 ~ "CD25 positive",
                              samples$HLA.DR.binary==1 & samples$CD25.binary==0 ~ "HLA DR Positive",
                              samples$HLA.DR.binary==0 & samples$CD25.binary==0 ~ "Double Negative")


samples$stage2<-case_when(samples$stage=="HDTime1" | samples$stage=="HDTime2" ~ "Uninfected", 
                          samples$stage=="Viremic" ~ "Viremic",
                          samples$stage=="Suppressed" ~ "Suppressed")

#save as the new "master file" for downstream analysis 

saveRDS(samples, "210623masterstimulated.rds")
