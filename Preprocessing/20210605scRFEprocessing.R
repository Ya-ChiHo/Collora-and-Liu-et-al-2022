library(dplyr)
setwd("scratch60/scRFEresting/")
#goal of this script is to process scRFE replicates, first read in files 
todo<-list.files()
results<-list()
for(i in todo){
  results[[i]]<-read.csv(i)
}

#helper function that takes the genes in the top 100, divides by the total replicates to get prop in top 100
top100<-function(results){
  genes<-c()
  for(i in results){
    cur<-i
    genes<-c(genes, cur$resting[1:100])
  }
  genes<-genes[!is.na(genes)]
  genes<-table(genes)/length(results)
  return(genes)
}


#running this function
top100res<-top100(results)
#writing results
write.csv(top100res, "~/project/JointAnalyses/sklearn/20210605top100scRFE.csv")
