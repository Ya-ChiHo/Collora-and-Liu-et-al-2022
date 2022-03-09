library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(ggpubr)


setwd("~/project/JointAnalyses/")

meta<-readRDS("Objects/metadata with both clone timepoint quantifications")

HIVcells<-readRDS("Objects/20210601CellsincloneswithHIV.rds")
meta$HIV <- rownames(meta)%in%HIVcells
pubfig169p<-function(filename){pdf(file  = paste(filename,".pdf", sep=""),height = 9, width = 16)}

#goal is to figure out if the clones infected are bigger or smaller than expected. 
#approach is to take the total metadata of captured cells, sample it with a n equivalent to the number of HIV-1 positive clones in that group
#then wilcox test and find how often the clones are "larger" than expected (one sided wilcoxon) 
#repeat 10k times 
#must be without replacement to be a proper permutation test 

#doing each antigen independently, recording p value and the mean of each sample for plotting
#must do both viermic and suppressed, but viremic is the more relevant since that's where the integration is occurring

Clonal_permutation<-function(data, antigen, nperm=10000, stage="Viremic", freq_col="bulkfreq_1"){
  datacur<-data[data$antigen2==antigen & data$infection2==stage,]
  dataHIV<-datacur[datacur$HIV==TRUE,]
  datacur<-datacur[datacur$bulk>1,]
  nHIV<-nrow(dataHIV)
  dataHIV<- group_by(dataHIV, junction2)%>%slice_max(n=1, !!freq_col,with_ties = FALSE)
  HIV_dist<-dataHIV[[freq_col]]
  HIVmean<-mean(HIV_dist)
  datacur<-datacur[,c("junction2","bulkfreq_1")]
  
  res_mean<-c()
  res_p<-c()
  
  for (i in 1:10000){
    #for each of these, step one is going to be to take a sample
    cur<-datacur[sample(1:nrow(datacur), nHIV),]
    #then collapse that sample into clones (ie remove secondary junctions)
    cur<-group_by(cur, junction2)%>%slice_max(n=1, !!freq_col,with_ties = FALSE)
    #then perform the test, add the mean
    res_mean<-c(res_mean, mean(cur[[freq_col]]))
    res_p<-c(res_p,wilcox.test(HIV_dist, cur[[freq_col]], alternative = "greater")$p.value)
  }
  return(list(HIV_dist, HIVmean,res_mean, res_p))
}

#completing it all for viremic

#within CMV and HIV, the p value are not significant during viremia 
resultsCMV<-Clonal_permutation(meta, "CMV")
resultsHIV<-Clonal_permutation(meta, "HIV")

#within resting and memory, the p value are highly significant (100% of theses performed are significant)
resultsMemory<-Clonal_permutation(meta, "Memory")
resultsresting<-Clonal_permutation(meta, "resting")

#completing it all for supp
resultsCMVS<-Clonal_permutation(meta, "CMV", stage="Suppressed")
resultsHIVS<-Clonal_permutation(meta, "HIV", stage="Suppressed")
resultsMemoryS<-Clonal_permutation(meta, "Memory", stage="Suppressed")
resultsrestingS<-Clonal_permutation(meta, "resting", stage="Suppressed")

current<-as.data.frame(resultsCMV[[3]])
colnames(current)<-"means"
current$ID<-"perms"
current<-rbind(current,c(resultsCMV[[2]], "Obs"))
current$means<-as.numeric(current$means)

p1<-ggplot(current, aes(x=means, fill=ID))+geom_histogram(aes(y=0.5*..density..), alpha=0.5,position='identity', bins = 100)+scale_x_log10(limits = c(5e-5,5e-3)) +theme_classic()

current<-as.data.frame(resultsHIV[[3]])
colnames(current)<-"means"
current$ID<-"perms"
current<-rbind(current,c(resultsHIV[[2]], "Obs"))
current$means<-as.numeric(current$means)

p2<-ggplot(current, aes(x=means, fill=ID))+geom_histogram(aes(y=0.5*..density..), alpha=0.5,position='identity', bins = 100)+scale_x_log10(limits = c(5e-5,5e-3)) +theme_classic()

current<-as.data.frame(resultsMemory[[3]])
colnames(current)<-"means"
current$ID<-"perms"
current<-rbind(current,c(resultsMemory[[2]], "Obs"))
current$means<-as.numeric(current$means)

p3<-ggplot(current, aes(x=means, fill=ID))+geom_histogram(aes(y=0.5*..density..), alpha=0.5,position='identity', bins = 100)+scale_x_log10(limits = c(5e-5,5e-3)) +theme_classic()

current<-as.data.frame(resultsresting[[3]])
colnames(current)<-"means"
current$ID<-"perms"
current<-rbind(current,c(resultsresting[[2]], "Obs"))
current$means<-as.numeric(current$means)

p4<-ggplot(current, aes(x=means, fill=ID))+geom_histogram(aes(y=0.5*..density..), alpha=0.5,position='identity', bins = 100)+scale_x_log10(limits = c(5e-5,5e-3)) +theme_classic()
                                                      

current<-as.data.frame(resultsCMVS[[3]])
colnames(current)<-"means"
current$ID<-"perms"
current<-rbind(current,c(resultsCMVS[[2]], "Obs"))
current$means<-as.numeric(current$means)

p5<-ggplot(current, aes(x=means, fill=ID))+geom_histogram(aes(y=0.5*..density..), alpha=0.5,position='identity', bins = 100)+scale_x_log10(limits = c(5e-5,5e-3)) +theme_classic()

current<-as.data.frame(resultsHIVS[[3]])
colnames(current)<-"means"
current$ID<-"perms"
current<-rbind(current,c(resultsHIVS[[2]], "Obs"))
current$means<-as.numeric(current$means)

p6<-ggplot(current, aes(x=means, fill=ID))+geom_histogram(aes(y=0.5*..density..), alpha=0.5,position='identity', bins = 100)+scale_x_log10(limits = c(5e-5,5e-3)) +theme_classic()

current<-as.data.frame(resultsMemoryS[[3]])
colnames(current)<-"means"
current$ID<-"perms"
current<-rbind(current,c(resultsMemoryS[[2]], "Obs"))
current$means<-as.numeric(current$means)

p7<-ggplot(current, aes(x=means, fill=ID))+geom_histogram(aes(y=0.5*..density..), alpha=0.5,position='identity', bins = 100)+scale_x_log10(limits = c(5e-5,5e-3)) +theme_classic()

current<-as.data.frame(resultsrestingS[[3]])
colnames(current)<-"means"
current$ID<-"perms"
current<-rbind(current,c(resultsrestingS[[2]], "Obs"))
current$means<-as.numeric(current$means)

p8<-ggplot(current, aes(x=means, fill=ID))+geom_histogram(aes(y=0.5*..density..), alpha=0.5,position='identity', bins = 100)+scale_x_log10(limits = c(5e-5,5e-3)) +theme_classic()

pubfig169p("Figures/supplemental/CloneSizepermfixedX")
ggarrange(p1,p2,p3,p4, nrow=1)/ggarrange(p5,p6,p7,p8, nrow=1)
dev.off()
