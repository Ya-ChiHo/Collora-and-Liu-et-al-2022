#reliant on the utilities R script, generates all of the module correlation maps, including those shown in figure 1

setwd("~/project/JointAnalyses/")

todo<-list.files("Objects/cor/")

todo<-grep(".{7,}_unstimulated.rds", todo, value=TRUE)
todo<-grep(" ", todo, invert = TRUE, value = TRUE)

modules<-readRDS("Objects/20210803UnstimulatedConsensusmods.rds")
results<-list()
for(i in 1:length(todo)){
  results[[i]]<-readRDS(paste("Objects/cor",todo[[i]], sep="/"))
}

names(results)<-todo

for (i in 1:length(modules)){splitwrap(results, modules[[i]], prefix=paste("Figures/splitcorplots/", "module",i, sep=""))}
