library(dplyr)
library(ggplot2)
setwd("~/Documents/HoLab/JoinedAnalyses/Enrichr/")
scRFE<-read.table("GO_scRFEup.txt", sep="\t", header=1)
scRFE$value<-(-log10(scRFE$Adjusted.P.value))
scRFE<-scRFE[order(scRFE$value, decreasing = TRUE),]
scRFE<-scRFE[1:10,]
scRFE$Term<-factor(scRFE$Term, levels = rev(scRFE$Term))
pdf("scRFEGOupEnrichr.pdf")
ggplot(scRFE, aes(y=Term, x=value))+geom_bar(stat="identity")+theme_classic()
dev.off()
scRFE<-read.table("GO_HIVup.txt", sep="\t", header=1)
scRFE$value<-(-log10(scRFE$Adjusted.P.value))
scRFE<-scRFE[order(scRFE$value, decreasing = TRUE),]
scRFE<-scRFE[1:10,]
scRFE$Term<-factor(scRFE$Term, levels = rev(scRFE$Term))
pdf("HIVUnstimulatedGOupEnrichr.pdf")
ggplot(scRFE, aes(y=Term, x=value))+geom_bar(stat="identity")+theme_classic()
dev.off()
                                             