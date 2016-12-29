library(reshape2)
library(ggplot2)
setwd("~/Documents")
raw.data <- read.table("Fly_Atlas_Raw.txt",sep="\t",header=T)

annotations <- read.csv("Annotations.csv")
colnames(annotations)[1]="Oligo"
merged.raw <- merge(raw.data,annotations,by="Oligo")

gene.list <- as.character(read.csv("Top_Annotations_Summary_No_Intergenic.csv")$Nearest.Gene)
relevant.columns <- merged.raw[which(merged.raw$Gene.Symbol%in%gene.list),c(150,97,2)] 
relevant.columns <- unique(relevant.columns)
length(which(relevant.columns$Feeded.larvae.central.nerve.system.vs.whole.flies....T.Test_Change.Direction=="None"))
