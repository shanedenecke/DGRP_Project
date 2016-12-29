library(dplyr)
library(ggplot2)
setwd("/home/sdenecke/Dropbox/PhD_Projects/DGRP_Paper/Results/Post_GWAS/Transcriptome_Associations")
dir.create(format(Sys.Date(),format="%d.%b.%y"))
phen <- read.csv("Mean_Combined.Doses.csv")
allele <- read.csv("CYP6G1_GENO.csv",header=T)
Doses <- c("25ppm","100ppm")
f.tran <- read.csv("F185TRANS_ANNOT.csv")
setwd(paste("./",format(Sys.Date(),format="%d.%b.%y"),sep=""))


## Edit Column Names

## Find common genotypes between 
## transcriptome dataset and phenotype dataset
trans.gen <- colnames(f.tran[-1:-3])
phen.gen <- as.character(phen$Genotype)
common.gen<- intersect(trans.gen,phen.gen)

## Add allele column for each genotype used 
allele.sub <- allele %>% filter(Genotype%in%common.gen)

## Find transcriptional data for each gene required
gene.list <- c("Cyp6g1","Cyp6g2","GstD1","GstD2")
rmv <- c("FBID","LOCATION")
f.tran.gen <- f.tran %>% select(one_of(c("SYMBOL",common.gen))) %>%
      filter(SYMBOL %in% gene.list) %>% 
      gather(Genotype, val, 2:length(c("SYMBOL",common.gen))) %>%  
      spread(SYMBOL, val) %>% arrange(Genotype)

## Modify phenotypic dataset
phen.sub <- phen %>% filter(Genotype %in% common.gen) %>% 
            arrange(Dose,Genotype)

## Create Summary plot data
summary.data <- full_join(f.tran.gen,phen.sub,by="Genotype") %>% 
            full_join(allele.sub,by="Genotype")




d="25ppm"
plot.data <- filter(summary.data,Dose==d)
pdf(file=paste(d,"Cyp6g1","pdf",sep="."))
gp <- ggplot(plot.data,aes(x=Cyp6g1,y=Mean_Wiggle))
gp <- gp+geom_point(size=4,aes(shape=Allele,colour=Allele))
gp <- gp+scale_shape_identity()
gp <- gp+geom_smooth(method="lm")
gp <- gp+ggtitle(d)
gp <- gp+ylab("RMR Value\n")
gp <- gp+xlab("\nCyp6g1 Expression")
gp <- gp+scale_colour_manual(values=c("red2","blue","green4","black"))
gp <- gp+theme_bw(base_size=16,base_family="serif")
gp <- gp+theme(panel.background=element_rect(fill=NA),
            text=element_text(face="bold"),
            panel.border=element_rect(colour="black",fill=NA),
            panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())
print(gp) 
dev.off()


## Post Plot analysis

cor.test(plot.data$Mean_Wiggle,plot.data$Cyp6g1,method="pearson")

sum <- summary(with(plot.data,lm(Mean_Wiggle~Cyp6g1)))
plot.data.resid <- mutate(plot.data,Cyp6g1.residual=sum$residual)
gen.sin.ral <- as.numeric(sapply(plot.data.resid$Genotype,function(x){gsub("RAL","",x)}))
write.table(cbind(gen.sin.ral,as.numeric(plot.data.resid$Cyp6g1.residual)),
      file="Cyp6g1_Residual_submit.csv",row.names=F,col.names=F,sep=",")
      
