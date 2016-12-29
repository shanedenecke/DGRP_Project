## Phenotype Correlation

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpmisc)
setwd("/home/sdenecke/Dropbox/PhD_Projects/DGRP_Paper/Results/Post_GWAS/Insecticide_Response_Correlations/")
data <- read.csv("Summary_Clean.csv")

## Data manipulations to rearrange the data frame
data.wide <- spread(data,key=Phenotype,value=Value)
data.sub <- data.wide[complete.cases(data.wide),]
data.clean <- data.sub %>% 
      gather(key=Other_Insecticide_Type,Other_Insecticide_Value,Azinphos_Larvae:Imi_Adult) %>% 
      gather(key=Wiggle_Dose,value=RMR,Wiggle_100:Wiggle_25)
     

pdf(file="Phenotypic correlations plot.pdf")
gp <- ggplot(data.clean,aes(x=Other_Insecticide_Value,y=RMR))
gp <- gp+facet_grid(Wiggle_Dose~Other_Insecticide_Type,scales="free")
gp <- gp+geom_point(size=2)
gp <- gp+geom_smooth(method="lm")
#gp <- gp+stat_poly_eq(formula = Other_Insecticide_Value~RMR, 
#                    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
#                   parse = TRUE) 
#gp <- gp+ylim(c(0,1.5))
gp <- gp+theme(panel.background=element_rect(fill=NA),
            text=element_text(size=16,face="bold",family="serif"),
            panel.border=element_rect(colour="black",fill=NA),
            panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())
print(gp) 
dev.off()

Emat <- matrix(nrow=length(unique(data.clean$Wiggle_Dose)),ncol=length(unique(data.clean$Other_Insecticide_Type)))
rownames(Emat) <- unique(data.clean$Wiggle_Dose)
colnames(Emat) <- unique(data.clean$Other_Insecticide_Type)
for (I in unique(data.clean$Wiggle_Dose)){
      for (O in unique(data.clean$Other_Insecticide_Type)){
            sub <- filter(data.clean,Wiggle_Dose==I & Other_Insecticide_Type==O)
            Emat[I,O] <- cor.test(sub$Other_Insecticide_Value,sub$RMR,method="pearson")$p.val
      }}

write.csv(Emat,"Insecticide_Correlation_Summary.csv")
