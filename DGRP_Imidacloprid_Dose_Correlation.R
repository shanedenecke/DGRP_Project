library(dplyr)
library(ggplot2)
library(tidyr)
setwd("/home/sdenecke/Dropbox/PhD_Projects/DGRP_Paper/Results/Pre_GWAS")
data <- read.csv("Completed_Data_Frame.csv")

data.sub <- filter(data,Time==60)

data.sum <- data %>% filter(Time==60) %>% group_by(Genotype,Dose) %>% summarise(RMR=mean(Modified))
data.long <- spread(data=data.sum,key=Dose,value=RMR)


pdf(file="25.100ppm_Corr_Plot.pdf")
gp <- ggplot(data.long,aes(x=`25ppm`,y=`100ppm`))
gp <- gp+geom_point(size=2)
gp <- gp+geom_smooth(method="lm")
gp <- gp+ggtitle("Correlation of 25ppm and 100ppm \nImidacloprid Response")
gp=gp+ylab("100ppm Response \n")
gp=gp+xlab("\n25ppm Response")
gp=gp+theme(panel.background=element_rect(fill=NA),
            text=element_text(size=16,face="bold",family="serif"),
            panel.border=element_rect(colour="black",fill=NA),
            panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())
print(gp) 
dev.off()
