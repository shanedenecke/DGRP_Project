library(tcltk)
library(ggplot2)
library(dplyr)

file <- tk_choose.files(default="/home/sdenecke/Dropbox/PhD_Projects/DGRP_Paper/Results/Cyp6g1KO/Individual_Well_Data/",multi=F)
Data <- read.csv(file)
setwd(dirname(file))
dir.create("./GraphDir")
setwd("./GraphDir")

sub <- subset(Data,Time==60)
for (d in unique(sub$Dose)){       
            
        pdf(file="RAL517-KO_Roberto.pdf",width=10)
        gp <- ggplot(sub,aes(x=Genotype,y=Modified,fill=Genotype))
        gp <- gp+stat_summary(geom="errorbar",colour='black',size=1,width=.4,fun.data="mean_cl_normal")
        gp <- gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.4,colour='black',size=.8)
        gp <- gp+geom_dotplot(binaxis='y',stackdir="center",binwidth=.05)
        gp <- gp+scale_fill_manual(values=c("darkblue","lightblue")) 
        gp <- gp+ylab("Relative Motility\n")
        gp <- gp+xlab("\nGenotype")
        gp <- gp+ylim(c(0,1))
        #gp <- gp+scale_y_continuous(breaks=(c(0,.2,.4,.6,.8,1.0)))
        gp <- gp+ggtitle(expression(paste(italic("Cyp6g1"),"Knockout")))
        gp <- gp+theme_bw(base_size=16,base_family="serif") 
        gp <- gp+theme(title=element_text(size=20,face="bold"),
                       axis.text.x=element_text(angle=30,hjust=1),strip.text=element_text(size=18,face="bold"),
                       strip.background=element_rect("grey95"),panel.grid.minor=element_blank(),
                       panel.grid.major=element_blank(),legend.position="none",title=element_text(size=24)) 
        
        print(gp)
        dev.off()
}


#################### Combined

library(tcltk)
library(ggplot2)
library(dplyr)
file <- tk_choose.files(default="/home/sdenecke/Dropbox/PhD_Projects/DGRP_Paper/Results/Cyp6g1KO/Individual_Well_Data/",multi=F)
Data <- read.csv(file)
setwd(dirname(file))
dir.create("./GraphDir")
setwd("./GraphDir")

sub <- subset(Data,Time==60)

sub$Genotype <- factor(sub$Genotype,levels=c("Canton-S","Cs-6g1KO","Wxac","Wx-6g1KO","RAL517","RAL-6g1KO"))
sub$Set <- factor(sub$Set,levels=c("Canton-S (5ppm)","Wxac (5ppm)","RAL517 (25ppm)","RAL517 (100ppm)"))
pdf("Cyp6g1KO_Comparison.pdf")
gp <- ggplot(sub,aes(x=Genotype,y=Modified,fill=Genotype))
gp <- gp+facet_wrap(facets=~Set,scales="free_x",nrow=2,ncol=2)
gp <- gp+stat_summary(geom="errorbar",colour='black',size=1,width=.4,fun.data="mean_cl_normal")
gp <- gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.4,colour='black',size=.8)
gp <- gp+geom_dotplot(binaxis='y',stackdir="center",binwidth=.05)
gp <- gp+scale_fill_manual(values=c("darkblue","lightblue","darkred","pink","darkgreen","lightgreen"))
gp <- gp+ylab("Relative Motility\n")
gp <- gp+xlab("\nGenotype")
gp <- gp+scale_y_continuous(breaks=(c(0,.2,.4,.6,.8,1.0)))
gp <- gp+ggtitle(expression(paste(italic("Cyp6g1"),"Knockouts")))
gp <- gp+theme_bw(base_size=16,base_family="serif") 
gp <- gp+theme(title=element_text(size=20,face="bold"),
               axis.text.x=element_text(angle=30,hjust=1),strip.text=element_text(size=18,face="bold"),
               strip.background=element_rect("grey95"),panel.grid.minor=element_blank(),
               panel.grid.major=element_blank(),legend.position="none",title=element_text(size=24)) 

print(gp)
dev.off()
